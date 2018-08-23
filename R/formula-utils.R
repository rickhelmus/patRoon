splitFormulaToList <- function(formula)
{
    # split string in pairs of elements+element counts (and optionally isotopic info), e.g.: { "C30", "^13C2" }
    spltform <- unlist(regmatches(formula, gregexpr("(\\^[[:digit:]]+)?[[:upper:]]{1}[[:lower:]]?[[:digit:]]*", formula)))

    # add '1' to pairs without a count, e.g. "Br" --> "Br1"
    reglist <- regexpr("[[:digit:]]+$", spltform)
    spltform[reglist == -1] <- paste0(spltform[reglist == -1], "1")

    # extract all element counts
    ret <- as.numeric(unlist(regmatches(spltform, gregexpr("[[:digit:]]+$", spltform))))

    # extract all elements (ie remove all counts)
    names(ret) <- gsub("[[:digit:]]+$", "", spltform)

    if (any(duplicated(names(ret))))
    {
        unel <- unique(names(ret))
        sortedr <- numeric()
        sortedr[unel] <- sapply(unel, function(e) sum(ret[names(ret) == e]))
        ret <- sortedr
    }

    return(ret)
}

getElements <- function(formula, elements)
{
    ret <- sapply(formula, function(f)
    {
        fl <- splitFormulaToList(f)
        df <- as.data.frame(t(as.data.frame(fl)))
        df[setdiff(elements, names(fl))] <- 0
        return(df[, elements, drop=F])
    }, simplify = F, USE.NAMES = F)

    ret <- as.data.frame(do.call(rbind, ret))
    rownames(ret) <- NULL
    return(ret)
}

formulaListToString <- function(formlist)
{
    el <- names(formlist)

    # make element - element count pairs (don't put element counts between 0-1)
    l <- lapply(el, function(e) { if (formlist[e] < 0 || formlist[e] > 1) paste0(e, formlist[e]) else e })
    return(do.call(paste0, l))
}

subtractFormula <- function(formula1, formula2)
{
    if (is.na(formula1) || is.na(formula2))
        return(NA)

    if (formula1 == formula2)
        return("")

    fl1 <- splitFormulaToList(formula1)
    fl2 <- splitFormulaToList(formula2)

    newfl <- fl1
    newfl[setdiff(names(fl2), names(newfl))] <- 0 # add missing elements with zero count
    ind <- match(names(fl2), names(newfl))
    newfl[ind] <- newfl[ind] - fl2
    newfl <- newfl[newfl != 0] # filter out zero elements

    if (length(newfl) == 0)
        return("")

    return(formulaListToString(newfl))
}

addFormula <- function(formula1, formula2)
{
    if (is.na(formula1) || is.na(formula2))
        return(NA)

    fl1 <- splitFormulaToList(formula1)
    fl2 <- splitFormulaToList(formula2)

    newfl <- fl1
    newfl[setdiff(names(fl2), names(newfl))] <- 0 # add missing elements with zero count
    ind <- match(names(fl2), names(newfl))
    newfl[ind] <- newfl[ind] + fl2

    return(formulaListToString(newfl))
}

calculateIonFormula <- function(formula, adduct)
{
    if (grepl("+H", adduct, fixed = TRUE))
        sapply(formula, addFormula, formula2 = "H", USE.NAMES = FALSE)
    else if (grepl("+Na", adduct, fixed = TRUE))
        sapply(formula, addFormula, formula2 = "Na", USE.NAMES = FALSE)
    else if (grepl("+K", adduct, fixed = TRUE))
        sapply(formula, addFormula, formula2 = "K", USE.NAMES = FALSE)
    else if (grepl("-H", adduct, fixed = TRUE))
        sapply(formula, subtractFormula, formula2 = "H", USE.NAMES = FALSE)
}

sortFormula <- function(formula)
{
    fl <- splitFormulaToList(formula)
    el <- names(fl)
    hasC <- "C" %in% el; hasH <- "H" %in% el
    el <- el[!el %in% c("C", "H")]
    el <- el[order(el)]
    if (hasH)
        el <- c("H", el)
    if (hasC)
        el <- c("C", el)
    return(formulaListToString(fl[el]))
}

formulaScoringColumns <- function() c("score", "MS_match", "treeScore", "isoScore",
                                      "frag_score", "MSMS_match", "comb_match")

consensusForFormulaList <- function(formList, fGroups, formThreshold, maxFormulas,
                                    maxFragFormulas, minIntensity, maxIntensity, minPreferredFormulas,
                                    minPreferredIntensity)
{
    printf("Generating consensus for formula lists (%s)... ", algorithm(formList))

    ftind <- groupFeatIndex(fGroups)
    gInfo <- groupInfo(fGroups)
    anaInfo <- analysisInfo(fGroups)
    gTable <- groups(fGroups)
    forms <- formulaTable(formList)

    # Add analysis & group names
    for (a in names(forms))
    {
        if (length(forms[[a]]) > 0)
        {
            for (g in names(forms[[a]]))
            {
                if (nrow(forms[[a]][[g]]) > 0)
                {
                    # BUG: need to copy() otherwise table will not be updated (because it's nested?)
                    forms[[a]][[g]] <- copy(forms[[a]][[g]])
                    forms[[a]][[g]][, c("analysis", "group") := .(a, g)]
                }
            }
        }
    }

    # merge all together
    formTable <- rbindlist(lapply(forms, rbindlist, fill = TRUE), fill = TRUE)
    haveMSMS <- "frag_formula" %in% colnames(formTable)

    gTable <- groups(fGroups)
    gInfo <- groupInfo(fGroups)

    # Filter irrelevant feature groups
    formTable <- formTable[group %in% colnames(gTable)]

    # add some general info and put to the front
    formTable[, c("ret", "mz") := .(gInfo[group, "rts"], gInfo[group, "mzs"])]
    setcolorder(formTable, c((ncol(formTable)-1):ncol(formTable), 1:(ncol(formTable)-2)))

    if (nrow(formTable) == 0)
        warning("No (relevant) formulas!")
    else
    {
        # (temporarily) add intensities
        formTable[, intensity := gTable[[group]][match(analysis, anaInfo$analysis)], by = .(group, analysis)]

        # filter formulas above/below intensity thresholds
        if (!is.null(minIntensity))
            formTable <- formTable[intensity >= minIntensity]
        if (!is.null(maxIntensity))
            formTable <- formTable[intensity <= maxIntensity]
        
        if (nrow(formTable) == 0)
            warning("Filtered all formulas!")
        else
        {
            # check if we can remove some more to stay in optimal intensity range
            if (!is.null(minPreferredIntensity))
            {
                formTable[, prefAnaCount := length(unique(.SD[intensity >= minPreferredIntensity]$analysis)), by = group]
                formTable <- formTable[prefAnaCount < minPreferredFormulas | intensity >= minPreferredIntensity]
                formTable[, prefAnaCount := NULL]
            }
            
            formTable[, min_intensity := min(intensity), by = .(group, byMSMS)]
            formTable[, max_intensity := max(intensity), by = .(group, byMSMS)]
            
            formTable[, ana_min_intensity := .SD[which.min(intensity), analysis], by = .(group, byMSMS)]
            formTable[, ana_max_intensity := .SD[which.max(intensity), analysis], by = .(group, byMSMS)]
            
            # number of analyses searched for formulas per group
            formTable[, anaCount := length(unique(analysis)), by = group]
            
            byCols <- c("group", "formula")
            if (haveMSMS)
                byCols <- c(byCols, "frag_formula")
            
            # Determine coverage of formulas within analyses
            formTable[, anaCoverage := .N / anaCount, by = byCols]
            if (formThreshold > 0)
                formTable <- formTable[anaCoverage >= formThreshold] # Apply coverage filter
            
            # Remove duplicate entries (do this after coverage!)
            formTable <- unique(formTable, by = byCols)
            
            # order from best to worst (important for max unique formula filter)
            colorder <- c("group", "byMSMS", intersect(names(formTable), formulaScoringColumns()))
            setorderv(formTable, colorder, c(1, 1, rep(-1, length(colorder)-2)))
            
            # filter max unique formulas
            formTable[, form_unique := match(formula, unique(.SD$formula)), by = group]
            formTable <- formTable[form_unique <= maxFormulas]
            if (haveMSMS)
            {
                formTable[, form_unique_frag := match(frag_formula, unique(frag_formula)), by=.(group, byMSMS, formula)]
                formTable <- formTable[form_unique_frag <= maxFragFormulas]
            }
            
            # Remove some uninteresting columns
            formTable[, c("analysis", "intensity", "anaCount", "form_unique") := NULL]
            if (haveMSMS)
                formTable[, form_unique_frag := NULL]
        }
    }
        
    cat("Done!\n")

    return(formulaConsensus(formulas = formTable, algorithm = algorithm(formList)))
}

formConsensusColOrder <- function(fConsTable)
{
    currentCols <- colnames(fConsTable)

    # all possible columns, depending on algorithm(s) used and their settings
    allCols <- c("group", "ret", "mz", "neutral_formula", "formula", "adduct", "formula_mz", "error",
                 "mSigma", "dbe", "rank", "score", "MS_match", "treeScore", "isoScore", "anaCoverage",
                 "listCoverage", "byMSMS", "frag_formula", "frag_intensity", "frag_mz", "frag_formula_mz",
                 "frag_error", "frag_mSigma", "neutral_loss", "frag_dbe", "frag_score", "MSMS_match",
                 "comb_match", "explainedPeaks", "explainedIntensity", "min_intensity", "max_intensity",
                 "ana_min_intensity", "ana_max_intensity")

    # add algorithm specific scoring/anaCoverage columns that may have been created during merging
    scorePos <- which(allCols == "score")
    curScoreCols <- currentCols[grepl("score-", currentCols)]
    if (length(curScoreCols) > 0)
        allCols <- append(allCols, curScoreCols, scorePos)

    covPos <- which(allCols == "anaCoverage")
    curCovCols <- currentCols[grepl("anaCoverage-", currentCols)]
    if (length(curCovCols) > 0)
        allCols <- append(allCols, curCovCols, covPos)

    ret <- currentCols[match(allCols, currentCols, nomatch = 0)] # re-order

    return(ret)
}

getFragmentInfoFromForms <- function(spec, fragFormTable)
{
    if (nrow(fragFormTable) == 0)
        return(data.table(mz = numeric(0), formula = character(0), intensity = numeric(0), PLIndex = numeric(0)))

    fi <- data.table(mz = fragFormTable$frag_mz, formula = fragFormTable$frag_formula)
    fi[, PLIndex := sapply(mz, function(omz) which.min(abs(omz - spec$mz)))] # UNDONE: is this always correct?
    fi[, intensity := spec$intensity[PLIndex]]
}

getFormInfoList <- function(formConsensus, precursor, groupName)
{
    formTable <- formulaTable(formConsensus)[group == groupName & byMSMS == TRUE & formula == precursor]

    if (nrow(formTable) == 0)
        return(NULL)

    precInfo <- formTable[1] # precursor info is duplicated over all fragment rows

    valText <- function(fmt, value) if (!is.null(value) && !is.na(value)) sprintf(fmt, value) else NULL
    ret <- character()
    ret <- c(ret, valText("Ion formula: %s", precInfo$formula))
    ret <- c(ret, valText("Neutral formula: %s", precInfo$neutral_formula))
    ret <- c(ret, valText("Error: %.2f ppm", precInfo$error))
    ret <- c(ret, valText("mSigma: %.1f", precInfo$mSigma))
    ret <- c(ret, valText("dbe: %f", precInfo$error))
    ret <- c(ret, valText("Rank: %f", precInfo$rank))
    ret <- c(ret, valText("Score: %.2f", precInfo$score))
    ret <- c(ret, valText("Score SIRIUS: %.2f", precInfo[["score-SIRIUS"]]))
    ret <- c(ret, valText("Score DataAnalysis: %.2f", precInfo[["score-Bruker_DataAnalysis"]]))
    ret <- c(ret, valText("MS match: %.2f", precInfo$MS_match))
    ret <- c(ret, valText("MSMS match: %.2f", precInfo$MSMS_match))
    ret <- c(ret, valText("comb match: %.2f", precInfo$comb_match))
    ret <- c(ret, valText("treeScore: %.2f", precInfo$treeScore))
    ret <- c(ret, valText("isoScore: %.2f", precInfo$isoScore))

    return(ret)
}
