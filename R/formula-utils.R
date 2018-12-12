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

rankFormulaTable <- function(formTable)
{
    # order from best to worst
    # do grap to handle merged columns
    scCols <- grep(paste0("^(", paste0(formulaScoringColumns(), collapse = "|"), ")"),
                   names(formTable), value = TRUE)
    colorder <- c("byMSMS", scCols)
    setorderv(formTable, colorder, c(1, rep(-1, length(colorder)-1)))
    return(formTable)
}

generateFormConsensusForGroup <- function(formAnaList, formThreshold)
{
    # merge all together
    formTable <- rbindlist(formAnaList, fill = TRUE, idcol = "analysis")
    haveMSMS <- "frag_formula" %in% colnames(formTable)
    
    if (nrow(formTable) > 0)
    {
        # number of analyses searched for formulas per group
        anaCount <- length(formAnaList)
        
        byCols <- "formula"
        if (haveMSMS)
            byCols <- c(byCols, "frag_formula")
        
        # Determine coverage of formulas within analyses
        formTable[, anaCoverage := .N / anaCount, by = byCols]
        if (formThreshold > 0)
            formTable <- formTable[anaCoverage >= formThreshold] # Apply coverage filter
        
        # Remove duplicate entries (do this after coverage!)
        formTable <- unique(formTable, by = byCols)
        
        formTable <- rankFormulaTable(formTable)

        formTable[, "analysis" := NULL]
    }
    
    return(formTable)
}

generateGroupFormulasByConsensus <- function(formList, formThreshold)
{
    cat("Generating feature group formula consensus...\n")
    
    hash <- makeHash(formList, formThreshold)
    formCons <- loadCacheData("formCons", hash)
    
    # figure out feature groups 
    gNames <- unique(unlist(sapply(formList, names, simplify = FALSE), use.names = FALSE))
    gCount <- length(gNames)
    
    if (gCount == 0)
        formCons <- list()
    else if (is.null(formCons))
    {
        prog <- txtProgressBar(0, gCount, style = 3)
        
        formCons <- lapply(seq_len(gCount), function(grpi)
        {
            fAnaList <- lapply(formList, "[[", gNames[[grpi]])
            fAnaList <- fAnaList[!sapply(fAnaList, is.null)]
            
            ret <- generateFormConsensusForGroup(fAnaList, formThreshold)
            setTxtProgressBar(prog, grpi)
            return(ret)
        })
        names(formCons) <- gNames
        
        setTxtProgressBar(prog, gCount)
        close(prog)
        
        saveCacheData("formCons", formCons, hash)
    }
    else
        cat("Done!\n")
    
    return(formCons)
}

formConsensusColOrder <- function(fConsTable)
{
    currentCols <- colnames(fConsTable)

    # all possible columns, depending on algorithm(s) used and their settings
    allCols <- c("group", "ret", "mz", "neutral_formula", "formula", "adduct", "formula_mz", "error",
                 "mSigma", "dbe", "rank", "score", "MS_match", "treeScore", "isoScore", "anaCoverage",
                 "listCoverage", "byMSMS", "frag_neutral_formula", "frag_formula", "frag_intensity", "frag_mz", "frag_formula_mz",
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

getFormInfoList <- function(formTable, precursor)
{
    formTable <- formTable[byMSMS == TRUE & formula == precursor]

    if (nrow(formTable) == 0)
        return(NULL)

    precInfo <- formTable[1] # precursor info is duplicated over all fragment rows

    addValText <- function(curText, fmt, col)
    {
        # get all columns matching value of 'col' as prefix: merged names after
        # consensus will have format 'col-X'.
        cols <- grep(paste0("^", col), names(precInfo), value = TRUE)
        
        ret <- character()
        for (cl in cols)
        {
            if (!is.null(precInfo[[cl]]) && !is.na(precInfo[[cl]]) &&
                (!is.character(precInfo[[cl]]) || nzchar(precInfo[[cl]])))
            {
                fm <- sprintf("%s: %s", cl, fmt)
                ret <- c(ret, sprintf(fm, precInfo[[cl]]))
            }
        }
        
        return(c(curText, ret))
    }
    
    ret <- character()
    
    ret <- addValText(ret, "%s", "formula")
    ret <- addValText(ret, "%s", "neutral_formula")
    ret <- addValText(ret, "%.2f ppm", "error")
    ret <- addValText(ret, "%.1f", "mSigma")
    ret <- addValText(ret, "%f", "error")
    ret <- addValText(ret, "%f", "rank")
    ret <- addValText(ret, "%.2f", "score")
    ret <- addValText(ret, "%.2f", "MS_match")
    ret <- addValText(ret, "%.2f", "MSMS_match")
    ret <- addValText(ret, "%.2f", "comb_match")
    ret <- addValText(ret, "%.2f", "treeScore")
    ret <- addValText(ret, "%.2f", "isoScore")

    return(ret)
}
