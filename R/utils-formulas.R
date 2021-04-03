#' @include utils.R
NULL

splitFormulaToList <- memoise(function(formula)
{
    if (!nzchar(formula))
        return(numeric())
        
    # NOTE: dash ('-') added to [:digit:] to allow minus sign (which MF manages to report once in a while...)
    
    # split string in pairs of elements+element counts (and optionally isotopic info), e.g.: { "C30", "^13C2" }
    spltform <- unlist(regmatches(formula, gregexpr("(\\^[[:digit:]-]+)?[[:upper:]]{1}[[:lower:]]?[[:digit:]-]*", formula)))
    
    # add '1' to pairs without a count, e.g. "Br" --> "Br1"
    reglist <- regexpr("[[:digit:]-]+$", spltform)
    spltform[reglist == -1] <- paste0(spltform[reglist == -1], "1")
    
    # extract all element counts
    ret <- as.numeric(unlist(regmatches(spltform, gregexpr("[[:digit:]-]+$", spltform))))
    
    # extract all elements (ie remove all counts)
    names(ret) <- gsub("[[:digit:]-]+$", "", spltform)
    
    if (anyDuplicated(names(ret)) != 0)
    {
        unel <- unique(names(ret))
        sortedr <- numeric()
        sortedr[unel] <- sapply(unel, function(e) sum(ret[names(ret) == e]))
        ret <- sortedr
    }
    
    return(ret)
})

getElements <- function(formula, elements)
{
    ret <- sapply(formula, function(f)
    {
        fl <- splitFormulaToList(f)
        df <- as.data.frame(t(as.data.frame(fl)))
        df[setdiff(elements, names(fl))] <- 0
        return(df[, elements, drop = FALSE])
    }, simplify = F, USE.NAMES = F)

    ret <- as.data.frame(do.call(rbind, ret))
    rownames(ret) <- NULL
    return(ret)
}

formulaListToString <- function(formlist)
{
    el <- names(formlist)

    # make element - element count pairs (don't put element counts between 0-1)
    return(paste0(fifelse(formlist < 0 | formlist > 1, paste0(el, formlist), el), collapse = ""))
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

sortFormula <- function(formula)
{
    if (!nzchar(formula))
        return("")
    
    fl <- splitFormulaToList(formula)
    el <- names(fl)
    hasC <- "C" %in% el; hasH <- "H" %in% el
    # Hill order: C and H first then alphabetical. If no C: all alphabetical
    if (hasC)
    {
        el <- el[!el %in% c("C", "H")]
        el <- el[order(el)]
        if (hasH)
            el <- c("H", el)
        el <- c("C", el)
    }
    else
        el <- el[order(el)]
    return(formulaListToString(fl[el]))
}

simplifyFormula <- function(formula) sortFormula(formula)

# Based on chemistry2expression() from ReSOLUTION package
# (authors: Emma Schymanski / Steffen Neumann). See https://github.com/schymane/ReSOLUTION/
subscriptFormula <- function(formulas, over = NULL, formulas2 = NULL, parse = TRUE)
{
    exprs <- sub("\\*$", "", gsub("([[:digit:]-]+)", "[\\1]*", formulas))
    if (!is.null(formulas2))
    {
        exprs2 <- sub("\\*$", "", gsub("([[:digit:]-]+)", "[\\1]*", formulas2))
        exprs <- paste0(exprs, "*'/'*", exprs2)
    }
    
    if (parse)
    {
        if (!is.null(over))
            return(parse(text = sprintf('atop("%s", %s)', over, exprs)))
        return(parse(text = exprs))
    }
    return(exprs)
}

# as above, but for HTML
subscriptFormulaHTML <- function(formulas) gsub("([0-9]+)", "<sub>\\1</sub>", formulas)

averageFormulas <- function(formulas)
{
    fltab <- rbindlist(lapply(formulas, function(f) as.list(splitFormulaToList(f))), fill = TRUE)
    for (j in seq_along(fltab))
        set(fltab, which(is.na(fltab[[j]])), j, 0)
    fl <- round(colMeans(fltab))
    return(formulaListToString(fl))
}

countUniqueFormulas <- function(fList) sum(unlist(lapply(fList, function(ft) length(unique(ft$neutral_formula)))))

addElementInfoToFormTable <- function(formTable, elements, fragElements, OM)
{
    # ensure CHNOPS counts are present
    if (OM)
        elements <- unique(c(if (is.null(elements)) c() else elements, c("C", "H", "N", "O", "P", "S")))

    if (!is.null(elements) && length(elements) > 0)
    {
        # Retrieve element lists from formulas
        el <- getElements(formTable$neutral_formula, elements)
        formTable[, names(el) := el]
    }
    if (!is.null(fragElements) && !is.null(formTable[["frag_formula"]]) &&
        length(fragElements) > 0)
    {
        el <- getElements(formTable$frag_formula, fragElements)
        formTable[, (paste0("frag_", names(el))) := el]
    }

    if (OM)
    {
        # add element ratios commonly used for plotting
        elrat <- function(el1, el2) ifelse(el2 == 0, 0, el1 / el2)
        formTable[, c("OC", "HC", "NC", "PC", "SC") :=
                      .(elrat(O, C), elrat(H, C), elrat(N, C), elrat(P, C), elrat(S, C))]

        # aromaticity index and related DBE (see Koch 2016, 10.1002/rcm.7433)
        formTable[, DBE_AI := 1 + C - O - S - 0.5 * (N + P + H)]
        getAI <- function(dbe, cai) ifelse(cai == 0, 0, dbe / cai)
        formTable[, AI := getAI(DBE_AI, (C - O - N - S - P))]

        formTable[, classification := Vectorize(classifyFormula)(OC, HC, NC, AI)]
    }

    return(formTable)
}

# classification according to Abdulla 2013 (10.1021/ac303221j)
classifyFormula <- function(OC, HC, NC, AI)
{
    if (OC <= 0.2 && HC >= 1.7 && HC <= 2.2)
        return("lipid")
    if (OC > 0.2 && OC <= 0.6 && HC >= 1.7 && HC <= 2.2 && NC > 0.05)
        return("protein")
    if (OC > 0.6 && OC <= 1.2 && HC >= 1.5 && HC <= 2.2)
        return("carbohydrate")
    if (OC > 0.1 && OC <= 0.6 && HC >= 0.6 && HC <= 1.7 && AI < 0.67)
        return("lignin_CRAM")
    if (OC > 0.6 && OC <= 1.2 && HC >= 0.5 && HC <= 1.5 && AI < 0.67)
        return("tannin")
    if (OC <= 0.1 && HC >= 0.7 && HC <= 1.5)
        return("unsat_hydrocarbon")
    if (OC <= 0.1 && HC >= 0.3 && HC <= 0.7 && AI >= 0.67)
        return("condensed_aromatic")

    return("other")
}

checkFormula <- function(formula, elementsVec, negate)
{
    for (elements in elementsVec)
    {
        # any ranges specified?
        if (grepl("[0-9]+\\-[0-9]+", elements))
        {
            minElements <- gsub("([0-9]+)\\-[0-9]+", "\\1", elements)
            maxElements <- gsub("[0-9]+\\-([0-9]+)", "\\1", elements)
            minElFL <- splitFormulaToList(minElements)
            maxElFL <- splitFormulaToList(maxElements)
        }
        else
            minElFL <- maxElFL <- splitFormulaToList(elements)

        formlist <- splitFormulaToList(formula)

        OK <- TRUE

        missingElements <- setdiff(names(minElFL), names(formlist))
        if (length(missingElements) > 0 &&
            any(sapply(missingElements, function(mel) minElFL[mel] > 0)))
            OK <- FALSE
        else
        {
            for (el in names(formlist))
            {
                elc <- formlist[el]

                if (el %in% names(minElFL))
                {
                    if (elc < minElFL[el] || elc > maxElFL[el])
                    {
                        OK <- FALSE
                        break
                    }
                }
            }
        }

        if (OK != negate)
            return(TRUE)
    }

    return(FALSE)
}

#' @details \code{formulaScorings} returns a \code{data.frame} with information
#'   on which scoring terms are used and what their algorithm specific name is.
#' @rdname formula-generation
#' @export
formulaScorings <- function()
{
    data.frame(name = c("combMatch", "frag_mSigma", "frag_score", "isoScore", "mSigma", "MSMSScore", "score"),
               genform = c("comb_match", "-", "-",  "MS_match", "-", "MSMS_match", "-"),
               sirius = c("-", "-", "-", "isoScore", "-", "treeScore", "score"),
               bruker = c("-", "mSigma (SmartFormula3D)", "Score (SmartFormula3D)", "-", "mSigma", "-", "Score"),
               description = c("MS and MS/MS combined match value", "Deviation of isotopic pattern of fragment",
                               "MS/MS fragment score", "How well the isotopic pattern matches", "Deviation of the isotopic pattern",
                               "How well MS/MS data matches", "Overall MS formula score"),
               stringsAsFactors = FALSE)
}

formulaRankingColumns <- function() c("byMSMS", "score", "combMatch", "isoScore", "mSigma", "MSMSScore")

rankFormulaTable <- function(formTable, mConsNames, rankCols = NULL)
{
    # order from best to worst

    if (is.null(rankCols))
        rankCols <- formulaRankingColumns()
    
    rankCols <- getAllMergedConsCols(rankCols, names(formTable), mConsNames)

    colorder <- rep(-1, length(rankCols))

    # low mSigma values are best
    colorder[grepl("^mSigma", names(rankCols))] <- 1

    setorderv(formTable, rankCols, colorder)
    return(formTable)
}

getPrecursorFormData <- function(formTable, cols)
{
    # formula tables may contain multiple lines for a precursor candidate: each
    # for every annotated fragment. The scorings, however, should be the same
    # for each candidate. Nevertheless, when dealing with consensus results,
    # scorings are not homogeneous if a fragment is not annotated by all
    # algorithms. For this reason try to find the first non-NA value from all
    # rows.
    return(formTable[, lapply(.SD, function(x) { xn <- x[!is.na(x)]; if (length(xn) > 0) xn[1] else x[1] }),
                     by = "neutral_formula", .SDcols = cols])
}

normalizeFormScores <- function(formResults, scoreRanges, mConsNames, minMaxNormalization, exclude = NULL)
{
    columns <- names(formResults)
    scoreCols <- getAllMergedConsCols(formulaScorings()$name, columns, mConsNames)
    
    if (!is.null(exclude))
        scoreCols <- scoreCols[!scoreCols %in% getAllMergedConsCols(exclude, columns, mConsNames)]
    
    if (length(scoreCols) > 0)
    {
        scoreRanges <- scoreRanges[scoreCols]
        sc <- getPrecursorFormData(formResults, scoreCols)
        sc[, (scoreCols) := mapply(.SD, scoreRanges, SIMPLIFY = FALSE,
                                   FUN = function(sc, scr) normalize(sc, minMaxNormalization, scr)),
           .SDcols = scoreCols]
        # add/merge normalized columns and restore original column order
        formResults <- formResults[, setdiff(columns, scoreCols), with = FALSE][sc, on = "neutral_formula"][, columns, with = FALSE]
    }
    
    return(formResults)
}

calculateFormScoreRanges <- function(formTable, mConsNames)
{
    scoreCols <- getAllMergedConsCols(formulaScorings()$name, names(formTable), mConsNames)
    return(setNames(lapply(scoreCols,
                           function(sc) if (any(!is.na(formTable[[sc]]))) range(formTable[[sc]], na.rm = TRUE) else c(NA, NA)),
                    scoreCols))
}

generateFormConsensusForGroup <- function(formList, mergeCount, formThreshold, formThresholdAnn, fromCol, mergedAllCol,
                                          mergeCovCol, mergeCovAnnCol, mConsNames, rankCols)
{
    # merge all together
    formTable <- rbindlist(formList, fill = TRUE, idcol = fromCol)
    haveMSMS <- "frag_formula" %in% colnames(formTable)

    if (nrow(formTable) > 0)
    {
        if (!is.null(rankCols)) # HACK: needed by set consensus
        {
            # get rank scores. NOTE: do this before removing any candidates as otherwise ranking will be messed up
            
            rcols <- intersect(rankCols, names(formTable))
            formTable[, rankScore := {
                allRanks <- unlist(getPrecursorFormData(.SD, rcols)[, rcols, with = FALSE])
                invRanks <- (.NGRP - (allRanks - 1)) / .NGRP
                invRanks[is.na(invRanks)] <- 0
                mean(invRanks)
            }, .SDcols = c(rcols, "neutral_formula"), by = "neutral_formula"]
        }
        
        # Determine coverage of precursor formulas.
        formTable[, (mergeCovCol) := uniqueN(get(fromCol)) / mergeCount, by = "neutral_formula"]
        formTable[, (mergeCovAnnCol) := uniqueN(get(fromCol)) / length(formList), by = "neutral_formula"]

        # add column that specifies of which datasets its merged from
        formTable[, (mergedAllCol) := paste0(unique(get(fromCol)), collapse = ","), by = "neutral_formula"]
        
        # Apply coverage filters
        if (formThreshold > 0 || formThresholdAnn > 0)
            formTable <- formTable[get(mergeCovCol) >= formThreshold & get(mergeCovAnnCol) >= formThresholdAnn]

        # remove MS only formulas if MS/MS candidate is also present (do after coverage filter).
        MSMSForms <- unique(formTable[byMSMS == TRUE, neutral_formula])
        formTable <- formTable[byMSMS == TRUE | !neutral_formula %in% MSMSForms]

        # average scorings
        
        # frag_error is unique per fragment, while other scorings are per precursor candidate
        fragAvgCols <- getAllMergedConsCols("frag_error", names(formTable), mConsNames)
        if (haveMSMS && length(fragAvgCols) > 0)
            formTable[, (fragAvgCols) := lapply(.SD, mean, na.rm = TRUE),
                      by = c("neutral_formula", "frag_formula", fromCol),
                      .SDcols = fragAvgCols]
        
        avgCols <- getAllMergedConsCols(c(formulaScorings()$name, "error", "error_median"),
                                        names(formTable), mConsNames)
        formTable[, (avgCols) := lapply(unique(.SD, by = fromCol)[, avgCols, with = FALSE], mean, na.rm = TRUE),
                  by = "neutral_formula", .SDcols = c(avgCols, fromCol)]
        
        # remove NaNs that may have been introduced due to mean(..., na.rm=TRUE)
        numCols <- names(which(sapply(formTable, is.numeric)))
        formTable[, (numCols) := lapply(.SD, nafill), .SDcols = numCols]
        
        # Remove duplicate entries (do this after coverage!)
        formTable <- unique(formTable, by = intersect(c("neutral_formula", "frag_formula"), names(formTable)))
        
        if (!is.null(rankCols))
        {
            setorderv(formTable, "rankScore", order = -1)
            formTable[, rankScore := NULL]
        }
        else
            formTable <- rankFormulaTable(formTable, mConsNames)
    }

    return(formTable)
}

generateGroupFormulasByConsensus <- function(formList, mergeCounts, formThreshold, formThresholdAnn, origGNames,
                                             fromCol, mergedAllCol, mergeCovCol, mergeCovAnnCol, MSPeakLists,
                                             absAlignMzDev, mConsNames, rankCols = NULL)
{
    cat("Generating formula consensus...\n")

    hash <- makeHash(formList, formThreshold, formThresholdAnn, origGNames, fromCol, mergedAllCol, mergeCovCol,
                     mergeCovAnnCol, MSPeakLists, absAlignMzDev, mConsNames, rankCols)
    formCons <- loadCacheData("formulasFGroupConsensus", hash)

    # figure out names
    resNames <- unique(unlist(lapply(formList, names)))
    resCount <- length(resNames)

    if (resCount == 0)
        formCons <- list()
    else if (is.null(formCons))
    {
        prog <- openProgBar(0, resCount)

        formCons <- lapply(seq_len(resCount), function(i)
        {
            fl <- pruneList(lapply(formList, "[[", resNames[[i]]))
            ret <- generateFormConsensusForGroup(fl, mergeCounts[[resNames[[i]]]], formThreshold, formThresholdAnn,
                                                 fromCol, mergedAllCol, mergeCovCol, mergeCovAnnCol, mConsNames,
                                                 rankCols)
            setTxtProgressBar(prog, i)
            return(ret)
        })
        names(formCons) <- resNames
        formCons <- pruneList(formCons, checkZeroRows = TRUE)
        
        # fix order
        formCons <- formCons[intersect(origGNames, names(formCons))]

        # sync with group averaged peak lists: the fragments may either have a slightly different m/z than what was used
        # to get the feature formula, or it may simply be removed.
        if (!is.null(MSPeakLists)) # NULL for sets, not needed there
        {
            formCons <- Map(formCons, lapply(averagedPeakLists(MSPeakLists)[names(formCons)], "[[", "MSMS"), f = function(fc, spec)
            {
                # NOTE: ignore non MS/MS formulae
                
                if (is.null(fc[["frag_mz"]]))
                    return(fc)
                
                fc <- fc[byMSMS == FALSE | sapply(frag_mz, function(mz)
                {
                    wh <- which(numLTE(abs(mz - spec$mz), absAlignMzDev))
                    if (length(wh) > 1)
                        warning("Found multiple MS/MS peak list m/z values that may correspond to formula fragment m/z. ",
                                "Consider lowering absAlignMzDev.", call. = FALSE)
                    return(length(wh) > 0)
                })]
                
                # align remaining mzs
                fc[byMSMS == TRUE, frag_mz := spec$mz[sapply(frag_mz, function(mz) which.min(abs(mz - spec$mz)))]]
                
                return(fc)
            })
            formCons <- pruneList(formCons, checkZeroRows = TRUE)
        }
        
        setTxtProgressBar(prog, resCount)
        close(prog)

        saveCacheData("formulasFGroupConsensus", formCons, hash)
    }
    else
        cat("Done!\n")

    return(formCons)
}

getFragmentInfoFromForms <- function(spec, fragFormTable)
{
    if (nrow(fragFormTable) == 0)
        return(data.table(mz = numeric(0), formula = character(0), neutral_loss = character(0),
                          intensity = numeric(0), PLIndex = numeric(0)))

    fi <- data.table(mz = fragFormTable$frag_mz, formula = fragFormTable$frag_formula, neutral_loss = fragFormTable$neutral_loss)
    fi[, PLIndex := sapply(mz, function(omz) which.min(abs(omz - spec$mz)))] # UNDONE: is this always correct?
    fi[, intensity := spec$intensity[PLIndex]]

    if (!is.null(fragFormTable[["mergedBy"]]))
        fi[, mergedBy := fragFormTable$mergedBy]

    return(fi)
}

getFormInfoList <- function(formTable, precursor, mConsNames, useHTML = FALSE)
{
    formTable <- formTable[neutral_formula == precursor]

    if (nrow(formTable) == 0)
        return(NULL)

    precInfo <- getPrecursorFormData(formTable, names(formTable))

    if (useHTML)
        precInfo[, neutral_formula := subscriptFormulaHTML(neutral_formula)]

    addValText <- function(curText, fmt, col)
    {
        cols <- getAllMergedConsCols(col, names(precInfo), mConsNames)

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

    ret <- addValText(ret, "%s", "formula_mz")
    ret <- addValText(ret, "%s", "neutral_formula")
    ret <- addValText(ret, "%.2f ppm", "error")
    ret <- addValText(ret, "%.2f ppm", "error_median")
    ret <- addValText(ret, "%.1f", "mSigma")
    ret <- addValText(ret, "%f", "rank")
    ret <- addValText(ret, "%.2f", "score")
    ret <- addValText(ret, "%.2f", "MSMSScore")
    ret <- addValText(ret, "%.2f", "combMatch")
    ret <- addValText(ret, "%.2f", "isoScore")

    return(ret)
}

getFormulaMass <- memoise(function(f, c = 0) rcdk::get.formula(f, c)@mass)
