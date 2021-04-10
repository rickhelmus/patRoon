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

checkFormula <- function(formula, elementsVec)
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

        if (OK)
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
    data.frame(name = c("combMatch", "isoScore", "mSigma", "MSMSScore", "score"),
               genform = c("comb_match", "MS_match", "-", "MSMS_match", "-"),
               sirius = c("-", "isoScore", "-", "treeScore", "score"),
               bruker = c("-", "-", "mSigma", "-", "Score"),
               description = c("MS and MS/MS combined match value", "How well the isotopic pattern matches",
                               "Deviation of the isotopic pattern", "How well MS/MS data matches",
                               "Overall MS formula score"),
               stringsAsFactors = FALSE)
}

rankFormulaTable <- function(formTable)
{
    # order from best to worst

    rankCols <- intersect(c("score", "combMatch", "isoScore", "mSigma", "MSMSScore"), names(formTable))

    colorder <- rep(-1, length(rankCols))

    # low mSigma values are best
    colorder[grepl("^mSigma", names(rankCols))] <- 1

    setorderv(formTable, rankCols, colorder)
    return(formTable)
}

calculateFormScoreRanges <- function(formTable, mConsNames)
{
    scoreCols <- getAllMergedConsCols(formulaScorings()$name, names(formTable), mConsNames)
    return(setNames(lapply(scoreCols,
                           function(sc) if (any(!is.na(formTable[[sc]]))) range(formTable[[sc]], na.rm = TRUE) else c(NA, NA)),
                    scoreCols))
}

generateFormConsensusForGroup <- function(formList, mergeCount, formThreshold, formThresholdAnn)
{
    # merge all together
    formTable <- rbindlist(formList, fill = TRUE, idcol = "analysis")
    
    if (nrow(formTable) > 0)
    {
        # Determine coverage of precursor formulas.
        formTable[, featCoverage := uniqueN(analysis) / mergeCount, by = "neutral_formula"]
        formTable[, featCoverageAnn := uniqueN(analysis) / length(formList), by = "neutral_formula"]
        
        # add column that specifies of which datasets its merged from
        formTable[, analyses := paste0(unique(analysis), collapse = ","), by = "neutral_formula"]
        
        # Apply coverage filters
        if (formThreshold > 0 || formThresholdAnn > 0)
            formTable <- formTable[featCoverage >= formThreshold & featCoverageAnn >= formThresholdAnn]
        
        # average scorings
        avgCols <- intersect(c(formulaScorings()$name, "error", "error_median"), names(formTable))
        formTable[, (avgCols) := lapply(unique(.SD, by = "analysis")[, avgCols, with = FALSE], mean, na.rm = TRUE),
                  by = "neutral_formula", .SDcols = c(avgCols, "analysis")]
        
        # merge fragInfos: average scorings and keep one row per fragment
        formTable[, fragInfo := {
            mergedFI <- rbindlist(fragInfo)
            if (nrow(mergedFI) > 0)
            {
                for (col in c("error", "mSigma", "score"))
                {
                    if (!is.null(mergedFI[[col]]))
                        mergedFI[, (col) := mean(get(col)), by = "ion_formula"]
                }
                mergedFI <- unique(mergedFI, by = "ion_formula")
                setorderv(mergedFI, "mz")
            }
            list(list(mergedFI)) # double wrap seems necessary?
        }, by = "neutral_formula"]
        formTable[, explainedPeaks := sapply(fragInfo, nrow)]
        
        # Remove duplicate entries (do this after coverage!)
        formTable <- unique(formTable, by = "neutral_formula")
        
        formTable <- rankFormulaTable(formTable)
        
        formTable[, analysis := NULL]
    }
    
    return(formTable)
}

generateGroupFormulasByConsensus <- function(formList, mergeCounts, formThreshold, formThresholdAnn, origGNames)
{
    cat("Generating formula consensus...\n")
    
    hash <- makeHash(formList, mergeCounts, formThreshold, formThresholdAnn, origGNames)
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
            ret <- generateFormConsensusForGroup(fl, mergeCounts[[resNames[[i]]]], formThreshold, formThresholdAnn)
            setTxtProgressBar(prog, i)
            return(ret)
        })
        names(formCons) <- resNames
        formCons <- pruneList(formCons, checkZeroRows = TRUE)
        
        # fix order
        formCons <- formCons[intersect(origGNames, names(formCons))]
        
        setTxtProgressBar(prog, resCount)
        close(prog)
        
        saveCacheData("formulasFGroupConsensus", formCons, hash)
    }
    else
        cat("Done!\n")
    
    return(formCons)
}

getFormInfoList <- function(formTable, index, mConsNames, useHTML)
{
    resultRow <- formTable[index]

    if (useHTML)
        resultRow[, neutral_formula := subscriptFormulaHTML(neutral_formula)]

    addValText <- function(curText, fmt, col)
    {
        cols <- getAllMergedConsCols(col, names(resultRow), mConsNames)

        ret <- character()
        for (cl in cols)
        {
            if (!is.null(resultRow[[cl]]) && !is.na(resultRow[[cl]]) &&
                (!is.character(resultRow[[cl]]) || nzchar(resultRow[[cl]])))
            {
                fm <- sprintf("%s: %s", cl, fmt)
                ret <- c(ret, sprintf(fm, resultRow[[cl]]))
            }
        }

        return(c(curText, ret))
    }

    ret <- character()

    ret <- addValText(ret, "%s", "neutral_formula")
    ret <- addValText(ret, "%s", "ion_formula")
    ret <- addValText(ret, "%.4f", "neutralMass")
    ret <- addValText(ret, "%.4f", "ion_formula_mz")
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

addMiscFormulaInfo <- function(formTable, adduct)
{
    formTable <- copy(formTable)
    
    if (is.null(formTable[["ion_formula"]]))
        formTable[, ion_formula := calculateIonFormula(neutral_formula, ..adduct)]
    if (is.null(formTable[["neutralMass"]]))
        formTable[, neutralMass := sapply(neutral_formula, getFormulaMass)]
    
    formTable[, fragInfo := Map(ion_formula, fragInfo, f = function(form, fi)
    {
        fi <- copy(fi)
        if (nrow(fi) == 0)
            fi[, neutral_loss := character()]
        else
            fi[, neutral_loss := sapply(ion_formula, subtractFormula, formula1 = form)]
        return(fi)
    })]
    
    if (is.null(formTable[["explainedPeaks"]]))
        formTable[, explainedPeaks := sapply(fragInfo, nrow)]
    
    return(formTable)
}

setFormulaPLIndex <- function(formList, MSPeakLists, absAlignMzDev)
{
    # sync fragInfos with group averaged peak lists: the fragments may either have a slightly different m/z than
    # what was used to get the feature formula, or it may simply be removed.
    
    formList <- Map(formList, lapply(averagedPeakLists(MSPeakLists)[names(formList)], "[[", "MSMS"), f = function(fc, spec)
    {
        fc <- copy(fc)
        fc[explainedPeaks > 0, fragInfo := lapply(fragInfo, function(fi)
        {
            # verify presence
            fi <- fi[sapply(mz, function(mz)
            {
                wh <- which(numLTE(abs(mz - spec$mz), absAlignMzDev))
                if (length(wh) > 1)
                    warning("Found multiple MS/MS peak list m/z values that may correspond to formula fragment m/z. ",
                            "Consider lowering absAlignMzDev.", call. = FALSE)
                return(length(wh) > 0)
            })]
            
            if (nrow(fi) > 0)
            {
                # align remaining mzs
                fi[, PLIndex := sapply(mz, function(x) which.min(abs(x - spec$mz)))]
                fi[, mz := spec$mz[PLIndex]]
            }
            
            return(fi)
        })]
        fc[, explainedPeaks := sapply(fragInfo, nrow)]
        fc[explainedPeaks == 0, fragInfo := lapply(fragInfo, function(fi) copy(fi)[, PLIndex := integer()])]
        
        return(fc)
    })
    formList <- pruneList(formList, checkZeroRows = TRUE)
    
    return(formList)
}
