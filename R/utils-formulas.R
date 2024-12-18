# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include utils.R
NULL

splitFormulaToList <- memoise(function(formula)
{
    if (!nzchar(formula))
        return(numeric())
        
    # NOTE: dash ('-') added to [:digit:] to allow minus sign (which MF manages to report once in a while...)
    
    # split string in pairs of elements+element counts (and optionally isotopic info), e.g.: { "C30", "[13]C2" }
    spltform <- unlist(regmatches(formula, gregexpr("(\\[[[:digit:]]+\\])?[[:upper:]]{1}[[:lower:]]?[[:digit:]-]*", formula)))
    
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
    return(rbindlist(lapply(formula, function(f)
    {
        if (is.na(f))
        {
            ret <- data.table()
            ret[, (elements) := NA_integer_]
            return(ret)
        }
            
        fl <- as.list(splitFormulaToList(f))
        fl[setdiff(elements, names(fl))] <- 0
        return(fl[elements])
    })))
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
    doSuperSub <- function(f)
    {
        # superscript isotopes
        f <- sub("\\*$", "", gsub("\\[([[:digit:]-]+)\\]", "^\\1*", f))
        
        # subscript element counts
        f <- sub("\\*$", "", gsub("([[:alpha:]]+)([[:digit:]-]+)", "\\1[\\2]*", f))
        
        # HACK: for isotopes specified directly after separated chunk (*) we need to add empty character ('') as it
        # seems we cannot start with superscripted text
        f <- gsub("\\*\\^", "*''^", f)
        return(f)
    }
    
    exprs <- doSuperSub(formulas)
    
    if (!is.null(formulas2))
    {
        exprs2 <- doSuperSub(formulas2)
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
subscriptFormulaHTML <- function(formulas, charges = TRUE)
{
    # isotopes
    formulas <- gsub("\\[([[:digit:]-]+)\\]", "<sup>\\1</sup>", formulas)
    
    if (charges)
        formulas <- gsub("(\\+|\\-)+", "<sup>\\1</sup>", formulas)
    
    # element counts
    formulas <- gsub("([[:alpha:]]+)([[:digit:]-]+)", "\\1<sub>\\2</sub>", formulas)
    
    return(formulas)
}

getFormulaDiffText <- function(form1, form2)
{
    sfl <- splitFormulaToList(subtractFormula(form1, form2))
    ret <- ""
    subfl <- sfl[sfl < 0]
    if (length(subfl) > 0)
        ret <- paste0("-", formulaListToString(abs(subfl)))
    addfl <- sfl[sfl > 0]
    if (length(addfl) > 0)
        ret <- if (nzchar(ret)) paste0(ret, " +", formulaListToString(addfl)) else paste0("+", formulaListToString(addfl))
    return(ret)
}

verifyFormulas <- function(formulas)
{
    data("isotopes", package = "enviPat", envir = environment())
    formulas[is.na(formulas)] <- "" # check_chemform() doesn't handle NAs
    
    # NOTE: check_chemform() errors when there are whitespaces in formulas. So far verifyFormulas() is only called on
    # data from MS libraries, which clearout spaces to avoid this. If verifyFormulas() is called from other code, this
    # should be dealt with...
    
    return(!enviPat::check_chemform(isotopes, formulas)$warning)
}

averageFormulas <- function(formulas)
{
    fltab <- rbindlist(lapply(formulas, function(f) as.list(splitFormulaToList(f))), fill = TRUE)
    for (j in seq_along(fltab))
        set(fltab, which(is.na(fltab[[j]])), j, 0)
    fl <- round(colMeans(fltab))
    return(formulaListToString(fl))
}

# based on RCDK PerformanceNote vignette (https://github.com/CDK-R/cdkr/blob/95f7e23591381c9f0aa75ecc3c0aa6da82f71ebf/rcdk/vignettes/PerformanceNotes.Rmd#L127)
formulaMW <- function(formula)
{
    mfm <- rJava::J("org/openscience/cdk/tools/manipulator/MolecularFormulaManipulator")
    sco <- rJava::J("org.openscience.cdk.silent.SilentChemObjectBuilder")
    
    bi <- rJava::.jcall(sco, "Lorg/openscience/cdk/interfaces/IChemObjectBuilder;", "getInstance")
    mform  <- rJava::.jcall(mfm, "Lorg/openscience/cdk/interfaces/IMolecularFormula;", "getMolecularFormula", formula,
                            bi)
    return(rJava::.jcall(mfm, "D", "getMass", mform))
}


countUniqueFormulas <- function(fList) sum(unlist(lapply(fList, function(ft) length(unique(ft$neutral_formula)))))

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
        if (length(missingElements) > 0 && any(sapply(missingElements, function(mel) minElFL[mel] > 0)))
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

#' Scorings terms for formula candidates
#'
#' Returns a \code{data.frame} with information on which scoring terms are used and what their algorithm specific name
#' is.
#' 
#' @seealso generateFormulas
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

formScoreNames <- function(onlyNums) formulaScorings()$name

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
    scoreCols <- getAllMergedConsCols(formScoreNames(FALSE), names(formTable), mConsNames)
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
        avgCols <- intersect(c(formScoreNames(FALSE), "error", "error_median"), names(formTable))
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
    {
        resultRow[, neutral_formula := subscriptFormulaHTML(neutral_formula)]
        ionFormCols <- getAllMergedConsCols("ion_formula", names(resultRow), mConsNames)
        resultRow[, (ionFormCols) := lapply(mget(ionFormCols), subscriptFormulaHTML)]
    }
    
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
    if (is.null(formTable[["ion_formula_mz"]]))
        formTable[, ion_formula_mz := sapply(neutralMass, calculateMasses, as.adduct(..adduct), type = "mz")]
    
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

setFormulaPLID <- function(formList, MSPeakLists, absAlignMzDev)
{
    # sync fragInfos with group averaged peak lists: the fragments may either have a slightly different m/z than
    # what was used to get the feature formula, or it may simply be removed.
    
    warnAlign <- warnDiffAnn <- FALSE
    
    formList <- Map(formList, lapply(averagedPeakLists(MSPeakLists)[names(formList)], "[[", "MSMS"), f = function(fc, spec)
    {
        fc <- copy(fc)
        fc[explainedPeaks > 0, fragInfo := lapply(fragInfo, function(fi)
        {
            # verify presence
            fi <- fi[sapply(mz, function(fimz)
            {
                wh <- which(numLTE(abs(fimz - spec$mz), absAlignMzDev))
                if (length(wh) > 1)
                    warnAlign <<- TRUE
                return(length(wh) > 0)
            })]
            
            if (nrow(fi) > 0)
            {
                # align remaining mzs
                fi[, PLID := sapply(mz, function(x) spec[which.min(abs(x - mz))]$ID)]
                fi[, mz := spec[match(PLID, ID)]$mz]
                if (anyDuplicated(fi$PLID))
                {
                    warnDiffAnn <<- TRUE
                    fi[, mzD := abs(mz - ion_formula_mz)]
                    setorderv(fi, c("PLID", "mzD"))
                    fi <- unique(fi, by = "PLID")[, -"mzD"]
                }
            }
            
            return(fi)
        })]
        fc[, explainedPeaks := sapply(fragInfo, nrow)]
        fc[explainedPeaks == 0, fragInfo := lapply(fragInfo, function(fi) copy(fi)[, PLID := integer()])]
        
        return(fc)
    })
    formList <- pruneList(formList, checkZeroRows = TRUE)
    
    if (warnAlign)
        warning("Found multiple MS/MS peak list m/z values that may correspond to formula fragment m/z. ",
                "Consider lowering absAlignMzDev.", call. = FALSE)
    if (warnDiffAnn)
        warning("Matched different formula fragment annotations to single MS/MS peak. ",
                "Taking annotation with lowest m/z deviation. ",
                "Consider lowering absAlignMzDev and/or mass tolerances.", call. = FALSE)
    
    return(formList)
}
