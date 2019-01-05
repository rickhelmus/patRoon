#' @include main.R
#' @include workflow-step.R
NULL

#' Formula lists class
#'
#' Contains data of generated chemical formulae for given feature groups.
#'
#' \code{formulas} objects are obtained from \link[=formula-generation]{formula
#' generators}.
#'
#' @slot formulas,featureFormulas Lists of all generated formulae. Use the
#'   \code{formulaTable} method for access.
#' @slot algorithm The algorithm that was used for generation of formulae. Use
#'   the \code{algorithm} method for access.
#'
#' @param obj,x,object,formulas The \code{formulas} object.
#' @param \dots \code{consensus}: One or more \code{formulas} objects that should be used to
#'   generate the consensus.
#'
#'   \code{as.data.frame}: Arguments passed to \code{as.data.table}.
#' @param OM For \code{as.data.table}/\code{as.data.frame}: if set to \code{TRUE} several columns
#'   with information relevant for organic matter (OM) characterization will be
#'   added (e.g. elemental ratios, classification). This will also make sure
#'   that \code{countElements} contains at least C, H, N, O, P and S.
#'
#'   For \code{filter}: If \code{TRUE} then several filters are applied to
#'   exclude unlikely formula candidates present in organic matter (OM). See
#'   Source section for details.
#'
#' @templateVar seli analyses
#' @templateVar selOrderi analyses()
#' @templateVar selj feature groups
#' @templateVar selOrderj groupNames()
#' @templateVar optionalji TRUE
#' @templateVar dollarOpName feature group
#' @template sub_op-args
#'
#' @section Source: Calculation of the aromaticity index (AI) and related double
#'   bond equivalents (DBE_AI) is performed as described in Koch 2015. Formula
#'   classification is performed by the rules described in Abdulla 2013.
#'   Filtering of OM related molecules is performed as described in Koch 2006
#'   and Kujawinski 2006. (see references).
#'
#' @references \insertRef{Koch2015}{patRoon} \cr\cr
#'   \insertRef{Abdulla2013}{patRoon} \cr\cr
#'   \insertRef{Koch2006}{patRoon} \cr\cr
#'   \insertRef{Kujawinski2006}{patRoon}
#'
#' @templateVar class formulas
#' @template class-hierarchy
#'
#' @export
formulas <- setClass("formulas", slots = c(formulas = "list", featureFormulas = "list"),
                     contains = "workflowStep")

#' @describeIn formulas Accessor method to obtain generated formulae.
#'
#' @param features If \code{TRUE} returns formula data for features, otherwise
#'   for feature groups.
#'
#' @return \code{formulaTable} returns a \code{list} containing for each feature
#'   group (or feature if \code{features=TRUE}) a \code{\link{data.table}}
#'   with an overview of all generated formulae and other data such as candidate
#'   scoring and MS/MS fragments.
#'
#' @aliases formulaTable
#' @export
setMethod("formulaTable", "formulas", function(obj, features) if (features) obj@featureFormulas else obj@formulas)

#' @describeIn formulas Accessor method for the algorithm (a character
#'   string) used to generate formulae.
#' @export
setMethod("algorithm", "formulas", function(obj) obj@algorithm)

#' @templateVar class formulas
#' @templateVar what analyses
#' @template strmethod
#' @export
setMethod("analyses", "formulas", function(obj) names(obj@featureFormulas))

#' @templateVar class formulas
#' @templateVar what feature groups
#' @template strmethod
#' @export
setMethod("groupNames", "formulas", function(obj) names(obj@formulas))

#' @describeIn formulas Obtain total number of formulae entries.
#' @export
setMethod("length", "formulas", function(x) sum(unlist(sapply(x@formulas, function(ft) length(unique(ft$formula))))))

#' @describeIn formulas Show summary information for this object.
#' @export
setMethod("show", "formulas", function(object)
{
    callNextMethod()

    ft <- formulaTable(object, TRUE)
    hasFeatForms <- length(ft) > 0
    ftcounts <- if (hasFeatForms) recursiveApplyDT(ft, function(x) length(unique(x$formula)), sapply) else 0
    ma <- mean(sapply(ftcounts, sum))
    mft <- mean(sapply(ftcounts, mean))
    printf("Formulas assigned to features:\n")
    printf("  - Total formula count: %d\n", sum(unlist(ftcounts)))
    printf("  - Average formulas per analysis: %.1f\n", ma)
    printf("  - Average formulas per feature: %.1f\n", mft)

    gft <- formulaTable(object)
    mfg <- if (length(gft) > 0) sapply(gft, nrow) else 0
    printf("Formulas assigned to feature groups:\n")
    printf("  - Total formula count: %d\n", length(object))
    printf("  - Average formulas per feature group: %.1f\n", mean(mfg))
})

#' @describeIn formulas Subset on feature groups.
#' @export
setMethod("[", c("formulas", "ANY", "missing", "missing"), function(x, i, j, ...)
{
    if (!missing(i))
        assertSubsetArg(i)

    if (!missing(i))
    {
        if (!is.character(i))
            i <- groupNames(x)[i]
        x@featureFormulas <- sapply(x@featureFormulas, function(a) pruneList(a[i]),
                                    simplify = FALSE)
        x@featureFormulas <- pruneList(x@featureFormulas, TRUE)

        x@formulas <- pruneList(x@formulas[i], TRUE)
    }

    return(x)
})

#' @describeIn formulas Extract a formula table. If both arguments (\code{i} and
#'   \code{j}) are specified, the feature specific formula table belonging to
#'   the analysis (\code{i})/feature group (\code{j}) is returned. Otherwise the
#'   formula table for the feature group specified by \code{j} is returned.
#' @export
setMethod("[[", c("formulas", "ANY", "ANY"), function(x, i, j)
{
    assertExtractArg(i)
    if (!missing(j))
        assertExtractArg(j)

    if (!missing(j))
    {
        # both arguments specified, return feature formula table

        if (length(x@featureFormulas) == 0)
            stop("This object does not contain formulas for features.")

        if (!is.character(i))
            i <- analyses(x)[i]

        if (!is.character(j))
            j <- groupNames(x)[j]

        return(x@featureFormulas[[c(i, j)]])
    }

    # else return regular feature group formulas

    if (!is.character(i))
        i <- groupNames(x)[i]
    return(x@formulas[[i]])
})

#' @describeIn formulas Extract a formula table for a feature group.
#' @export
setMethod("$", "formulas", function(x, name)
{
    eval(substitute(x@formulas[[NAME_ARG]], list(NAME_ARG = name)))
})

#' @describeIn formulas Generates a table with all candidate formulae for each
#'   feature group and other information such as element counts.
#'
#' @param average If set to \code{TRUE} an 'average formula' is generated for
#'   each feature group by combining all elements from all candidates and
#'   averaging their amounts. This obviously leads to non-existing formulae,
#'   however, this data may be useful to deal with multiple candidate formulae
#'   per feature group when performing elemental characterization.
#' @param countElements,countFragElements A \code{character} vector with
#'   elements that should be counted for each MS(/MS) formula candidate. For
#'   instance, \code{c("C", "H")} adds columns for both carbon and hydrogen
#'   amounts of each formula. Note that the neutral formula
#'   (\code{neutral_formula} column) is used to count elements of non-fragmented
#'   formulae, whereas the charged formula of fragments (\code{frag_formula}
#'   column) is used for fragments. Set to \code{NULL} to not count any
#'   elements.
#' @param maxFormulas,maxFragFormulas Maximum amount of unique candidate
#'   formulae (or fragment formulae) per feature group. Set to \code{NULL} to
#'   ignore.
#'
#' @template as_data_table-args
#' 
#' @return \code{as.data.table} returns a \code{\link{data.table}}.
#'
#' @export
setMethod("as.data.table", "formulas", function(x, fGroups = NULL, average = FALSE, countElements = NULL,
                                                countFragElements = NULL, OM = FALSE,
                                                maxFormulas = NULL, maxFragFormulas = NULL)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", null.ok = TRUE, add = ac)
    checkmate::assertFlag(average, add = ac)
    checkmate::assertCharacter(countElements, min.chars = 1, any.missing = FALSE, null.ok = TRUE, add = ac)
    checkmate::assertCharacter(countFragElements, min.chars = 1, any.missing = FALSE, null.ok = TRUE, add = ac)
    checkmate::assertFlag(OM, add = ac)
    checkmate::reportAssertions(ac)

    ret <- rbindlist(formulaTable(x), fill = TRUE, idcol = "group")
    if (length(ret) == 0)
        return(ret)

    if (!is.null(fGroups))
    {
        ret[, c("ret", "mz") := groupInfo(fGroups)[group, ]]
        setcolorder(ret, c("group", "ret", "mz"))
    }

    if (average)
    {
        # collapse byMSMS: will be TRUE if at least an MS/MS formula candidate was there
        ret[, byMSMS := any(byMSMS), by = "group"]

        ret[, formula_avg_count := length(unique(formula)), by = "group"]

        avgCols <- c("formula", "neutral_formula")
        ret[, (avgCols) := lapply(.SD, function(f) averageFormulas(unique(f))), .SDcols = avgCols, by = "group"]

        # remove columns which don't really make sense anymore
        rmCols <- c("neutral_loss", "error", "formula_mz", "dbe", "anaCoverage",
                    "adduct", "rank", "explainedPeaks", "explainedIntensity",
                    # add any fragment columns
                    grep("^frag_", names(ret), value = TRUE),
                    formulaScorings()$name)

        rmCols <- getAllFormulasCols(rmCols, names(ret))
        if (length(rmCols) > 0)
            ret[, (rmCols) := NULL]

        ret <- unique(ret, by = c("group", "formula"))
    }
    else
    {
        if (!is.null(maxFormulas))
        {
            ret[, unFormNr := match(formula, unique(.SD$formula)), by = "group"]
            ret <- ret[unFormNr <= maxFormulas][, unFormNr := NULL]
        }

        if (!is.null(maxFragFormulas) && any(ret$byMSMS))
        {
            ret[, unFormNr := match(frag_formula, unique(.SD$frag_formula)),
                by = c("group", "byMSMS", "formula")]
            ret <- ret[unFormNr <= maxFragFormulas][, unFormNr := NULL]
        }
    }

    ret <- addElementInfoToFormTable(ret, countElements, countFragElements, OM)

    return(ret[])
})

#' @describeIn formulas Performs rule based filtering on formula results.
#'
#' @param minExplainedFragPeaks Minimum number of fragment peaks that are
#'   explained. Setting this to \samp{1} will remove any MS only formula
#'   results. Set to \code{NULL} to ignore.
#' @param topMost Only retain no more than this amount of best ranked candidates
#'   for each feature group.
#' @param scoreLimits Filter results by their scores. Should be a named
#'   \code{list} that contains two-sized numeric vectors with the
#'   minimum/maximum value of a score (use \code{-Inf}/\code{Inf} for no
#'   limits). The names of each element should follow the values returned by
#'   \code{\link{formulaScorings}()$name}. For instance,
#'   \code{scoreLimits=list(isoScore=c(0.5, Inf))} specifies that the isotopic
#'   match score should be at least \samp{0.5}. More details of scorings can be
#'   obtained with \code{\link{formulaScorings}}. Note that a result without a
#'   specified scoring is never removed. Set to \code{NULL} to skip this filter.
#'
#' @templateVar withLoss TRUE
#' @template element-args
#'
#' @note \code{filter} does not modify any formula results for features (if
#'   present).
#'
#' @return \code{filter} returns a filtered \code{\link{formulas}} object.
#'
#' @export
setMethod("filter", "formulas", function(obj, minExplainedFragPeaks = NULL, elements = NULL,
                                         fragElements = NULL, lossElements = NULL,
                                         topMost = NULL, scoreLimits = NULL,
                                         OM = NULL)
{
    scCols <- formulaScorings()$name

    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertCount, . ~ topMost + minExplainedFragPeaks,
           positive = c(TRUE, FALSE), null.ok = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertCharacter, . ~ elements + fragElements + lossElements,
           min.chars = 1, min.len = 1, null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertList(scoreLimits, null.ok = TRUE, types = "numeric", add = ac)
    if (!is.null(scoreLimits))
    {
        checkmate::assertNames(names(scoreLimits), type = "unique", subset.of = scCols, add = ac)
        checkmate::qassertr(scoreLimits, "N2")
    }
    checkmate::assertLogical(OM, null.ok = TRUE)
    checkmate::reportAssertions(ac)

    cat("Filtering formulas... ")

    oldn <- length(obj)
    obj@formulas <- pruneList(sapply(groupNames(obj), function(grp)
    {
        formTable <- obj[[grp]]
        if (!is.null(minExplainedFragPeaks) && minExplainedFragPeaks > 0)
        {
            formTable <- formTable[byMSMS == TRUE]
            if (nrow(formTable) == 0)
                return(formTable)
            fragCounts <- formTable[, ifelse(byMSMS, length(frag_formula), 0L), by = "formula"][[2]]
            formTable <- formTable[fragCounts >= minExplainedFragPeaks]
        }

        if (!is.null(elements))
            formTable <- formTable[sapply(neutral_formula, checkFormula, elements)]
        if ((!is.null(fragElements) || !is.null(lossElements)))
        {
            formTable <- formTable[byMSMS == TRUE]
            if (nrow(formTable) == 0)
                return(formTable)
            if (!is.null(fragElements))
                formTable <- formTable[formTable[, rep(any(sapply(frag_formula, checkFormula, fragElements)), .N),
                                                 by = "formula"][[2]]]
            if (!is.null(lossElements))
                formTable <- formTable[formTable[, rep(any(sapply(neutral_loss, checkFormula, lossElements)), .N),
                                                 by = "formula"][[2]]]
        }

        if (!is.null(scoreLimits))
        {
            for (sc in names(scoreLimits))
            {
                cols <- getAllFormulasCols(sc, names(formTable))
                if (length(cols) == 0)
                    next
                formTable <- formTable[formTable[, do.call(pmin, c(.SD, list(na.rm = TRUE))) >= scoreLimits[[sc]][1] &
                                                   do.call(pmax, c(.SD, list(na.rm = TRUE))) <= scoreLimits[[sc]][2],
                                                 .SDcols = cols]]
            }
        }

        if (!is.null(OM) && OM)
        {
            fElTable <- addElementInfoToFormTable(copy(formTable), NULL, NULL, OM = TRUE)
            keep <- fElTable[,
                         # rules from Kujawinski & Behn, 2006 (10.1021/ac0600306)
                         H >= 1/3 * C &
                         H <= ((2 * C) + N + 2) &
                         (H + N) %% 2 == 0 &
                         N <= C &
                         O <= C &
                         P <= 2 &
                         S <= 2 &

                         # rules from Koch & dittmar 2006 (10.1002/rcm.2386)
                         sapply(DBE_AI, checkmate::checkInt) &
                         HC <= 2.2 &
                         OC <= 1.2 &
                         NC <= 0.5]
            formTable <- formTable[keep]
        }
        
        if (!is.null(topMost))
        {
            unFormNrs <- formTable[, match(formula, unique(.SD$formula))]
            formTable <- formTable[unFormNrs <= topMost]
        }

        return(formTable)
    }, simplify = FALSE), checkZeroRows = TRUE)

    newn <- length(obj)
    printf("Done! Filtered %d (%.2f%%) formulas. Remaining: %d\n", oldn - newn, if (oldn == 0) 0 else (1-(newn/oldn))*100, newn)
    return(obj)
})

#' @describeIn formulas Plots an annotated spectrum for a given candidate
#'   formula of a feature or feature group.
#'
#' @param precursor The formula of the precursor (in ionic form, \emph{i.e.} as
#'   detected by the MS).
#' @param analysis A \code{character} specifying the analysis for which the
#'   annotated spectrum should be plotted. If \code{NULL} then annotation
#'   results for the complete feature group will be plotted.
#'
#' @template plotSpec-args
#'
#' @template useGGplot2
#'
#' @return \code{plotSpec} will return a \code{\link[=ggplot2]{ggplot object}}
#'   if \code{useGGPlot2} is \code{TRUE}.
#'
#' @export
setMethod("plotSpec", "formulas", function(obj, precursor, groupName, analysis = NULL, MSPeakLists, useGGPlot2 = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertString(precursor, min.chars = 1, add = ac)
    checkmate::assertString(groupName, min.chars = 1, add = ac)
    checkmate::assertString(analysis, min.chars = 1, null.ok = TRUE, add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertFlag(useGGPlot2, add = ac)
    checkmate::reportAssertions(ac)

    if (!is.null(analysis))
    {
        formTable <- obj[[analysis, groupName]]
        spec <- MSPeakLists[[analysis, groupName]][["MSMS"]]
    }
    else
    {
        formTable <- obj[[groupName]]
        spec <- MSPeakLists[[groupName]][["MSMS"]]
    }

    formTable <- formTable[byMSMS == TRUE & formula == precursor]

    if (nrow(formTable) == 0 || is.null(spec))
        return(NULL)

    fi <- getFragmentInfoFromForms(spec, formTable)

    if (useGGPlot2)
        return(makeMSPlotGG(spec, fi) + ggtitle(precursor))

    makeMSPlot(spec, fi, main = precursor)
})

#' @describeIn formulas Generates a consensus of results from multiple
#'   \code{formulas} objects.
#'
#' @param formThreshold Fractional minimum amount (0-1) of which a formula
#'   candidate should be present within all objects. For instance, a value of
#'   \samp{0.5} means that a particular formula should be present in at least
#'   \samp{50\%} of all objects.
#'
#' @return \code{consensus} returns a \code{formulas} object that is produced by
#'   merging results from multiple \code{formulas} objects.
#'
#' @export
setMethod("consensus", "formulas", function(obj, ..., formThreshold = 0)
{
    allFormulas <- c(list(obj), list(...))

    ac <- checkmate::makeAssertCollection()
    checkmate::assertList(allFormulas, types = "formulas", min.len = 2, any.missing = FALSE,
                          unique = TRUE, .var.name = "...", add = ac)
    checkmate::assertNumber(formThreshold, lower = 0, finite = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    allFormNames <- sapply(allFormulas, algorithm)
    allFormNames <- make.unique(allFormNames)

    allFormulasLists <- sapply(seq_along(allFormulas), function(fi)
    {
        return(lapply(formulaTable(allFormulas[[fi]]), function(ft)
        {
            ret <- copy(ft)
            ret[, mergedBy := allFormNames[fi]]
            setnames(ret, paste0(names(ret), "-", allFormNames[fi]))
            return(ret)
        }))

    }, simplify = FALSE)

    # UNDONE: remove old style columns?
    uniqueCols <- c("neutral_formula", "formula_mz", "error", "dbe", "frag_mz", "frag_neutral_formula",
                    "frag_formula_mz", "frag_error", "neutral_loss", "frag_dbe", "min_intensity", "max_intensity",
                    "ana_min_intensity", "ana_max_intensity")

    consFormulaList <- allFormulasLists[[1]]
    leftName <- allFormNames[1]
    for (righti in seq(2, length(allFormulasLists)))
    {
        rightName <- allFormNames[righti]

        printf("Merging %s with %s... ", paste0(allFormNames[seq_len(righti-1)], collapse = ","), rightName)

        rightFList <- allFormulasLists[[righti]]

        for (grp in union(names(consFormulaList), names(rightFList)))
        {
            if (is.null(rightFList[[grp]]))
                next # nothing to merge
            else if (is.null(consFormulaList[[grp]])) # not yet present
            {
                mTable <- rightFList[[grp]]

                # rename columns that should be unique from right to left
                unCols <- c(uniqueCols, c("formula", "byMSMS", "frag_formula", "mergedBy"))
                unCols <- unCols[sapply(unCols, function(uc) !is.null(mTable[[paste0(uc, "-", rightName)]]))]
                setnames(mTable, paste0(unCols, "-", rightName), paste0(unCols, "-", leftName))
            }
            else
            {
                haveLeftMSMS <- paste0("frag_formula-", leftName) %in% names(consFormulaList[[grp]])
                haveRightMSMS <- paste0("frag_formula-", rightName) %in% names(rightFList[[grp]])

                mergeCols <- c("formula", "byMSMS") # put byMSMS in there anyway in case only left/right has MSMS
                if (haveLeftMSMS && haveRightMSMS)
                    mergeCols <- c(mergeCols, "frag_formula")
                mTable <- merge(consFormulaList[[grp]], rightFList[[grp]], all = TRUE,
                                by.x = paste0(mergeCols, "-", leftName),
                                by.y = paste0(mergeCols, "-", rightName))

                if (!haveLeftMSMS && haveRightMSMS)
                    setnames(mTable, paste0("frag_formula-", rightName), paste0("frag_formula-", leftName))

                # remove duplicate columns that shouldn't
                for (col in uniqueCols)
                {
                    colLeft <- paste0(col, "-", leftName)
                    colRight <- paste0(col, "-", rightName)
                    if (!is.null(mTable[[colRight]]))
                    {
                        if (is.null(mTable[[colLeft]]))
                            setnames(mTable, colRight, colLeft)
                        else
                        {
                            mTable[, (colLeft) := ifelse(!is.na(get(colLeft)), get(colLeft), get(colRight))]
                            mTable[, (colRight) := NULL]
                        }
                    }
                }

                # collapse mergedBy
                ml <- paste0("mergedBy-", leftName); mr <- paste0("mergedBy-", rightName)
                mTable[!is.na(get(ml)), (ml) := ifelse(!is.na(get(mr)), paste(get(ml), get(mr), sep = ","), get(ml))]
                mTable[is.na(get(ml)), (ml) := get(mr)]
                mTable[, (mr) := NULL]
            }

            consFormulaList[[grp]] <- mTable
        }

        cat("Done!\n")
    }

    printf("Determining coverage... ")

    for (grpi in seq_along(consFormulaList))
    {
        # fix up de-duplicated column names
        deDupCols <- c(uniqueCols, "formula", "byMSMS", "frag_formula", "mergedBy")
        leftCols <- paste0(deDupCols, "-", leftName)
        deDupCols <- deDupCols[leftCols %in% names(consFormulaList[[grpi]])]
        leftCols <- leftCols[leftCols %in% names(consFormulaList[[grpi]])]
        if (length(leftCols) > 0)
            setnames(consFormulaList[[grpi]], leftCols, deDupCols)

        consFormulaList[[grpi]][, coverage := (sapply(mergedBy, countCharInStr, ",") + 1) / length(allFormulas)]

        if (formThreshold > 0)
            consFormulaList[[grpi]] <- consFormulaList[[grpi]][coverage >= formThreshold]

        consFormulaList[[grpi]] <- rankFormulaTable(consFormulaList[[grpi]])
    }

    cat("Done!\n")

    return(formulas(formulas = consFormulaList, featureFormulas = list(),
                    algorithm = paste0(unique(sapply(allFormulas, algorithm)), collapse = ",")))
})


#' @templateVar func generateFormulas
#' @templateVar what generate formulae
#' @templateVar ex1 generateFormulasDA
#' @templateVar ex2 generateFormulasGenForm
#' @templateVar algos bruker,genform,sirius
#' @template generic-algo
#'
#' @param ... Any parameters to be passed to the selected formula generation
#'   algorithm.
#'
#' @rdname formula-generation
#' @aliases generateFormulas
#' @export
setMethod("generateFormulas", "featureGroups", function(fGroups, algorithm, ...)
{
    f <- switch(algorithm,
                bruker = generateFormulasDA,
                genform = generateFormulasGenForm,
                sirius = generateFormulasSirius,
                stop("Invalid algorithm! Should be: bruker, genform or sirius"))

    f(fGroups, ...)
})
