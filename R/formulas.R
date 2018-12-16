#' @include main.R
NULL

#' Formula lists class
#'
#' Contains data of generated chemical formulae for given feature groups.
#'
#' \code{formulas} objects are obtained from \link[=formula-generation]{formula
#' generators}.
#'
#' @slot formulas Lists of all generated formulae. Use the \code{formulaTable}
#'   method for access.
#' @slot algorithm The algorithm that was used for generation of formulae. Use
#'   the \code{algorithm} method for access.
#'
#' @param obj,x,object,formulas The \code{formulas} object.
#'
#' @templateVar seli analyses
#' @templateVar selOrderi analyses()
#' @templateVar selj feature groups
#' @templateVar selOrderj groupNames()
#' @templateVar optionalji TRUE
#' @template sub_op-args
#'
#' @export
formulas <- setClass("formulas",
                     slots = c(formulas = "list", featureFormulas = "list", algorithm = "character"))

#' @describeIn formulas Accessor method to obtain generated formulae.
#' @return \code{formulas} returns a \code{list} containing for each analysis
#'   and each feature group a \code{\link{data.table}} with an overview of all
#'   generated formulae and other data such as candidate scoring and MS/MS
#'   fragments.
#' @export
setMethod("formulaTable", "formulas", function(obj) obj@formulas)

# UNDONE
setMethod("featureFormulas", "formulas", function(obj) obj@featureFormulas)

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
    printf("A formulas object (%s)\n", class(object))
    printf("Algorithm: %s\n", algorithm(object))

    ft <- featureFormulas(object)
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

    showObjectSize(object)
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
        x@featureFormulas <- sapply(x@eatureFormulas, function(a) pruneList(a[i]),
                                    simplify = FALSE)
        x@eatureFormulas <- pruneList(x@eatureFormulas, TRUE)

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
#' @param fGroups The \code{\link{featureGroups}} object that was used to
#'   generate this \code{formulas} object.
#' @param elements,fragElements A \code{character} vector with elements that
#'   should be counted for each MS(/MS) formula candidate. For instance,
#'   \code{c("C", "H")} adds columns for both carbon and hydrogen amounts in
#'   each formula. Set to \code{NULL} to not count any elements.
#' @param OM If set to \code{TRUE} several columns with information relevant for
#'   organic matter (OM) characterization will be added (e.g. elemental ratios,
#'   classification). This will also make sure that \code{elements} contains at
#'   least C, H, N, O, P and S.
#' @param maxFormulas,maxFragFormulas Maximum amount of unique candidate
#'   formulae (or fragment formulae) per feature group. Set to \code{NULL} to
#'   ignore.
#'
#' @return \code{makeTable} returns a \code{\link{data.table}}.
#'
#' @export
setMethod("makeTable", "formulas", function(obj, fGroups, average = FALSE, elements = NULL,
                                            fragElements = NULL, OM = FALSE,
                                            maxFormulas = NULL, maxFragFormulas = NULL)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    checkmate::assertFlag(average, add = ac)
    checkmate::assertCharacter(elements, min.chars = 1, any.missing = FALSE, null.ok = TRUE, add = ac)
    checkmate::assertCharacter(fragElements, min.chars = 1, any.missing = FALSE, null.ok = TRUE, add = ac)
    checkmate::assertFlag(OM, add = ac)
    checkmate::reportAssertions(ac)

    gInfo <- groupInfo(fGroups)

    ret <- rbindlist(formulaTable(obj), fill = TRUE, idcol = "group")
    ret[, c("ret", "mz") := gInfo[group, ]]
    setcolorder(ret, c("group", "ret", "mz"))

    if (average)
    {
        # UNDONE: continue this?

        # collapse byMSMS: will be TRUE if at least an MS/MS formula candidate was there
        ret[, byMSMS := any(byMSMS), by = "group"]

        ret[, formula_avg_count := length(unique(formula)), by = "group"]

        avgCols <- c("formula", "neutral_formula")
        ret[, (avgCols) := lapply(.SD, function(f) averageFormulas(unique(f))), .SDcols = avgCols, by = "group"]

        # remove columns which don't really make sense anymore
        rmCols <- c("neutral_loss", "error", "formula_mz", "dbe", "anaCoverage",
                    "adduct", "mSigma", "rank", "explainedPeaks", "explainedIntensity",
                    # add any fragment columns
                    grep("^frag_", names(ret), value = TRUE),
                    formulaScoringColumns())

        rmCols <- getAllFormulasCols(rmCols, names(ret))
        if (length(rmCols) > 0)
            ret[, (rmCols) := NULL]

        # average scores
        # scCols <- intersect(names(ret), formulaScoringColumns())
        # ret[, (scCols) := as.list(colMeans(.SD)), .SDcols = scCols, by = "group"]

        ret <- unique(ret, by = "formula")
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

    # ensure CHNOPS counts are present
    if (OM)
        elements <- unique(c(if (is.null(elements)) c() else elements, c("C", "H", "N", "O", "P", "S")))

    if (!is.null(elements) && length(elements) > 0)
    {
        # Retrieve element lists from formulas
        el <- getElements(ret$formula, elements)
        ret[, names(el) := el]
    }
    if (!is.null(fragElements) && !is.null(ret[["frag_formula"]]) &&
        length(fragElements) > 0)
    {
        el <- getElements(ret$frag_formula, fragElements)
        ret[, (paste0("frag_", names(el))) := el]
    }

    if (OM)
    {
        # add element ratios commonly used for plotting
        elrat <- function(el1, el2) ifelse(el2 == 0, 0, el1 / el2)
        ret[, c("OC", "HC", "NC") := .(elrat(O, C), elrat(H, C), elrat(N, C))]

        # aromaticity index and related DBE (see Koch 2016, 10.1002/rcm.7433)
        ret[, DBE_AI := 1 + C - O - S - 0.5 * (N + P + H)]
        getAI <- function(dbe, cai) ifelse(cai == 0, 0, dbe / cai)
        ret[, AI := getAI(DBE_AI, (C - O - N - S - P))]

        ret[, classification := Vectorize(classifyFormula)(OC, HC, NC, AI)]
    }

    return(ret)
})

# NOTE: feature formulas untouched
setMethod("filter", "formulas", function(obj, minExplainedMSMSPeaks = NULL, neutralLoss = NULL,
                                         topMost = NULL)
{
    # UNDONE: also add (ranges for) element counts/ratios, AI, scorings, errors, classification?
    # UNDONE: handle merged column names

    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertCount, . ~ topMost + minExplainedMSMSPeaks,
           positive = c(TRUE, TRUE, FALSE), null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertCharacter(neutralLoss, min.chars = 1, min.len = 1, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    cat("Filtering formulas... ")

    oldn <- length(obj)
    obj@formulas <- pruneList(sapply(obj@formulas, function(formTable)
    {
        if (!is.null(minExplainedMSMSPeaks) && minExplainedMSMSPeaks > 0)
        {
            fragCounts <- formTable[, ifelse(byMSMS, length(frag_formula), 0), by = "formula"]
            formTable[fragCounts >= minExplainedMSMSPeaks]
        }

        if (!is.null(neutralLoss))
            formTable <- formTable[byMSMS == TRUE & neutral_loss %in% neutralLoss]

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

#' @describeIn formulas Generates a consensus and other summarizing data for
#'   formulae of feature groups that are present in multiple analyses.
#'
#' @param \dots Any further \code{formulas} objects on top of the object given
#'   for the \code{obj} parameter to which a consensus should be generated.
#' @param formThreshold Fractional minimum amount (0-1) to
#'   which a generated precursor formula should be present in all given formula objects
#'   (\code{formListThreshold}) containing the related feature group. For
#'   instance, a value of \samp{0.5} for \code{formAnaThreshold} means that a
#'   particular formula should be present in at least \samp{50\%} of all
#'   analyses containing the feature group that was used to generate the
#'   formula.
#'
#' @return \code{consensus} returns a \code{formulas} object that is produced
#'   by merging multiple \code{formulas} objects.
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

    allFormulas <- allFormulas[lengths(allFormulas) > 0]
    if (length(allFormulas) < 2)
        stop("Need at least two non-empty formulas objects")

    allFormNames <- sapply(allFormulas, algorithm)
    allFormNames <- make.unique(allFormNames)

    allFormulasLists <- sapply(seq_along(allFormulas), function(fi)
    {
        return(lapply(formulaTable(allFormulas[[fi]]), function(ft)
        {
            ret <- copy(ft)
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
                unCols <- c(uniqueCols, c("formula", "byMSMS", "frag_formula"))
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
            }

            consFormulaList[[grp]] <- mTable
        }

        cat("Done!\n")
    }

    printf("Determining coverage... ")

    for (grpi in seq_along(consFormulaList))
    {
        # fix up de-duplicated column names
        deDupCols <- c(uniqueCols, "formula", "byMSMS", "frag_formula")
        leftCols <- paste0(deDupCols, "-", leftName)
        deDupCols <- deDupCols[leftCols %in% names(consFormulaList[[grpi]])]
        leftCols <- leftCols[leftCols %in% names(consFormulaList[[grpi]])]
        if (length(leftCols) > 0)
            setnames(consFormulaList[[grpi]], leftCols, deDupCols)

        # match all that has a dash inbetween
        mergedCols <- getAllMergedFormulasCols(names(consFormulaList[[grpi]]))

        # figure out what was merged (i.e. name after dash)
        mergedColsWhat <- sub(".+\\-", "", mergedCols)

        if (length(unique(mergedColsWhat)) < 2) # nothing was merged
            consFormulaList[[grpi]][, coverage := 1 / length(allFormulas)]
        else
        {
            for (r in seq_len(nrow(consFormulaList[[grpi]])))
            {
                # count how many merged objects have a value in this line: each
                # object should at least have one non-NA value
                notNAs <- sapply(mergedCols, function(mc) !is.na(consFormulaList[[grpi]][[mc]][r]))
                mergedNotNACount <- length(unique(mergedColsWhat[notNAs]))

                set(consFormulaList[[grpi]], r, "coverage", mergedNotNACount / length(allFormulas))
            }
        }

        if (formThreshold > 0)
            consFormulaList[[grpi]] <- consFormulaList[[grpi]][coverage >= formThreshold]

        # setcolorder(consFormulaList[[grpi]], formConsensusColOrder(consFormulaList[[grpi]]))
        consFormulaList[[grpi]] <- rankFormulaTable(consFormulaList[[grpi]])
    }

    cat("Done!\n")

    return(formulas(formulas = consFormulaList, featureFormulas = list(),
                    algorithm = paste0(unique(sapply(allFormulas, algorithm)), collapse = ",")))
})


#' Formula consensus class
#'
#' Contains all formulae for given feature groups after a consensus was made.
#'
#' Objects of \code{formulaConsensus} are generated by the
#' \code{\link{consensus}} method.
#'
#' @slot formulas A \code{\link{data.table}} containing an overview of all
#'   consensus formulae. Use the \code{formulaTable} method for access.
#' @slot algorithm The algorithm that was used for generation of formulae. Use
#'   the \code{algorithm} method for access.
#'
#' @param obj,x,object The \code{formulaConsensus} object.
#'
#' @templateVar seli feature groups
#' @templateVar selOrderi groupNames()
#' @templateVar dollarOpName feature group
#' @template sub_op-args
#'
#' @seealso \code{\link{formulas}}
#'
#' @export
formulaConsensus <- setClass("formulaConsensus",
                             slots = c(formulas = "data.table", algorithm = "character"),
                             prototype = list(formulas = data.table(), algorithm = "none"))

#' @describeIn formulaConsensus Accessor method to obtain generated formulae.
#' @return \code{formulaTable} returns a \code{\link{data.table}} containing all
#'   generated formulae, scoring information, MS/MS fragments, etc.
#' @export
setMethod("formulaTable", "formulaConsensus", function(obj) obj@formulas)

#' @describeIn formulaConsensus Accessor method for the algorithm (a character
#'   string) used to generate formulae.
#' @export
setMethod("algorithm", "formulaConsensus", function(obj) obj@algorithm)

#' @templateVar class formulaConsensus
#' @templateVar what feature groups
#' @template strmethod
#' @export
setMethod("groupNames", "formulaConsensus", function(obj) unique(obj@formulas$group))

#' @describeIn formulaConsensus Obtain total number of formulae entries.
#' @export
setMethod("length", "formulaConsensus", function(x) nrow(x@formulas))

#' @describeIn formulaConsensus Show summary information for this object.
#' @export
setMethod("show", "formulaConsensus", function(object)
{
    printf("A formulaConsensus object (%s)\n", class(object))
    printf("Algorithm: %s\n", paste0(algorithm(object), collapse = "/"))

    printf("Total formula count: %d (all), %d (MS) and %d (MS/MS)\n", length(object),
           nrow(object@formulas[byMSMS == FALSE]), nrow(object@formulas[byMSMS == TRUE]))

    fCounts <- object@formulas[, .(byMSMS, .N), by = group]
    printf("Average formulas per feature group: %.1f (all), %.1f (MS), %.1f (MS/MS)\n", mean(fCounts$N),
           mean(fCounts[byMSMS == FALSE, N]), mean(fCounts[byMSMS == TRUE, N]))

    showObjectSize(object)
})

#' @describeIn formulaConsensus Subset on feature groups.
#' @export
setMethod("[", c("formulaConsensus", "ANY", "missing", "missing"), function(x, i, ...)
{
    if (!missing(i))
        assertSubsetArg(i)

    if (!missing(i))
    {
        if (!is.character(i))
            i <- groupNames(x)[i]
        x@formulas <- x@formulas[group %in% i]
    }

    return(x)
})

#' @describeIn formulaConsensus Extract formula information a feature group.
#' @export
setMethod("[[", c("formulaConsensus", "ANY", "missing"), function(x, i, j)
{
    assertExtractArg(i)
    if (!is.character(i))
        i <- groupNames(x)[i]
    return(x@formulas[group == i])
})

#' @describeIn formulaConsensus Extract formula information a feature group.
#' @export
setMethod("$", "formulaConsensus", function(x, name)
{
    eval(substitute(x@formulas[group == NAME_ARG], list(NAME_ARG = name)))
})

#' @describeIn formulaConsensus Plots an annotated spectrum for a given
#'   candidate formula of a feature group.
#'
#' @param precursor The formula of the precursor (in ionic form, \emph{i.e.} as
#'   detected by the MS).
#'
#' @template plotSpec-args
#'
#' @template useGGplot2
#'
#' @return \code{plotSpec} will return a \code{\link[=ggplot2]{ggplot object}}
#'   if \code{useGGPlot2} is \code{TRUE}.
#'
#' @export
setMethod("plotSpec", "formulaConsensus", function(obj, precursor, groupName, MSPeakLists, useGGPlot2 = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertString(precursor, min.chars = 1, add = ac)
    checkmate::assertString(groupName, min.chars = 1, add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertFlag(useGGPlot2, add = ac)
    checkmate::reportAssertions(ac)

    formTable <- formulaTable(obj)[group == groupName & byMSMS == TRUE & formula == precursor]

    if (nrow(formTable) == 0)
        return(NULL)

    spec <- MSPeakLists[[groupName]][["MSMS"]]
    if (is.null(spec))
        return(NULL)

    fi <- getFragmentInfoFromForms(spec, formTable)

    if (useGGPlot2)
        return(makeMSPlotGG(spec, fi) + ggtitle(precursor))

    makeMSPlot(spec, fi, main = precursor)
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
