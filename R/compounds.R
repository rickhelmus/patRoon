#' @include main.R
#' @include workflow-step.R
NULL

#' Compound lists class
#'
#' Contains data of generated chemical compounds for given feature groups.
#'
#' \code{compounds} objects are obtained from
#' \link[=compound-generation]{compound generators}.
#'
#' @slot compounds Lists of all generated compounds. Use the \code{compounds}
#'   method for access.
#' @slot scoreTypes A \code{character} with all the score types that were used
#'   when generating the compounds.
#' @slot scoreRanges The original min/max values of all scorings when candidate
#'   results were generated. This is used for normalization.
#'
#' @param formulas The \code{\link{formulas}} object that should be used for
#'   scoring/annotation. For \code{plotSpec}: set to \code{NULL} to ignore.
#' @param obj,object,x,compounds The \code{compound} object.
#' @param index The numeric index of the candidate structure. Multiple indices
#'   (\emph{i.e.} vector with length >=2) may be specified for
#'   \code{plotStructure} and are mandatory for \code{getMCS}. Alternatively,
#'   \samp{-1} may be specified to these methods to select all candidates. When
#'   multiple indices are specified for \code{plotStructure}, their maximum
#'   common substructure will be drawn.
#' @param \dots For \code{plotSpec}: Further arguments passed to
#'   \code{\link[graphics]{plot}}.
#'
#'   Others: Any further (and unique) \code{compounds} objects.
#'
#' @templateVar seli feature groups
#' @templateVar selOrderi groupNames()
#' @templateVar dollarOpName feature group
#' @template sub_op-args
#'
#' @return \code{plotSpec} and \code{plotStructure} will return a
#'   \code{\link[=ggplot2]{ggplot object}} if \code{useGGPlot2} is \code{TRUE}.
#'
#' @template plotSpec-args
#'
#' @template useGGplot2
#'
#' @templateVar normParam normalizeScores
#' @templateVar excludeParam excludeNormScores
#' @template norm-args
#'
#' @templateVar class compounds
#' @template class-hierarchy
#'
#' @export
compounds <- setClass("compounds",
                      slots = c(compounds = "list", scoreTypes = "character", scoreRanges = "list"),
                      contains = "workflowStep")

setMethod("initialize", "compounds", function(.Object, ...)
{
    .Object <- callNextMethod(.Object, ...)
    .Object@compounds <- makeEmptyListNamed(.Object@compounds)
    return(.Object)
})

#' @rdname compounds-class
compoundsConsensus <- setClass("compoundsConsensus",
                               slots = c(mergedCompNames = "character"),
                               contains = "compounds")


#' @describeIn compounds Accessor method to obtain generated compounds.
#' @return \code{compoundTable} returns a \code{list} containing for each feature
#'   group a \code{\link{data.table}} with an overview of all candidate
#'   compounds and other data such as candidate scoring, matched MS/MS
#'   fragments, etc.
#' @aliases compoundTable
#' @export
setMethod("compoundTable", "compounds", function(obj) obj@compounds)

#' @describeIn compounds Accessor method for the algorithm (a character
#'   string) used to generate compounds.
#' @export
setMethod("algorithm", "compounds", function(obj) obj@algorithm)

#' @templateVar class compounds
#' @templateVar what feature groups
#' @template strmethod
#' @export
setMethod("groupNames", "compounds", function(obj) names(obj@compounds))

#' @describeIn compounds Obtain total number of candidate compounds.
#' @export
setMethod("length", "compounds", function(x) if (length(x@compounds) > 0) sum(sapply(x@compounds, nrow)) else 0)

#' @describeIn compounds Show summary information for this object.
#' @export
setMethod("show", "compounds", function(object)
{
    callNextMethod()

    mn <- mergedCompoundNames(object)
    if (length(mn) > 1)
        printf("Merged: %s\n", paste0(mn, collapse = ", "))

    printf("Number of feature groups with compounds in this object: %d\n", length(object@compounds))

    cCounts <- if (length(object) == 0) 0 else sapply(object@compounds, nrow)
    printf("Number of compounds: %d (total), %.1f (mean), %d - %d (min - max)\n",
           sum(cCounts), mean(cCounts), min(cCounts), max(cCounts))
})

#' @describeIn compounds Subset on feature groups.
#' @export
setMethod("[", c("compounds", "ANY", "missing", "missing"), function(x, i, ...)
{
    if (!missing(i))
    {
        i <- assertSubsetArgAndToChr(i, groupNames(x))
        x@compounds <- x@compounds[i]
        x@scoreRanges <- x@scoreRanges[i]
    }

    return(x)
})

#' @describeIn compounds Extract a compound table for a feature group.
#' @export
setMethod("[[", c("compounds", "ANY", "missing"), function(x, i, j)
{
    assertExtractArg(i)
    return(x@compounds[[i]])
})

#' @describeIn compounds Extract a compound table for a feature group.
#' @export
setMethod("$", "compounds", function(x, name)
{
    eval(substitute(x@compounds$NAME_ARG, list(NAME_ARG = name)))
})

#' @describeIn compounds Returns all MS peak list data in a table.
#'
#' @param fragments If \code{TRUE} then information on annotated fragments will
#'   be included.
#'
#' @template as_data_table-args
#'
#' @export
setMethod("as.data.table", "compounds", function(x, fGroups = NULL, fragments = FALSE, normalizeScores = "none",
                                                 excludeNormScores = c("score", "individualMoNAScore"))
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", null.ok = TRUE, add = ac)
    checkmate::assertFlag(fragments, add = ac)
    assertNormalizationMethod(normalizeScores, add = ac)
    checkmate::assertCharacter(excludeNormScores, min.chars = 1, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    mcn <- mergedCompoundNames(x)
    cTable <- compoundTable(x)
    if (normalizeScores != "none")
    {
        cTable <- mapply(cTable, x@scoreRanges, SIMPLIFY = FALSE, FUN = normalizeCompScores,
                         MoreArgs = list(mcn, normalizeScores == "minmax", excludeNormScores))
    }

    if (fragments)
    {
        ret <- rbindlist(lapply(cTable, function(ct)
        {
            ct <- copy(ct)
            ct[, row := seq_len(nrow(ct))]

            fragTab <- rbindlist(ct$fragInfo, idcol = "row", fill = TRUE)
            fragTab[, PLIndex := NULL]
            cnames <- setdiff(names(fragTab), "row")
            setnames(fragTab, cnames, paste0("frag_", cnames))

            return(merge(ct, fragTab, by = "row")[, -"row"])
        }), idcol = "group", fill = TRUE)
    }
    else
        ret <- rbindlist(cTable, idcol = "group", fill = TRUE)

    if (!is.null(fGroups))
    {
        ret[, c("ret", "group_mz") := groupInfo(fGroups)[group, c("rts", "mzs")]]
        setcolorder(ret, c("group", "ret", "group_mz"))
    }

    if (!is.null(ret[["fragInfo"]]))
        return(ret[, -"fragInfo"]) # not there if empty results
    return(ret)
})

#' @describeIn compounds Returns a list containing for each feature group a
#'   character vector with database identifiers for all candidate compounds. The
#'   list is named by feature group names, and is typically used with the
#'   \code{identifiers} option of \code{\link{generateCompoundsMetfrag}}.
#' @aliases identifiers
#' @export
setMethod("identifiers", "compounds", function(compounds)
{
    compTable <- compoundTable(compounds)
    return(sapply(names(compTable),
                  function(grp) unlist(strsplit(as.character(compTable[[grp]]$identifier), ";")), simplify = FALSE))
})

#' @describeIn compounds Provides rule based filtering for generated compounds.
#'   Useful to eliminate unlikely candidates and speed up further processing.
#'
#' @param minExplainedPeaks,minScore,minFragScore,minFormulaScore Minimum number
#'   of explained peaks, overall score, in-silico fragmentation score and
#'   formula score, respectively. Set to \code{NULL} to ignore. The
#'   \code{scoreLimits} argument allows for more advanced score filtering.
#' @param scoreLimits Filter results by their scores. Should be a named
#'   \code{list} that contains two-sized numeric vectors with the
#'   minimum/maximum value of a score (use \code{-Inf}/\code{Inf} for no
#'   limits). The names of each element should follow the values returned by
#'   \code{\link{compoundScorings}()$name}. For instance,
#'   \code{scoreLimits=list(numberPatents=c(10, Inf))} specifies that
#'   \code{numberPatents} should be at least \samp{10}. For more details of
#'   scorings see \code{\link{compoundScorings}}. Note that a result without a
#'   specified scoring is never removed. Set to \code{NULL} to skip this filter.
#' @param topMost Only keep a maximum of \code{topMost} candidates with highest
#'   score (or least highest if \code{negate=TRUE}). Set to \code{NULL} to ignore.
#' @param negate If \code{TRUE} then filters are applied in opposite manner.
#'
#' @template element-args
#'
#' @return \code{filter} returns a filtered \code{compounds} object.
#'
#' @export
setMethod("filter", "compounds", function(obj, minExplainedPeaks = NULL, minScore = NULL, minFragScore = NULL,
                                          minFormulaScore = NULL, scoreLimits = NULL, elements = NULL,
                                          fragElements = NULL, lossElements = NULL, topMost = NULL,
                                          negate = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertCount, . ~ minExplainedPeaks + topMost, positive = c(FALSE, TRUE),
           null.ok = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ minScore + minFragScore + minFormulaScore, finite = TRUE,
           null.ok = TRUE, fixed = list(add = ac)) # note: negative scores allowed for SIRIUS
    checkmate::assertList(scoreLimits, null.ok = TRUE, types = "numeric", add = ac)
    if (!is.null(scoreLimits))
    {
        scCols <- unique(compoundScorings()$name)
        checkmate::assertNames(names(scoreLimits), type = "unique", subset.of = scCols, add = ac)
        checkmate::qassertr(scoreLimits, "N2")
    }
    aapply(checkmate::assertCharacter, . ~ elements + fragElements + lossElements,
           min.chars = 1, min.len = 1, null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertFlag(negate, add = ac)
    checkmate::reportAssertions(ac)

    cat("Filtering compounds... ")

    mCompNames <- mergedCompoundNames(obj)
    filterMinCols <- function(cmpTable, col, minVal)
    {
        cols <- getAllCompCols(col, names(cmpTable), mCompNames)
        pred <- function(cl) is.na(cl) | cl >= minVal
        if (negate)
            pred <- Negate(pred)
        for (cl in cols)
            cmpTable <- cmpTable[pred(get(cl))]
        return(cmpTable)
    }

    fMinVals <- c(explainedPeaks = minExplainedPeaks, score = minScore, fragScore = minFragScore,
                  formulaScore = minFormulaScore)
    fMinVals <- fMinVals[!sapply(fMinVals, is.null)]

    oldn <- length(obj)
    obj@compounds <- sapply(obj@compounds, function(cmpTable)
    {
        for (cl in names(fMinVals))
            cmpTable <- filterMinCols(cmpTable, cl, fMinVals[cl])

        if (!is.null(scoreLimits))
        {
            for (sc in names(scoreLimits))
            {
                cols <- getAllCompCols(sc, names(cmpTable), mCompNames)
                if (length(cols) == 0)
                    next

                keep <- cmpTable[, do.call(pmin, c(.SD, list(na.rm = TRUE))) >= scoreLimits[[sc]][1] &
                                     do.call(pmax, c(.SD, list(na.rm = TRUE))) <= scoreLimits[[sc]][2],
                                 .SDcols = cols]
                if (negate)
                    keep <- !keep
                cmpTable <- cmpTable[keep]
            }
        }

        if (nrow(cmpTable) == 0)
            return(cmpTable)

        if (!is.null(elements))
        {
            keep <- sapply(cmpTable$formula, checkFormula, elements, negate = negate)
            cmpTable <- cmpTable[keep]
        }
        if (!is.null(fragElements) || !is.null(lossElements))
        {
            keep <- sapply(cmpTable$fragInfo, function(fi)
            {
                if (nrow(fi) == 0)
                    return(FALSE)
                if (!is.null(fragElements) && !any(sapply(fi$formula, checkFormula, fragElements, negate = negate)))
                    return(FALSE)
                if (!is.null(lossElements) && !any(sapply(fi$neutral_loss, checkFormula, lossElements, negate = negate)))
                    return(FALSE)
                return(TRUE)
            })
            cmpTable <- cmpTable[keep]
        }

        if (!is.null(topMost))
        {
            if (negate)
                cmpTable <- tail(cmpTable, topMost)
            else
                cmpTable <- head(cmpTable, topMost)
        }

        return(cmpTable)
    }, simplify = FALSE)

    if (length(obj) > 0)
        obj@compounds <- obj@compounds[sapply(obj@compounds, function(cm) !is.null(cm) && nrow(cm) > 0)]

    obj@scoreRanges <- obj@scoreRanges[names(obj@compounds)]
    
    newn <- length(obj)
    printf("Done! Filtered %d (%.2f%%) compounds. Remaining: %d\n", oldn - newn, if (oldn == 0) 0 else (1-(newn/oldn))*100, newn)
    return(obj)
})

#' @describeIn compounds Adds formula ranking data from a \code{\link{formulas}}
#'   object as an extra compound candidate scoring (\code{formulaScore} column).
#'   The formula score for each compound candidate is between \samp{0-1}, where
#'   \emph{zero} means no match with any formula candidates, and \emph{one}
#'   means that the compound candidate's formula is the highest ranked.
#'
#' @param updateScore If set to \code{TRUE} then the \code{score} column is
#'   updated by adding the normalized \option{formulaScore} (weighted by
#'   \option{formulaScoreWeight}). Currently, this \strong{only} makes sense for
#'   \command{MetFrag} results!
#' @param formulaScoreWeight Weight used to update scoring (see
#'   \code{updateScore} parameter).
#'
#' @return \code{addFormulaScoring} returns a \code{compounds} object updated
#'   with formula scoring.
#'
#' @aliases addFormulaScoring
#' @export
setMethod("addFormulaScoring", "compounds", function(compounds, formulas, updateScore,
                                                     formulaScoreWeight)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(formulas, "formulas", add = ac)
    checkmate::assertFlag(updateScore, add = ac)
    checkmate::assertNumber(formulaScoreWeight, lower = 0, finite = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    if (length(compounds) == 0)
        return(compounds)

    cTable <- compoundTable(compounds)
    cGNames <- names(cTable)

    # UNDONE?
    if (length(mergedCompoundNames(compounds)) > 0)
        stop("Currently formula scoring cannot be calculated for consensus results. Please add the scorings before calling consensus()")
    
    calculateScores <- function(cr, forms)
    {
        forms <- unique(forms, by = "neutral_formula") # ensure we have only one row per precursor (ranking remains)
        
        fCount <- nrow(forms)
        fRanks <- match(cr$formula, forms$neutral_formula, nomatch = fCount + 1) # add one to no-match results to make it the worst score
        fRanks <- (fCount - (fRanks - 1)) / fCount # convert to 0-1 where one is best ranked and zero a no match
        
        return(fRanks)
    }

    printf("Adding formula scoring...\n")
    prog <- openProgBar(0, length(cTable))

    cTable <- lapply(seq_along(cTable), function(grpi)
    {
        forms <- formulas[[cGNames[grpi]]]

        ct <- copy(cTable[[grpi]])

        if (nrow(ct) == 0)
            return(ct)

        if (length(forms) == 0 || nrow(forms) == 0)
            ct[, formulaScore := 0]
        else
            ct[, formulaScore := calculateScores(ct, forms)]

        # update overall scoring
        if (any(ct$formulaScore > 0))
        {
            normFormScores <- ct$formulaScore / max(ct$formulaScore)
            ct[, score := score + formulaScoreWeight * normFormScores]
        }

        setTxtProgressBar(prog, grpi)
        return(ct)
    })
    names(cTable) <- cGNames

    setTxtProgressBar(prog, length(cTable))
    close(prog)

    compounds@compounds <- cTable
    compounds@scoreRanges <- mapply(compounds@scoreRanges, cTable, SIMPLIFY = FALSE, FUN = function(sc, ct)
    {
        ret <- c(sc, list(formulaScore = range(ct$formulaScore)))
        # extend score range if necessary
        ret$score <- c(min(ret$score, ct$score, na.rm = TRUE),
                       max(ret$score, ct$score, na.rm = TRUE))
        return(ret)
    })
    compounds@scoreTypes <- union(compounds@scoreTypes, "formulaScore")

    return(compounds)
})

#' @describeIn compounds Calculates the maximum common substructure (MCS)
#'   for two or more candidate structures for a feature group. This method uses
#'   the \code{\link{get.mcs}} function from \CRANpkg{rcdk}.
#' @return \code{getMCS} returns an \CRANpkg{rcdk} molecule object
#'   (\code{IAtomContainer}).
#' @export
setMethod("getMCS", "compounds", function(obj, index, groupName)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assert(
        checkmate::checkIntegerish(index, lower = 1, any.missing = FALSE, min.len = 2, unique = TRUE),
        checkmate::checkTRUE(index == -1),
        .var.name = "index"
    )

    assertChoiceSilent(groupName, names(obj@compounds), add = ac)
    checkmate::reportAssertions(ac)

    if (length(index) == 1 && index == -1)
        index <- seq_len(nrow(compoundTable(obj)[[groupName]]))

    mols <- getMoleculesFromSMILES(compoundTable(obj)[[groupName]][["SMILES"]][index])
    mcons <- mols[[1]]
    if (length(mols) > 1)
    {
        for (i in seq(2, length(mols)))
        {
            if (!isValidMol(mols[[i]]))
                return(emptyMol())

            # might fail if there is no overlap...
            tryCatch(mcons <- rcdk::get.mcs(mcons, mols[[i]]), error = function(e) FALSE)
            if (mcons == FALSE)
                return(emptyMol())
        }
    }

    return(mcons)
})

# NOTE: argument docs 'borrowed' from plotSpec-args.R template

#' @describeIn compounds Plots a structure of a candidate compound using the
#'   \CRANpkg{rcdk} package. If multiple candidates are specified (\emph{i.e.}
#'   by specifying a \code{vector} for \code{index}) then the maximum common
#'   substructure (MCS) of the selected candidates is drawn.
#'
#' @param width,height The dimensions (in pixels) of the raster image that
#'   should be plotted.
#'
#' @references \addCitations{rcdk}{1}
#'
#' @export
setMethod("plotStructure", "compounds", function(obj, index, groupName, width = 500, height = 500, useGGPlot2 = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertIntegerish(index, lower = -1, any.missing = FALSE, min.len = 1, unique = TRUE, add = ac)
    checkmate::assertString(groupName, min.chars = 1, add = ac)
    aapply(checkmate::assertNumber, . ~ width + height, lower = 0, finite = TRUE, fixed = list(add = ac))
    checkmate::assertFlag(useGGPlot2, add = ac)
    checkmate::reportAssertions(ac)

    compTable <- compoundTable(obj)[[groupName]]

    if (is.null(compTable) || nrow(compTable) == 0)
        return(NULL)

    if (length(index) > 1 || index == -1)
        mol <- getMCS(obj, index, groupName)
    else
        mol <- getMoleculesFromSMILES(compTable$SMILES[index], emptyIfFails = TRUE)[[1]]

    img <- getRCDKStructurePlot(mol, width, height)
    if (useGGPlot2)
        cowplot::ggdraw() + cowplot::draw_image(img)
    else
        plot(img)
})

setMethod("plotStructureHash", "compounds", function(obj, index, groupName, width = 500,
                                                     height = 500, useGGPlot2 = FALSE)
{
    compTable <- compoundTable(obj)[[groupName]]
    SMI <- if (is.null(compTable) || nrow(compTable) == 0) NULL else compTable$SMILES[index]
    return(makeHash(SMI, width, height, useGGPlot2))
})

#' @describeIn compounds Plots a barplot with scoring of a candidate compound.
#'
#' @param onlyUsed If \code{TRUE} then only scorings are plotted that actually
#'   have been used to rank data (see the \code{scoreTypes} argument to
#'   \code{\link{generateCompoundsMetfrag}} for more details).
#'
#' @export
setMethod("plotScores", "compounds", function(obj, index, groupName, normalizeScores = "max",
                                              excludeNormScores = c("score", "individualMoNAScore"),
                                              onlyUsed = TRUE, useGGPlot2 = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertCount(index, positive = TRUE, add = ac)
    checkmate::assertString(groupName, min.chars = 1, add = ac)
    checkmate::assertChoice(normalizeScores, c("none", "max", "minmax"))
    checkmate::assertCharacter(excludeNormScores, min.chars = 1, null.ok = TRUE, add = ac)
    checkmate::assertFlag(onlyUsed, add = ac)
    checkmate::assertFlag(useGGPlot2, add = ac)
    checkmate::reportAssertions(ac)

    compTable <- compoundTable(obj)[[groupName]]

    if (is.null(compTable) || nrow(compTable) == 0 || index > nrow(compTable))
        return(NULL)

    mcn <- mergedCompoundNames(obj)

    if (normalizeScores != "none")
        compTable <- normalizeCompScores(compTable, obj@scoreRanges[[groupName]], mcn, normalizeScores == "minmax", excludeNormScores)

    scoreCols <- getAllCompCols(c(getCompScoreColNames(), getCompSuspectListColNames()), names(compTable), mcn)
    if (onlyUsed)
        scoreCols <- intersect(scoreCols, obj@scoreTypes)

    makeScoresPlot(compTable[index, scoreCols, with = FALSE], mcn, useGGPlot2)
})

setMethod("plotScoresHash", "compounds", function(obj, index, groupName, normalizeScores = "max",
                                                  excludeNormScores = c("score", "individualMoNAScore"),
                                                  onlyUsed = TRUE, useGGPlot2 = FALSE)
{
    compTable <- compoundTable(obj)[[groupName]]
    if (is.null(compTable) || nrow(compTable) == 0 || index > nrow(compTable))
        compTable <- NULL
    else if (normalizeScores == "none")
        compTable <- compTable[index]
    
    return(makeHash(index, compTable, normalizeScores, excludeNormScores, onlyUsed, useGGPlot2))
})

#' @describeIn compounds Returns an MS/MS peak list annotated with data from a
#'   given candidate compound for a feature group.
#'
#' @param onlyAnnotated Set to \code{TRUE} to filter out any peaks that could
#'   not be annotated.
#'
#' @export
setMethod("annotatedPeakList", "compounds", function(obj, index, groupName, MSPeakLists, formulas = NULL,
                                                     onlyAnnotated = FALSE)
{
    checkmate::assertClass(formulas, "formulas", null.ok = TRUE)

    allFGroups <- groupNames(obj)
    if (!is.null(formulas))
        allFGroups <- union(allFGroups, groupNames(formulas))

    ac <- checkmate::makeAssertCollection()
    checkmate::assertCount(index, positive = TRUE, add = ac)
    assertChoiceSilent(groupName, allFGroups, add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertFlag(onlyAnnotated, add = ac)
    checkmate::reportAssertions(ac)

    spec <- MSPeakLists[[groupName]][["MSMS"]]
    if (is.null(spec))
        return(NULL)

    spec <- copy(spec)
    spec[, PLIndex := seq_len(nrow(spec))] # for merging

    compTable <- compoundTable(obj)[[groupName]]
    fragInfo <- NULL
    if (!is.null(compTable) && nrow(compTable) > 0)
    {
        compr <- compTable[index, ]
        fragInfo <- compr$fragInfo[[1]]
    }

    formTable <- if (!is.null(formulas)) formulas[[groupName]] else NULL
    if (!is.null(formTable))
    {
        formTable <- formTable[byMSMS == TRUE & neutral_formula == compr$formula]
        if (nrow(formTable) > 0)
        {
            formFragInfo <- getFragmentInfoFromForms(spec, formTable)
            if (is.null(fragInfo))
            {
                fragInfo <- formFragInfo
                fragInfo[, mergedBy := algorithm(formulas)]
            }
            else
                fragInfo <- mergeFragInfo(fragInfo, formFragInfo, algorithm(obj), algorithm(formulas))
        }
    }

    if (!is.null(fragInfo) && nrow(fragInfo) > 0)
        spec <- merge(spec, fragInfo[, -c("intensity", "mz")], all.x = TRUE, by = "PLIndex")

    spec <- spec[, PLIndex := NULL]

    if (onlyAnnotated)
    {
        if (is.null(spec[["formula"]]))
            spec <- spec[0]
        else
            spec <- spec[!is.na(formula)]
    }

    return(spec[])
})

#' @describeIn compounds Plots an annotated spectrum for a given candidate
#'   compound for a feature group.
#'
#' @param plotStruct If \code{TRUE} then the candidate structure is drawn in the
#'   spectrum.
#' @param title The title of the plot. If \code{NULL} a title will be
#'   automatically made.
#'
#' @template fsubscript_source
#'
#' @template plot-lim
#'
#' @export
setMethod("plotSpec", "compounds", function(obj, index, groupName, MSPeakLists, formulas = NULL,
                                            plotStruct = TRUE, title = NULL, useGGPlot2 = FALSE, xlim = NULL,
                                            ylim = NULL, ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertCount(index, positive = TRUE, add = ac)
    checkmate::assertString(groupName, min.chars = 1, add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertClass(formulas, "formulas", null.ok = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~plotStruct + useGGPlot2, fixed = list(add = ac))
    assertXYLim(xlim, ylim, add = ac)
    checkmate::reportAssertions(ac)

    compTable <- compoundTable(obj)[[groupName]]

    if (is.null(compTable) || nrow(compTable) == 0)
        return(NULL)

    if (!is.null(formulas) && groupName %in% groupNames(formulas))
        fTable <- formulas[[groupName]][byMSMS == TRUE]
    else
        fTable <- NULL

    compr <- compTable[index, ]
    spec <- annotatedPeakList(obj, index, groupName, MSPeakLists, formulas)

    if (plotStruct)
        mol <- getMoleculesFromSMILES(compr$SMILES)

    if (is.null(title))
    {
        if (!is.null(compr$compoundName) && !is.na(compr$compoundName) && nzchar(compr$compoundName))
            title <- subscriptFormula(compr$formula, over = compr$compoundName) #subscriptFormula(compr$formula, prefix = paste0(compr$compoundName, "\n("), postfix = ")")
        else
            title <- subscriptFormula(compr$formula)
    }

    if (!useGGPlot2)
    {
        # oldp <- par(mar = par("mar") * c(1, 1, 0, 0))

        if (plotStruct && isValidMol(mol))
        {
            molHInch <- 1.5
            makeMSPlot(spec, xlim, ylim, main = title, ..., extraHeightInch = molHInch)
        }
        else
            makeMSPlot(spec, xlim, ylim, main = title, ...)

        # draw structure
        if (plotStruct && isValidMol(mol))
        {
            img <- getRCDKStructurePlot(mol[[1]], 100, 100)

            dpi <- (par("cra")/par("cin"))[1]

            startx <- par("usr")[1]
            xlim <- par("usr")[2]
            ylim <- par("usr")[4]
            imgInfo <- magick::image_info(img)

            imgPlotW <- xinch(imgInfo$width / dpi)
            imgPlotH <- yinch(imgInfo$height / dpi)

            maxW <- 0.2 * xlim
            if (imgPlotW > maxW)
            {
                hresize <- imgPlotW / maxW
                imgPlotH <- imgPlotH / hresize
                imgPlotW <- maxW
            }

            maxH <- yinch(molHInch)
            if (imgPlotH > maxH)
            {
                wresize <- imgPlotH / maxH
                imgPlotW <- imgPlotW / wresize
                imgPlotH <- maxH
            }

            # offset a little
            startx <- startx + 0.01 * xlim
            ylim <- ylim * 0.99

            rasterImage(img, startx, ylim - imgPlotH, startx + imgPlotW, ylim)
        }

        # par(oldp)
    }
    else
    {
        MSPlot <- makeMSPlotGG(spec) + ggtitle(title)

        if (plotStruct && isValidMol(mol))
        {
            img <- getRCDKStructurePlot(mol[[1]], 100, 100, transparent = FALSE)

            # positioning if legend is on top... doesn't work too well :(
            # mar <- 7
            # pos <- grid::convertUnit(grid::unit(0.8, "npc") - grid::unit(mar, "pt"), "npc", valueOnly = TRUE)
            # MSPlot <- MSPlot + theme(plot.margin = margin(mar, mar, mar, mar)) +
            #     ylim(0, max(fi$intensity) * (1.2 + (0.8 - pos)))

            pos <- 0.8; size <- 1 - pos
            MSPlot <- MSPlot +
                ylim(0, max(spec$intensity) * (1 + size + 0.1)) # add a bit of space between most intense point+label
            MSPlot <- cowplot::ggdraw(MSPlot) +
                cowplot::draw_image(img, pos, pos, size, size)
        }

        return(MSPlot)
        # mcn <- mergedCompoundNames(obj)
        # scoreCols <- getAllCompCols(getCompScoreColNames(), names(compTable), mcn)
        # scores <- setnames(transpose(compr[, scoreCols, with = FALSE]), "score")
        # scores[, type := scoreCols]
        # scores <- scores[!is.na(score)]
        #
        # if (length(mcn) > 1)
        # {
        #     scores[, merged := "both"]
        #     for (n in mcn)
        #     {
        #         withM <- which(grepl(paste0("-", n), scores[["type"]], fixed = TRUE))
        #         set(scores, withM, "merged", n)
        #         set(scores, withM, "type", gsub(paste0("-", n), "", scores[["type"]][withM]))
        #     }
        # }

        # scorePlot <- ggplot(scores, aes_string(x = "type", y = "score")) +
        #     theme_cowplot(font_size = 12) +
        #     theme(axis.title.y = element_blank(), axis.title.x = element_blank(), # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        #           legend.position = "top", legend.title = element_blank()) +
        #     guides(colour = guide_legend(nrow = 3, ncol = 2, byrow = TRUE))
        #
        #
        # if (length(mcn) > 1)
        #     scorePlot <- scorePlot + geom_bar(stat = "identity", position = "dodge",
        #                                       aes_string(colour = "merged", fill = "merged"))
        # else
        #     scorePlot <- scorePlot + geom_bar(stat = "identity", aes_string(colour = "type", fill = "type"))
        #
        # scorePlot <- scorePlot + coord_flip()

        # ap <- align_plots(MSPlot, scorePlot, align = "h")
        # plot_grid(ap[[1]], plot_grid(structPlot, ap[[2]], nrow = 2, rel_heights = c(1, 2)), rel_widths = c(2, 1))
        # gridExtra::grid.arrange(MSPlot, structPlot, scorePlot, layout_matrix = matrix(c(1, 2, 1, 3, 1, 3), ncol = 2, byrow = TRUE), widths = c(2, 1))

    }
})

setMethod("plotSpecHash", "compounds", function(obj, index, groupName, MSPeakLists, formulas = NULL,
                                                plotStruct = TRUE, title = NULL, useGGPlot2 = FALSE, xlim = NULL,
                                                ylim = NULL, ...)
{
    compTable <- compoundTable(obj)[[groupName]]
    cRow <- if (is.null(compTable) || nrow(compTable) == 0) NULL else compTable[index, ]
    return(makeHash(cRow, annotatedPeakList(obj, index, groupName, MSPeakLists, formulas),
                    plotStruct, title, useGGPlot2, xlim, ylim, ...))
})

#' @describeIn compounds plots a Venn diagram (using \pkg{\link{VennDiagram}})
#'   outlining unique and shared compound candidates of up to five different
#'   \code{compounds} objects. Comparison is made on \code{InChIKey1}.
#'
#' @inheritParams plotVenn,formulas-method
#'
#' @template plotvenn-ret
#'
#' @export
setMethod("plotVenn", "compounds", function(obj, ..., labels = NULL, vennArgs = NULL)
{
    allCompounds <- c(list(obj), list(...))

    ac <- checkmate::makeAssertCollection()
    checkmate::assertList(allCompounds, types = "compounds", min.len = 2, any.missing = FALSE,
                          unique = TRUE, .var.name = "...", add = ac)
    checkmate::assertCharacter(labels, min.chars = 1, len = length(allCompounds), null.ok = TRUE, add = ac)
    checkmate::assertList(vennArgs, names = "unique", null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    if (is.null(labels))
        labels <- make.unique(sapply(allCompounds, algorithm))
    if (is.null(vennArgs))
        vennArgs <- list()

    allCompoundTabs <- lapply(allCompounds, as.data.table)
    do.call(makeVennPlot, c(list(allCompoundTabs, labels, lengths(allCompounds), function(obj1, obj2)
    {
        if (length(obj1) == 0 || length(obj2) == 0)
            return(data.table())
        fintersect(obj1[, c("group", "InChIKey1")], obj2[, c("group", "InChIKey1")])
    }, nrow), vennArgs))
})

#' @describeIn compounds plots an UpSet diagram (using the
#'   \code{\link[UpSetR]{upset}} function) outlining unique and shared compound
#'   candidates between different \code{compounds} objects. Comparison is made
#'   on \code{InChIKey1}.
#'
#' @inheritParams plotUpSet,formulas-method
#'
#' @references \insertRef{Conway2017}{patRoon} \cr\cr
#'   \insertRef{Lex2014}{patRoon}
#'
#' @export
setMethod("plotUpSet", "compounds", function(obj, ..., labels = NULL, nsets = length(list(...)) + 1,
                                             nintersects = NA, upsetArgs = NULL)
{
    allCompounds <- c(list(obj), list(...))

    ac <- checkmate::makeAssertCollection()
    checkmate::assertList(allCompounds, types = "compounds", min.len = 2, any.missing = FALSE,
                          unique = TRUE, .var.name = "...", add = ac)
    checkmate::assertCharacter(labels, min.chars = 1, len = length(allCompounds), null.ok = TRUE, add = ac)
    checkmate::assertList(upsetArgs, names = "unique", null.ok = TRUE, add = ac)
    checkmate::assertCount(nsets, positive = TRUE)
    checkmate::assertCount(nintersects, positive = TRUE, na.ok = TRUE)
    checkmate::reportAssertions(ac)

    if (is.null(labels))
        labels <- make.unique(sapply(allCompounds, algorithm))

    allCompsTabs <- mapply(allCompounds, labels, SIMPLIFY = FALSE, FUN = function(f, l)
    {
        ret <- as.data.table(f)
        if (length(ret) == 0)
            ret <- data.table(group = character(), InChIKey1 = character())
        ret <- unique(ret[, c("group", "InChIKey1")])[, (l) := 1]
    })

    compTab <- Reduce(function(f1, f2)
    {
        merge(f1, f2, by = c("group", "InChIKey1"), all = TRUE)
    }, allCompsTabs)

    compTab <- compTab[, labels, with = FALSE]
    for (j in seq_along(compTab))
        set(compTab, which(is.na(compTab[[j]])), j, 0)

    if (sum(sapply(compTab, function(x) any(x>0))) < 2)
        stop("Need at least two non-empty objects to plot")

    do.call(UpSetR::upset, c(list(compTab, nsets = nsets, nintersects = nintersects), upsetArgs))
})


setMethod("mergedCompoundNames", "compounds", function(compounds) character(0))
setMethod("mergedCompoundNames", "compoundsConsensus", function(compounds) compounds@mergedCompNames)

#' @templateVar what compounds
#' @template consensus-form_comp
#'
#' @param minMaxNormalization Set to \code{TRUE} to apply min-max normalization
#'   of (merged) scoring columns. \code{FALSE} will apply normalization to the
#'   maximum value. Scorings with negative values will always be min-max
#'   normalized.
#'
#' @templateVar what compounds
#' @template consensus-common-args
#'
#' @return \code{consensus} returns a \code{compounds} object that is produced
#'   by merging multiple specified \code{compounds} objects.
#'
#' @export
setMethod("consensus", "compounds", function(obj, ..., absMinAbundance = NULL,
                                             relMinAbundance = NULL,
                                             uniqueFrom = NULL, uniqueOuter = FALSE,
                                             minMaxNormalization = FALSE,
                                             rankWeights = 1, labels = NULL)
{
    allCompounds <- c(list(obj), list(...))

    ac <- checkmate::makeAssertCollection()
    checkmate::assertList(allCompounds, types = "compounds", min.len = 2, any.missing = FALSE,
                          unique = TRUE, .var.name = "...", add = ac)
    checkmate::assertNumeric(rankWeights, lower = 0, finite = TRUE, add = ac)
    checkmate::assertCharacter(labels, min.chars = 1, len = length(allCompounds), null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    rankWeights <- rep(rankWeights, length.out = length(allCompounds))
    compNames <- if (!is.null(labels)) labels else sapply(allCompounds, algorithm)
    if (anyDuplicated(compNames))
    {
        # duplicate algorithms used, try to form unique names by adding library
        dbs <- lapply(allCompounds, function(cmp) # UNDONE: make this a method
        {
            for (res in compoundTable(cmp))
            {
                if (nrow(res) > 0)
                    return(res$database[[1]])
            }
        })

        compNames <- sapply(seq_along(allCompounds), function(cmpi) paste0(substr(algorithm(allCompounds[[cmpi]]), 1, 3), "-",
                                                                           substr(dbs[[cmpi]], 1, 3)))

        # in case names are still duplicated
        compNames <- make.unique(compNames)
    }

    assertConsCommonArgs(absMinAbundance, relMinAbundance, uniqueFrom, uniqueOuter, compNames)

    relMinAbundance <- max(NULLToZero(absMinAbundance) / length(allCompounds), NULLToZero(relMinAbundance))

    # initialize all compound objects for merge: copy them, rename columns to
    # avoid duplicates and set merged by field of fragInfo.
    allCompTables <- lapply(seq_along(allCompounds), function(cmpi)
    {
        mergedBy <- compNames[[cmpi]]

        return(lapply(compoundTable(allCompounds[[cmpi]]), function(ct)
        {
            ret <- copy(ct)

            for (r in seq_len(nrow(ret)))
            {
                fi <- ret[["fragInfo"]][[r]]
                if (!is.null(fi) && nrow(fi) > 0 && is.null(fi[["mergedBy"]]))
                {
                    fi <- copy(fi)
                    set(fi, j = "mergedBy", value = mergedBy)
                    set(ret, r, "fragInfo", list(list(fi)))
                }
            }

            ret[, mergedBy := compNames[cmpi]]
            ret[, rank := seq_len(.N)]
            ret[, rankscore := (.N - (rank - 1)) / .N * rankWeights[cmpi]]

            setnames(ret, paste0(names(ret), "-", compNames[[cmpi]]))

            return(ret)
        }))
    })

    # columns that should be unique (fragInfo and InChIKey1 are dealt separately)
    uniqueCols <- c("SMILES", "formula", "InChI", "InChIKey2", "InChIKey", "neutralMass")

    leftName <- compNames[[1]]
    mCompList <- allCompTables[[1]]
    for (compIndex in seq(2, length(allCompTables)))
    {
        rightName <- compNames[[compIndex]]

        printf("Merging %s with %s... ", paste0(compNames[seq_len(compIndex-1)], collapse = ","), rightName)

        rightTable <- allCompTables[[compIndex]]

        for (grp in union(names(mCompList), names(rightTable)))
        {
            if (is.null(rightTable[[grp]]))
                next # nothing to merge
            else if (is.null(mCompList[[grp]])) # not yet present
            {
                mCompounds <- rightTable[[grp]]

                # rename columns that should be unique from right to left
                unCols <- c(uniqueCols, "fragInfo", "InChIKey1", "mergedBy")
                unCols <- unCols[sapply(unCols, function(uc) !is.null(mCompounds[[paste0(uc, "-", rightName)]]))]
                setnames(mCompounds, paste0(unCols, "-", rightName), paste0(unCols, "-", leftName))
            }
            else
            {
                mCompounds <- merge(mCompList[[grp]], rightTable[[grp]],
                                    by.x = paste0("InChIKey1-", leftName),
                                    by.y = paste0("InChIKey1-", rightName),
                                    all = TRUE)

                # merge fragment info
                fiColLeft <- paste0("fragInfo-", leftName)
                fiColRight <- paste0("fragInfo-", rightName)

                if (!is.null(mCompounds[[fiColLeft]]))
                {
                    for (r in seq_len(nrow(mCompounds)))
                    {
                        # use copy as a workaround for buggy nested data.tables
                        fiRLeft <- copy(mCompounds[[fiColLeft]][[r]])
                        fiRRight <- copy(mCompounds[[fiColRight]][[r]])
                        hasLeft <- length(fiRLeft) > 0 && nrow(fiRLeft) > 0
                        hasRight <- length(fiRRight) > 0 && nrow(fiRRight) > 0

                        if (hasLeft && hasRight)
                        {
                            # both have fraginfo
                            fiMerged <- mergeFragInfo(fiRLeft, fiRRight, leftName, rightName)
                            set(mCompounds, r, fiColLeft, list(list(fiMerged)))
                        }
                        else if (hasRight) # only right
                            set(mCompounds, r, fiColLeft, list(list(fiRRight)))
                    }

                    mCompounds[, (fiColRight) := NULL]
                }

                # remove duplicate columns that shouldn't
                for (col in uniqueCols)
                {
                    colLeft <- paste0(col, "-", leftName)
                    colRight <- paste0(col, "-", rightName)
                    if (!is.null(mCompounds[[colRight]]))
                    {
                        if (is.null(mCompounds[[colLeft]]))
                            setnames(mCompounds, colRight, colLeft)
                        else
                        {
                            mCompounds[, (colLeft) := ifelse(!is.na(get(colLeft)), get(colLeft), get(colRight))]
                            mCompounds[, (colRight) := NULL]
                        }
                    }
                }

                # collapse mergedBy
                ml <- paste0("mergedBy-", leftName); mr <- paste0("mergedBy-", rightName)
                mCompounds[!is.na(get(ml)), (ml) := ifelse(!is.na(get(mr)), paste(get(ml), get(mr), sep = ","), get(ml))]
                mCompounds[is.na(get(ml)), (ml) := get(mr)]
                mCompounds[, (mr) := NULL]
            }

            mCompList[[grp]] <- mCompounds
        }

        cat("Done!\n")
    }

    printf("Determining coverage and final scores... ")

    # rename & merge score types and ranges
    scoreTypes <- Reduce(union, mapply(allCompounds, compNames, FUN = function(cmp, cn)
    {
        paste0(cmp@scoreTypes, "-", cn)
    }))

    scRanges <- Reduce(modifyList, mapply(allCompounds, compNames, SIMPLIFY = FALSE, FUN = function(cmp, cn)
    {
        lapply(cmp@scoreRanges, function(scrg) setNames(scrg, paste0(names(scrg), "-", cn)))
    }))

    # Determine coverage of compounds between objects and the merged score. The score column can be
    # used for the former as there is guaranteed to be one for each merged object.
    for (grpi in seq_along(mCompList))
    {
        # fix up de-duplicated column names
        deDupCols <- c(uniqueCols, c("fragInfo", "InChIKey1", "mergedBy"))
        leftCols <- paste0(deDupCols, "-", leftName)
        deDupCols <- deDupCols[leftCols %in% names(mCompList[[grpi]])]
        leftCols <- leftCols[leftCols %in% names(mCompList[[grpi]])]
        if (length(leftCols) > 0)
            setnames(mCompList[[grpi]], leftCols, deDupCols)

        mCompList[[grpi]][, coverage := sapply(mergedBy, function(mb) (countCharInStr(mb, ",") + 1) / length(allCompounds))]

        if (relMinAbundance > 0)
            mCompList[[grpi]] <- mCompList[[grpi]][coverage >= relMinAbundance]
        else if (!is.null(uniqueFrom))
        {
            if (!is.character(uniqueFrom))
                uniqueFrom <- compNames[uniqueFrom]

            keep <- function(mergedBy)
            {
                mbs <- unlist(strsplit(mergedBy, ","))
                return(all(mbs %in% uniqueFrom) && (!uniqueOuter || length(mbs) == 1))
            }

            mCompList[[grpi]] <- mCompList[[grpi]][mCompList[[grpi]][, sapply(mergedBy, keep)]]
        }

        rnames <- getAllCompCols("rankscore", names(mCompList[[grpi]]), compNames)
        mCompList[[grpi]][, rankscore := rowSums(.SD, na.rm = TRUE) / length(rnames), .SDcols = rnames]
        setorderv(mCompList[[grpi]], "rankscore", order = -1)
        mCompList[[grpi]][, c(rnames, "rankscore") := NULL]
    }

    cat("Done!\n")

    # prune empty/NULL results
    if (length(mCompList) > 0)
        mCompList <- mCompList[sapply(mCompList, function(r) !is.null(r) && nrow(r) > 0, USE.NAMES = FALSE)]

    return(compoundsConsensus(compounds = mCompList, scoreTypes = scoreTypes, scoreRanges = scRanges,
                              algorithm = paste0(unique(sapply(allCompounds, algorithm)), collapse = ", "),
                              mergedCompNames = compNames))
})

#' @templateVar func generateCompounds
#' @templateVar what generate compounds
#' @templateVar ex1 generateCompoundsMetfrag
#' @templateVar ex2 generateCompoundsSIRIUS
#' @templateVar algos metfrag,sirius
#' @template generic-algo
#'
#' @param ... Any parameters to be passed to the selected compound generation
#'   algorithm.
#'
#' @rdname compound-generation
#' @aliases generateCompounds
#' @export
setMethod("generateCompounds", "featureGroups", function(fGroups, MSPeakLists, algorithm, ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertChoice(algorithm, c("metfrag", "sirius"), add = ac)
    checkmate::reportAssertions(ac)
    
    f <- switch(algorithm,
                metfrag = generateCompoundsMetfrag,
                sirius = generateCompoundsSIRIUS)

    f(fGroups, MSPeakLists, ...)
})
