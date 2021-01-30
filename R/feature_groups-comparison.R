#' @include main.R
#' @include features.R
#' @include feature_groups.R
NULL

#' Feature groups comparison class
#'
#' This class is used for comparing different \code{\link{featureGroups}}
#' objects.
#'
#' Objects from this class are returned by \code{\link{comparison}}.
#'
#' @slot fGroupsList A \code{list} of \code{\link{featureGroups}} object that
#'   were compared
#' @slot comparedFGroups A \emph{pseudo} \code{featureGroups} object containing
#'   grouped feature groups.
#'
#' @param x A \code{featureGroupsComparison} object.
#'
#' @templateVar seli labels
#' @templateVar selOrderi names()
#' @templateVar dollarOpName label
#' @template sub_op-args

#' @export
featureGroupsComparison <- setClass("featureGroupsComparison",
                                    slots = c(fGroupsList = "list", comparedFGroups = "featureGroups"))

#' @describeIn featureGroupsComparison Obtain the labels that were given to each compared feature group.
#' @export
setMethod("names", "featureGroupsComparison", function(x) names(x@fGroupsList))

#' @describeIn featureGroupsComparison Number of feature groups objects that were compared.
#' @export
setMethod("length", "featureGroupsComparison", function(x) length(x@fGroupsList))

#' @describeIn featureGroupsComparison Subset on labels that were assigned to compared feature groups.
#' @param \dots Ignored.
#' @export
setMethod("[", c("featureGroupsComparison", "ANY", "missing", "missing"), function(x, i, ...)
{
    if (!missing(i))
    {
        i <- assertSubsetArgAndToChr(i, names(x))
        x@fGroupsList <- x@fGroupsList[i]
        x@comparedFGroups <- x@comparedFGroups[i]
    }

    return(x)
})

#' @describeIn featureGroupsComparison Extract a \code{\link{featureGroups}} object by its label.
#' @export
setMethod("[[", c("featureGroupsComparison", "ANY", "missing"), function(x, i, j)
{
    assertExtractArg(i)
    return(x@fGroupsList[[i]])
})

#' @describeIn featureGroupsComparison Extract a compound table for a feature group.
#' @export
setMethod("$", "featureGroupsComparison", function(x, name)
{
    eval(substitute(x@fGroupsList$NAME_ARG, list(NAME_ARG = name)))
})

#' Comparing feature groups
#'
#' Functionality to compare feature groups and make a consensus.
#'
#' Feature groups objects originating from differing feature finding and/or
#' grouping algorithms (or their parameters) may be compared to assess their
#' output and generate a consensus.
#'
#' @param \dots For \code{comparison}: \code{featureGroups} objects that should
#'   be compared. If the arguments are named (\emph{e.g.} \code{myGroups =
#'   fGroups}) then these are used for labelling, otherwise objects are
#'   automatically labelled by their \code{\link{algorithm}}.
#'
#'   For \code{plot}, \code{plotVenn}, \code{plotChord}: further options passed
#'   to \code{plot}, \pkg{\link{VennDiagram}} plotting functions (\emph{e.g.}
#'   \code{\link{draw.pairwise.venn}}) and \code{\link{chordDiagram}}
#'   respectively.
#'
#'   For \code{plotUpSet}: any further arguments passed to the \code{plotUpSet}
#'   method defined for \code{\link{featureGroups}}.
#' @param which A character vector specifying one or more labels of compared
#'   feature groups. For \code{plotVenn}: if \code{NULL} then all compared
#'   groups are used.

#' @name featureGroups-compare
NULL

#' @rdname featureGroups-compare
featuresFromFeatGroups <- setClass("featuresFromFeatGroups", contains = "features")
setMethod("initialize", "featuresFromFeatGroups",
          function(.Object, ...) callNextMethod(.Object, algorithm = "collapsed_feature_groups", ...))

#' @rdname featureGroups-compare
#' @export
featuresConsensus <- setClass("featuresConsensus", contains = "features")

#' @rdname featureGroups-compare
#' @export
featureGroupsConsensus <- setClass("featureGroupsConsensus", contains = "featureGroups")


# pseudo features: feature groups are converted to features by averaging their
# intensities. Given feature group objects are then considered as an 'analysis'.
convertFeatureGroupsToFeatures <- function(fGroupsList)
{
    fGNames <- names(fGroupsList)

    flist <- list()
    for (fgi in seq_along(fGroupsList))
    {
        if (length(fGroupsList[[fgi]]) == 0)
        {
            flist[[fGNames[fgi]]] <- data.table(ID = character(), ret = numeric(), mz = numeric(), intensity = numeric(),
                                                area = numeric(), retmin = numeric(), retmax = numeric(), mzmin = numeric(),
                                                mzmax = numeric())
            next
        }

        gt <- groupTable(fGroupsList[[fgi]])
        gi <- groupInfo(fGroupsList[[fgi]])

        # use group info as basis
        ft <- as.data.table(gi)
        setnames(ft, c("rts", "mzs"), c("ret", "mz"))
        ft[, ID := colnames(gt)]

        if (nrow(ft) == 0)
            ft[, intensity := 0]
        else
        {
            avgit <- transpose(gt[, lapply(.SD, mean)])
            setnames(avgit, 1, "intensity")
            ft <- cbind(ft, avgit)
        }

        # dummy ranges
        ft[, retmin := ret - 3]
        ft[, retmax := ret + 3]
        ft[, mzmin := mz - 0.0025]
        ft[, mzmax := mz + 0.0025]

        # dummy area
        ft[, area := intensity * 2.5]

        flist[[fGNames[fgi]]] <- ft
    }

    anaInfo <- data.frame(analysis = fGNames, group = fGNames, blank = "", path = ".", stringsAsFactors = FALSE)

    return(featuresFromFeatGroups(features = flist, analysisInfo = anaInfo))
}

#' @details The \code{comparison} method generates a
#'   \code{\link{featureGroupsComparison}} object from given feature groups
#'   objects, which in turn may be used for (visually) comparing presence of
#'   feature groups and generating a consensus. Internally, this function will
#'   collapse each feature groups object to \emph{pseudo} features objects by
#'   averaging their retention times, \emph{m/z} values and intensities, where
#'   each original feature groups object becomes an 'analysis'. All
#'   \emph{pseudo} features are then grouped using
#'   \link[=feature-grouping]{regular feature grouping algorithms} so that a
#'   comparison can be made.
#'
#' @param groupAlgo The \code{\link[=feature-grouping]{feature grouping
#'   algorithm}} that should be used for grouping \emph{pseudo} features (see
#'   details). Valid values are: \code{"xcms"}, \code{xcms3} or \code{"openms"}.
#' @param groupArgs A \code{list} containing further parameters for
#'   \code{\link[=feature-grouping]{feature grouping}}.
#' @param x,obj The \code{featureGroupsComparison} object.
#'
#' @return \code{comparison} returns a \code{\link{featureGroupsComparison}}
#'   object.
#'
#' @rdname featureGroups-compare
#' @aliases comparison
#' @export
setMethod("comparison", "featureGroups", function(..., groupAlgo, groupArgs = list(rtalign = FALSE))
{
    fGroupsList <- list(...)

    ac <- checkmate::makeAssertCollection()
    checkmate::assertChoice(groupAlgo, c("xcms", "xcms3", "openms"), add = ac)
    checkmate::assertList(groupArgs, any.missing = FALSE, names = "unique", add = ac)
    checkmate::assertList(fGroupsList, types = "featureGroups", min.len = 2, any.missing = FALSE,
                          unique = TRUE, .var.name = "...", add = ac)
    checkmate::reportAssertions(ac)

    n <- getArgNames(..., def = sapply(fGroupsList, algorithm))
    names(fGroupsList) <- make.unique(n)

    # convert feature groups to features
    featsFromGroups <- convertFeatureGroupsToFeatures(fGroupsList)

    if (groupAlgo == "xcms" || groupAlgo == "xcms3")
        groupArgs <- c(list(exportedData = FALSE), groupArgs)
    compGroups <- do.call(groupFeatures, c(list(featsFromGroups, groupAlgo), groupArgs))

    return(featureGroupsComparison(fGroupsList = fGroupsList, comparedFGroups = compGroups))
})

#' @details \code{plot} generates an \emph{m/z} \emph{vs} retention time plot.
#' @param retMin If \code{TRUE} retention times are plotted as minutes (seconds otherwise).
#' @rdname featureGroups-compare
#' @export
setMethod("plot", "featureGroupsComparison", function(x, retMin = TRUE, ...) plot(x@comparedFGroups, retMin, ...))

#' @details \code{plotVenn} plots a Venn diagram outlining unique and shared
#'   feature groups between up to five compared feature groups.
#' @template plotvenn-ret
#' @rdname featureGroups-compare
#' @export
setMethod("plotVenn", "featureGroupsComparison", function(obj, which = NULL, ...) plotVenn(obj@comparedFGroups, which, ...))

#' @details \code{plotUpSet} plots an UpSet diagram outlining unique and shared
#'   feature groups.
#' @rdname featureGroups-compare
#' @export
setMethod("plotUpSet", "featureGroupsComparison", function(obj, which = NULL, ...) plotUpSet(obj@comparedFGroups, which, ...))

#' @details \code{plotChord} plots a chord diagram to visualize the distribution
#'   of feature groups.
#'
#' @template plotChord-args
#'
#' @rdname featureGroups-compare
#' @export
setMethod("plotChord", "featureGroupsComparison",
          function(obj, addSelfLinks, addRetMzPlots, ...) plotChord(obj@comparedFGroups, addSelfLinks, addRetMzPlots, average = TRUE, ...))

#' @details \code{consensus} combines all compared feature groups and averages
#'   their retention, \emph{m/z} and intensity data.
#'
#' @templateVar what feature groups
#' @template consensus-common-args
#'
#' @return \code{consensus} returns a \code{\link{featureGroups}} object with a
#'   consensus from the compared feature groups.
#'
#' @rdname featureGroups-compare
#' @export
setMethod("consensus", "featureGroupsComparison", function(obj, absMinAbundance = NULL,
                                                           relMinAbundance = NULL,
                                                           uniqueFrom = NULL, uniqueOuter = FALSE)
{
    # available info:
    # - grouped feature groups --> these are the new consensus feature groups
    # - original feature groups --> from pseudo features
    # - original features --> from original feature groups

    ac <- checkmate::makeAssertCollection()
    assertConsCommonArgs(absMinAbundance, relMinAbundance, uniqueFrom, uniqueOuter, names(obj), add = ac)
    checkmate::reportAssertions(ac)

    allAnaInfos <- lapply(obj@fGroupsList, analysisInfo)

    # UNDONE: is this a limitation?
    if (!all(sapply(allAnaInfos[-1], identical, allAnaInfos[[1]]))) # from https://stackoverflow.com/a/30850654
        stop("This function only works with feature groups with equal analyses")

    anaInfo <- allAnaInfos[[1]]

    # synchronize analyses
    # allAnalyses <- lapply(obj@fGroupsList, function(fg) analysisInfo(fg)$analysis)
    # allAnalyses <- Reduce(intersect, allAnalyses)
    # fGroupsList <- sapply(obj@fGroupsList, function(fg) fg[allAnalyses], simplify = FALSE)
    fGroupsList <- obj@fGroupsList
    comparedFGroups <- obj@comparedFGroups
    
    if (!is.null(absMinAbundance) || !is.null(relMinAbundance))
        comparedFGroups <- minAnalysesFilter(comparedFGroups, absMinAbundance, relMinAbundance, verbose = FALSE)

    if (!is.null(uniqueFrom))
    {
        if (!is.character(uniqueFrom))
            uniqueFrom <- names(obj)[uniqueFrom]
        comparedFGroups <- unique(comparedFGroups, which = uniqueFrom,
                                  outer = uniqueOuter)
        fGroupsList <- fGroupsList[replicateGroups(comparedFGroups)]
    }

    compFeats <- featureTable(comparedFGroups)
    compFeatInds <- groupFeatIndex(comparedFGroups)
    candidates <- names(fGroupsList)
    
    # Add consensus feature groups to original features
    adjFeatures <- sapply(candidates, function(ca)
    {
        adjF <- featureTable(fGroupsList[[ca]])
        adjF <- sapply(names(adjF),
                       function(ana) { ret <- copy(adjF[[ana]]); ret[, consFGroup := NA_character_] }, simplify = FALSE)

        ftind <- groupFeatIndex(fGroupsList[[ca]])
        cai <- match(ca, candidates)

        for (compfgi in seq_along(compFeatInds))
        {
            cfind <- compFeatInds[[compfgi]][cai]
            if (cfind == 0)
                next

            oldFGName <- as.character(compFeats[[ca]][["ID"]][cfind])

            ftinds <- ftind[[oldFGName]]
            for (anai in seq_along(ftinds))
            {
                if (ftinds[anai] != 0)
                    set(adjF[[anaInfo$analysis[anai]]], as.integer(ftinds[anai]), "consFGroup", colnames(compFeatInds)[compfgi])
            }
        }

        # clearout unassigned features
        adjF <- sapply(names(adjF), function(ana) adjF[[ana]][!is.na(consFGroup)], simplify = FALSE)

        return(adjF)
    }, simplify = FALSE)

    # Generate consensus features by merging and averaging all feature data one by one
    cat("Generating consensus features...")
    consFeatures <- sapply(anaInfo$analysis, function(ana)
    {
        # collect all features from all candidates of this analysis
        fts <- lapply(adjFeatures, function(af) af[[ana]])

        fts <- Reduce(function(f1, f2)
        {
            ret <- merge(f1, f2, by = "consFGroup", all = TRUE)

            # make new IDs
            ret[, c("ID.x", "ID.y") := NULL]
            ret[, ID := seq_len(nrow(ret))]

            colsToAvg <- c("ret", "mz", "area", "retmin", "retmax", "mzmin", "mzmax", "intensity")

            for (col in colsToAvg)
                set(ret, j = col, value = rowMeans(ret[, c(paste0(col, ".x"), paste0(col, ".y")), with = FALSE], na.rm = TRUE))

            return(ret[, c("ID", "consFGroup", colsToAvg), with = FALSE]) # only keep relevant columns
        }, fts)

        return(fts)
    }, simplify = FALSE)
    cat("Done!\n")

    allAlgos <- paste0(sapply(fGroupsList, algorithm), collapse = ",")

    retFeatures <- featuresConsensus(features = consFeatures, analysisInfo = anaInfo,
                                     algorithm = allAlgos)

    if (nrow(compFeatInds) == 0) # all input were empty feature groups
        return(featureGroupsConsensus(groups = data.table(), groupInfo = data.frame(),
                                      analysisInfo = anaInfo, features = retFeatures,
                                      ftindex = data.table(),
                                      algorithm = allAlgos))

    # initialize new feature group tables

    consFeatInds <- data.table(matrix(0, nrow = nrow(anaInfo), ncol = ncol(compFeatInds)))
    consFGNames <- colnames(compFeatInds)
    setnames(consFeatInds, consFGNames)
    consGroups <- copy(consFeatInds)

    # Generate consensus feature group tables from consensus features from the consFGroup assignment
    cat("Generating consensus feature groups...\n")
    prog <- openProgBar(0, nrow(anaInfo))
    for (anai in seq_len(nrow(anaInfo)))
    {
        cfts <- consFeatures[[anaInfo$analysis[anai]]]
        for (ftrow in seq_len(nrow(cfts)))
        {
            ft <- cfts[ftrow]
            set(consFeatInds, anai, ft$consFGroup, ftrow)
            set(consGroups, anai, ft$consFGroup, ft$intensity)
        }
        setTxtProgressBar(prog, anai)
    }

    setTxtProgressBar(prog, nrow(anaInfo))
    close(prog)

    return(featureGroupsConsensus(groups = consGroups, analysisInfo = anaInfo,
                                  groupInfo = groupInfo(comparedFGroups), features = retFeatures,
                                  ftindex = consFeatInds, algorithm = allAlgos))
})
