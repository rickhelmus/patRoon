#' @include main.R
#' @include workflow-step.R
NULL

printFeatStats <- function(fList)
{
    fCounts <- sapply(fList, nrow)
    fTotCount <- sum(fCounts)
    printf("Feature statistics:\n")
    printf("%s: %d (%.1f%%)\n", names(fList), fCounts, if (fTotCount == 0) 0 else fCounts * 100 / fTotCount)
    printf("Total: %d\n", fTotCount)
}

#' Base features class
#'
#' Holds information for all features present within a set of analysis.
#'
#' This class provides a way to store intensity, retention times, \emph{m/z} and
#' other data for all features in a set of analyses. The class is \code{virtual}
#' and derived objects are created by 'feature finders' such as
#' \code{findFeaturesOpenMS}, \code{findFeaturesXCMS} and
#' \code{findFeaturesBruker}.
#'
#' @param obj,x,object \code{features} object to be accessed
#'
#' @seealso \code{\link{feature-finding}}
#'
#' @slot features List of features per analysis file. Use the
#'   \code{featureTable} method for access.
#' @slot analysisInfo Analysis group information. Use the \code{analysisInfo} method
#'   for access.
#'
#' @templateVar seli analyses
#' @templateVar selOrderi analyses()
#' @templateVar dollarOpName analysis
#' @template sub_op-args
#'
#' @templateVar class features
#' @template class-hierarchy
#'
#' @export
features <- setClass("features",
                     slots = c(features = "list", analysisInfo = "data.frame"),
                     contains = c("VIRTUAL", "workflowStep"))

#' @describeIn features Obtain total number of features.
#' @export
setMethod("length", "features", function(x) if (length(x@features) > 0) sum(sapply(x@features, nrow)) else 0)

#' @describeIn features Shows summary information for this object.
#' @export
setMethod("show", "features", function(object)
{
    callNextMethod(object)
    ftcounts <- if (length(object@features) > 0) sapply(object@features, nrow) else 0
    printf("Total feature count: %d\n", sum(ftcounts))
    printf("Average feature count/analysis: %.0f\n", if (length(object) > 0) sum(ftcounts) / nrow(analysisInfo(object)) else 0)
    printf("Least features: %s\n", names(object)[which.min(ftcounts)])
    printf("Most features: %s\n", names(object)[which.max(ftcounts)])
    showAnaInfo(analysisInfo(object))
})

#' @describeIn features Get table with feature information
#'
#' @return \code{featureTable}: A \code{list} containing a
#'   \code{\link{data.table}} for each analysis with feature data
#'
#' @export
setMethod("featureTable", "features", function(obj) obj@features)

#' @describeIn features Get analysis information
#' @return \code{analysisInfo}: A \code{data.frame} containing a column with
#'   analysis name (\code{analysis}), its path (\code{path}), and other columns
#'   such as replicate group name (\code{group}) and blank reference
#'   (\code{ref}).
#' @export
setMethod("analysisInfo", "features", function(obj) obj@analysisInfo)

#' @templateVar class features
#' @templateVar what analyses
#' @template strmethod
#' @export
setMethod("analyses", "features", function(obj) analysisInfo(obj)$analysis)

#' @templateVar class features
#' @templateVar what replicate groups
#' @template strmethod
#' @export
setMethod("replicateGroups", "features", function(obj) unique(analysisInfo(obj)$group))

#' @describeIn features Returns all feature data in a table.
#' @export
setMethod("as.data.table", "features", function(x) rbindlist(featureTable(x), idcol = "analysis", fill = TRUE))

#' @describeIn features Performs common rule based filtering of features. Note
#'   that this (and much more) functionality is also provided by the
#'   \code{filter} method defined for \code{\link{featureGroups}}. However,
#'   filtering a \code{features} object may be useful to avoid grouping large
#'   amounts of features.
#' @templateVar feat TRUE
#' @template feat-filter-args
#' @export
setMethod("filter", "features", function(obj, absMinIntensity = NULL, relMinIntensity = NULL,
                                         retentionRange = NULL, mzRange = NULL, mzDefectRange = NULL,
                                         chromWidthRange = NULL, negate = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertNumber, . ~ absMinIntensity + relMinIntensity, lower = 0, finite = TRUE,
           null.ok = TRUE, fixed = list(add = ac))
    aapply(assertRange, . ~ retentionRange + mzRange + mzDefectRange + chromWidthRange, null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertFlag(negate, add = ac)
    checkmate::reportAssertions(ac)

    if (length(obj) == 0)
        return(obj)

    oldn <- length(obj)

    hash <- makeHash(obj, absMinIntensity, relMinIntensity, retentionRange, mzRange, mzDefectRange, chromWidthRange, negate)
    cache <- loadCacheData("filterFeatures", hash)
    if (!is.null(cache))
        obj <- cache
    else
    {
        anaInfo <- analysisInfo(obj)

        absIntPred <- if (!negate) function(x) x >= absMinIntensity else function(x) x < absMinIntensity
        relIntPred <- if (!negate) function(x, m) (x / m) >= relMinIntensity else function(x, m) (x / m) < relMinIntensity
        rangePred <- function(x, range) numGTE(x, range[1]) & numLTE(x, range[2])

        if (negate)
            rangePred <- Negate(rangePred)

        for (ana in analyses(obj))
        {
            if (!is.null(absMinIntensity))
                obj@features[[ana]] <- obj@features[[ana]][absIntPred(intensity)]

            if (!is.null(relMinIntensity))
            {
                maxInt <- max(obj@features[[ana]]$intensity)
                obj@features[[ana]] <- obj@features[[ana]][relIntPred(intensity, maxInt)]
            }
            
            if (!is.null(retentionRange))
                obj@features[[ana]] <- obj@features[[ana]][rangePred(ret, retentionRange)]

            if (!is.null(mzRange))
                obj@features[[ana]] <- obj@features[[ana]][rangePred(mz, mzRange)]

            if (!is.null(mzDefectRange))
                obj@features[[ana]] <- obj@features[[ana]][rangePred(mz - floor(mz), mzDefectRange)]
            
            if (!is.null(chromWidthRange))
                obj@features[[ana]] <- obj@features[[ana]][rangePred(retmax - retmin, chromWidthRange)]
        }

        saveCacheData("filterFeatures", obj, hash)
    }

    newn <- length(obj)
    printf("Done! Filtered %d (%.2f%%) features. Remaining: %d\n", oldn - newn, if (oldn == 0) 0 else (1-(newn/oldn))*100, newn)

    return(obj)
})

#' @describeIn features Subset on analyses.
#' @param \dots Ignored.
#' @export
setMethod("[", c("features", "ANY", "missing", "missing"), function(x, i, ...)
{
    if (!missing(i))
    {
        i <- assertSubsetArgAndToChr(i, analyses(x))
        x@features <- x@features[i]
        x@analysisInfo <- x@analysisInfo[x@analysisInfo$analysis %in% i, ]
    }

    return(x)
})

#' @describeIn features Extract a feature table for an analysis.
#' @export
setMethod("[[", c("features", "ANY", "missing"), function(x, i)
{
    assertExtractArg(i)
    return(x@features[[i]])
})

#' @describeIn features Extract a feature table for an analysis.
#' @export
setMethod("$", "features", function(x, name)
{
    eval(substitute(x@features$NAME_ARG, list(NAME_ARG = name)))
})

#' @rdname target-screening
#' @export
setMethod("screenTargets", "features", function(obj, targets, rtWindow, mzWindow)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertDataFrame(targets, any.missing = FALSE, min.rows = 1, add = add)
    assertHasNames(targets, c("name", "mz"), add = ac)
    aapply(checkmate::assertNumber, . ~ rtWindow + mzWindow, lower = 0, finite = TRUE, fixed = list(add = ac))
    checkmate::reportAssertions(ac)

    fTable <- featureTable(obj)
    anaInfo <- analysisInfo(obj)

    retlist <- lapply(seq_len(nrow(targets)), function(ti)
    {
        hasRT <- !is.null(targets$rt) && !is.na(targets$rt[ti])

        rbindlist(lapply(names(fTable), function(ana)
        {
            if (hasRT)
                fts <- fTable[[ana]][numLTE(abs(ret - targets$rt[ti]), rtWindow) & numLTE(abs(mz - targets$mz[ti]), mzWindow), ]
            else
                fts <- fTable[[ana]][numLTE(abs(mz - targets$mz[ti]), mzWindow), ]

            if (nrow(fts) == 0) # no results? --> add NA result
                return(data.table(name = targets$name[ti], rt = if (hasRT) targets$rt[ti] else NA,
                                  mz = targets$mz[ti], analysis = ana,
                                  feature = NA, d_rt = NA, d_mz = NA, intensity = NA, area = NA))

            return(rbindlist(lapply(seq_len(nrow(fts)), function(i)
            {
                data.table(name = targets$name[ti], rt = if (hasRT) targets$rt[ti] else NA, mz = targets$mz[ti], analysis = ana,
                           feature = fts[["ID"]][i], d_rt = if (hasRT) fts[["ret"]][i] - targets$rt[ti] else NA,
                           d_mz = fts[["mz"]][i] - targets$mz[ti], intensity = fts[["intensity"]][i],
                           area = if (is.null(fts[["area"]][i])) 0 else fts[["area"]][i])
            })))

        }))
    })

    return(rbindlist(retlist, fill = TRUE))
})

#' @templateVar func findFeatures
#' @templateVar what find features
#' @templateVar ex1 findFeaturesOpenMS
#' @templateVar ex2 findFeaturesBruker
#' @templateVar algos bruker,openms,xcms,envipick
#' @template generic-algo
#'
#' @rdname feature-finding
#' @aliases findFeatures
#' @export
findFeatures <- function(analysisInfo, algorithm, ..., verbose = TRUE)
{
    assertAnalysisInfo(analysisInfo)

    f <- switch(algorithm,
                bruker = findFeaturesBruker,
                openms = findFeaturesOpenMS,
                xcms = findFeaturesXCMS,
                envipick = findFeaturesEnviPick,
                stop("Invalid algorithm! Should be: bruker, openms, xcms or envipick"))

    f(analysisInfo, ..., verbose = verbose)
}
