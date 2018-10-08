#' @include main.R
NULL

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
#' @export
features <- setClass("features",
                     slots = c(features = "list", analysisInfo = "data.frame"),
                     prototype = list(features = list(), analysisInfo = data.frame()),
                     contains = "VIRTUAL")

#' @describeIn features Obtain total number of features.
#' @export
setMethod("length", "features", function(x) if (length(x@features) > 0) sum(sapply(x@features, nrow)) else 0)

#' @describeIn features Shows summary information for this object.
#' @export
setMethod("show", "features", function(object)
{
    ftcounts <- if (length(object@features) > 0) sapply(object@features, nrow) else 0
    printf("A features object ('%s')\n", class(object))
    printf("Total feature count: %d\n", sum(ftcounts))
    printf("Average feature count/analysis: %.0f\n", sum(ftcounts) / nrow(analysisInfo(object)))
    printf("Least features: %s\n", names(object)[which.min(ftcounts)])
    printf("Most features: %s\n", names(object)[which.max(ftcounts)])
    showAnaInfo(analysisInfo(object))
    showObjectSize(object)
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

#' @describeIn features Performs common rule based filtering of features.
#' @param intensityThreshold Minimum intensity of a feature. Set to \code{NULL}
#'   to ignore.
#' @param retentionRange,mzRange,chromWidthRange Range of retention time (in
#'   seconds), \emph{m/z} or chromatographic peak width (in seconds),
#'   respectively. Features outside this range will be removed. Should be a
#'   numeric vector with length of two containing the min/max values. If the max
#'   value is set to a value below 0 then no maximum is assumed. Set to
#'   \code{NULL} to skip this step.
#' @param negate If set to \code{TRUE} then filtering operations are performed
#'   in opposite manner.
#' @export
setMethod("filter", "features", function(obj, intensityThreshold = NULL, retentionRange = NULL,
                                         mzRange = NULL, chromWidthRange = NULL, negate = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertNumber(intensityThreshold, lower = 0, finite = TRUE, null.ok = TRUE, add = ac)
    aapply(assertRange, . ~ retentionRange + mzRange + chromWidthRange, null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertFlag(negate, add = ac)
    checkmate::reportAssertions(ac)
    
    if (length(obj) == 0)
        return(obj)
    
    hash <- makeHash(obj, intensityThreshold, retentionRange, mzRange, chromWidthRange, negate)
    cache <- loadCacheData("filterFeatures", hash)
    if (!is.null(cache))
        return(cache)
    
    anaInfo <- analysisInfo(obj)
    
    intPred <- if (!negate) function(x) x >= intensityThreshold else function(x) x < intensityThreshold
    rangePred <- function(x, range)
    {
        if (range[2] < 0)
            numGTE(x, range[1])
        else
            numGTE(x, range[1]) & numLTE(x, range[2])
    }
    
    if (negate)
        rangePred <- Negate(rangePred)
    
    oldn <- length(obj)
    
    for (ana in analyses(obj))
    {
        if (!is.null(intensityThreshold))
            obj@features[[ana]] <- obj@features[[ana]][intPred(intensity)]
        
        if (!is.null(retentionRange))
            obj@features[[ana]] <- obj@features[[ana]][rangePred(ret, retentionRange)]
        
        if (!is.null(mzRange))
            obj@features[[ana]] <- obj@features[[ana]][rangePred(mz, mzRange)]
        
        if (!is.null(chromWidthRange))
            obj@features[[ana]] <- obj@features[[ana]][rangePred(retmax - retmin, chromWidthRange)]
    }
    
    saveCacheData("filterFeatures", obj, hash)
    
    newn <- length(obj)
    printf("Done! Filtered %d (%.2f%%) features. Remaining: %d\n", oldn - newn, if (oldn == 0) 0 else (1-(newn/oldn))*100, newn)
    
    return(obj)
})

#' @describeIn features Subset on analyses.
#' @export
setMethod("[", c("features", "ANY", "missing", "missing"), function(x, i, ...)
{
    if (!missing(i))
    {
        assertSubsetArg(i)
    
        if (!is.character(i))
            i <- analyses(x)[i]
        
        i <- i[i %in% analyses(x)]
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
findFeatures <- function(analysisInfo, algorithm, ...)
{
    assertAnalysisInfo(analysisInfo)
    
    f <- switch(algorithm,
                bruker = findFeaturesBruker,
                openms = findFeaturesOpenMS,
                xcms = findFeaturesXCMS,
                envipick = findFeaturesEnviPick,
                stop("Invalid algorithm! Should be: bruker, openms, xcms or envipick"))

    f(analysisInfo, ...)
}
