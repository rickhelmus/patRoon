#' @include features.R
NULL

#' @rdname features-class
#' @export
featuresXCMS <- setClass("featuresXCMS", slots = list(xs = "xcmsSet"), contains = "features")

setMethod("initialize", "featuresXCMS",
          function(.Object, ...) callNextMethod(.Object, algorithm = "xcms", ...))

#' @rdname features-class
#' @export
setMethod("[", c("featuresXCMS", "ANY", "missing", "missing"), function(x, i, j, ..., drop = TRUE)
{
    x <- callNextMethod(x, i, j, ..., drop)
    if (length(analyses(x)) > 0)
        x@xs <- x@xs[, analyses(x)]
    else
        warning("XCMS currently does not allow removing data from all analyses by subsetting.")
    return(x)
})

#' @rdname features-class
#' @export
setMethod("filter", "featuresXCMS", function(obj, ...)
{
    obj <- callNextMethod(obj, ...)

    # check if amount of features (peaks) changed (e.g. due to filtering), if so update
    if (length(obj) != nrow(peaks(obj@xs)))
    {
        cat("Updating xcmsSet...\n")
        # NOTE: use base method to force update as overloaded method simply returns @xs slot
        obj@xs <- selectMethod(getXCMSSet, "features")(obj, TRUE)
    }

    return(obj)
})


#' @details \code{findFeaturesXCMS} uses the \code{\link[xcms]{xcmsSet}}
#'   function from the \pkg{xcms} package to find features.
#'
#' @param method The method setting used by XCMS peak finding, see
#'   \code{\link[xcms:findPeaks-methods]{xcms::findPeaks}}
#'
#' @references \addCitations{xcms}{1} \cr\cr
#'    \addCitations{xcms}{2} \cr\cr
#'    \addCitations{xcms}{3}
#'
#' @rdname feature-finding
#' @export
findFeaturesXCMS <- function(analysisInfo, method = "centWave", ..., verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, c("mzXML", "mzML"), add = ac)
    checkmate::assertString(method, min.chars = 1, add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)

    files <- sapply(seq_len(nrow(analysisInfo)),
                    function(i) getMzMLOrMzXMLAnalysisPath(analysisInfo$analysis[i], analysisInfo$path[i]),
                    USE.NAMES = FALSE)

    hash <- makeHash(analysisInfo, do.call(makeFileHash, as.list(files)), method, list(...))
    cachef <- loadCacheData("featuresXCMS", hash)
    if (!is.null(cachef))
        return(cachef)

    if (verbose)
        printf("Finding features with XCMS for %d analyses ...\n", nrow(analysisInfo))

    if (verbose)
        xs <- xcmsSet(files, analysisInfo$analysis, analysisInfo$group, method = method, ...)
    else
        suppressMessages(xs <- xcmsSet(files, analysisInfo$analysis, analysisInfo$group, method = method, ...))

    ret <- importFeaturesXCMS(xs, analysisInfo)

    saveCacheData("featuresXCMS", ret, hash)

    if (verbose)
    {
        printf("Done!\n")
        printFeatStats(ret@features)
    }

    return(ret)
}

#' @details \code{importFeaturesXCMS} converts features from an existing
#'   \code{\link{xcmsSet}} object (obtained with the \pkg{xcms} package)
#'   to a new \code{\link{features}} object.
#'
#' @param xs An \code{\link{xcmsSet}} object.
#'
#' @rdname feature-finding
#' @export
importFeaturesXCMS <- function(xs, analysisInfo)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(xs, "xcmsSet", add = ac)
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, c("mzXML", "mzML"), add = ac)
    checkmate::reportAssertions(ac)

    feat <- importXCMSPeaks(xcms::peaks(xs), analysisInfo)

    return(featuresXCMS(xs = xs, features = feat, analysisInfo = analysisInfo))
}
