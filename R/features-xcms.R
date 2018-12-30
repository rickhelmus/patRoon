#' @include features.R
NULL

#' @rdname features-class
#' @export
featuresXCMS <- setClass("featuresXCMS", slots = list(xs = "xcmsSet"), contains = "features")

setMethod("initialize", "featuresXCMS",
          function(.Object, ...) callNextMethod(.Object, algorithm = "xcms", ...))

#' @export
setMethod("[", c("featuresXCMS", "ANY", "missing", "missing"), function(x, i, j, ..., drop = TRUE)
{
    x <- callNextMethod(x, i, j, ..., drop)
    x@xs <- x@xs[, analyses(x)]
    return(x)
})

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
#' @note The file format of analyses for \code{findFeaturesXCMS} must be
#'   \code{mzML} or \code{mzXML}.
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
    assertAnalysisInfo(analysisInfo, c("mzXML", "mzML"), add = ac)
    checkmate::assertString(method, min.chars = 1, add = ac)
    checkmate::reportAssertions(ac)

    hash <- makeHash(analysisInfo, list(...)) # UNDONE: better hash?
    cachef <- loadCacheData("featuresXCMS", hash)
    if (!is.null(cachef))
        return(cachef)

    if (verbose)
        cat("Finding features with XCMS...\n===========\n")

    files <- sapply(seq_len(nrow(analysisInfo)),
                    function(i) getMzMLOrMzXMLAnalysisPath(analysisInfo$analysis[i], analysisInfo$path[i]),
                    USE.NAMES = FALSE)
    if (verbose)
        xs <- xcmsSet(files, analysisInfo$analysis, analysisInfo$group, method = method, ...)
    else
        suppressMessages(xs <- xcmsSet(files, analysisInfo$analysis, analysisInfo$group, method = method, ...))

    ret <- importFeaturesXCMS(xs, analysisInfo)

    saveCacheData("featuresXCMS", ret, hash)

    if (verbose)
        cat("\n===========\nDone!\n")

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
    assertAnalysisInfo(analysisInfo, add = ac)
    checkmate::reportAssertions(ac)

    plist <- as.data.table(peaks(xs))
    snames <- sampnames(xs)
    feat <- list()

    feat <- lapply(seq_along(snames), function(sind)
    {
        ret <- plist[sample == sind]
        ret[, ID := seq_len(nrow(ret))]
        setnames(ret, c("rt", "rtmin", "rtmax", "maxo", "into"), c("ret", "retmin", "retmax", "intensity", "area"))
        return(ret[, c("mz", "mzmin", "mzmax", "ret", "retmin", "retmax", "intensity", "area", "ID")])
    })
    names(feat) <- snames

    return(featuresXCMS(xs = xs, features = feat, analysisInfo = analysisInfo))
}
