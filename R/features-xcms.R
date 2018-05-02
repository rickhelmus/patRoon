#' @include features.R
NULL

#' @rdname feature-finding
#' @export
featuresXCMS <- setClass("featuresXCMS", slots = list(xs = "xcmsSet"), contains = "features")

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
findFeaturesXCMS <- function(analysisInfo, method = "centWave", ...)
{
    ac <- checkmate::makeAssertCollection()
    assertAnalysisInfo(analysisInfo, add = ac)
    checkmate::assertString(method, min.chars = 1, add = ac)
    checkmate::reportAssertions(ac)
    
    hash <- makeHash(analysisInfo, list(...)) # UNDONE: better hash?
    cachef <- loadCacheData("featuresXCMS", hash)
    if (!is.null(cachef))
        return(cachef)

    cat("Finding features with XCMS...\n===========\n")

    xs <- xcmsSet(sapply(seq_len(nrow(analysisInfo)),
                         function(i) getMzMLOrMzXMLAnalysisPath(analysisInfo$analysis[i], analysisInfo$path[i]),
                         USE.NAMES = FALSE),
                  analysisInfo$analysis, analysisInfo$group, method = method, ...)

    ret <- importFeaturesXCMS(xs, analysisInfo)

    saveCacheData("featuresXCMS", ret, hash)

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
