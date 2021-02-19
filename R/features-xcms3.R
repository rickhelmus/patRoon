#' @include features.R
NULL

#' @rdname features-class
#' @export
featuresXCMS3 <- setClass("featuresXCMS3", slots = list(xdata = "ANY"), contains = "features")

setMethod("initialize", "featuresXCMS3",
          function(.Object, ...) callNextMethod(.Object, algorithm = "xcms3", ...))

#' @rdname features-class
#' @export
setMethod("delete", "featuresXCMS3", function(obj, i = NULL, j = NULL, ...)
{
    old <- obj
    obj <- callNextMethod()
    
    # simple ana subset
    if (is.null(j) && !setequal(analyses(old), analyses(obj)))
        obj@xdata <- xcms::filterFile(obj@xdata, which(analyses(old) %in% analyses(obj)))
    else if (!is.null(j)) # sync features
    {
        # UNDONE: ask for exported method...
        obj@xdata@msFeatureData <- xcms:::.filterChromPeaks(obj@xdata, getKeptXCMSPeakInds(old, obj))
    }
    return(obj)
})

#' @details \code{findFeaturesXCMS3} uses the new \code{xcms3} interface from
#'   the \pkg{xcms} package to find features.
#'
#' @param param The method parameters used by XCMS peak finding, see
#'   \code{\link[xcms:findChromPeaks]{xcms::findChromPeaks}}
#'
#' @references \addCitations{xcms}{1} \cr\cr \addCitations{xcms}{2} \cr\cr
#'   \addCitations{xcms}{3}
#'
#' @rdname feature-finding
#' @export
findFeaturesXCMS3 <- function(analysisInfo, param = xcms::CentWaveParam(), ..., verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, c("mzXML", "mzML"), add = ac)
    assertS4(param, add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)

    files <- sapply(seq_len(nrow(analysisInfo)),
                    function(i) getMzMLOrMzXMLAnalysisPath(analysisInfo$analysis[i], analysisInfo$path[i]),
                    USE.NAMES = FALSE)

    hash <- makeHash(analysisInfo, do.call(makeFileHash, as.list(files)), param)
    cachef <- loadCacheData("featuresXCMS3", hash)
    if (!is.null(cachef))
        return(cachef)

    if (verbose)
        printf("Finding features with XCMS for %d analyses ...\n", nrow(analysisInfo))

    if (verbose)
        printf("Loading raw data...\n")
    rawData <- readMSDataForXCMS3(analysisInfo)

    if (verbose)
        xdata <- xcms::findChromPeaks(rawData, param = param, ...)
    else
        suppressMessages(xdata <- xcms::findChromPeaks(rawData, param = param, ...))

    ret <- importFeaturesXCMS3(xdata, analysisInfo)

    saveCacheData("featuresXCMS3", ret, hash)

    if (verbose)
    {
        printf("Done!\n")
        printFeatStats(ret@features)
    }

    return(ret)
}

#' @details \code{importFeaturesXCMS3} converts features from an existing
#'   \code{\link{XCMSnExp}} object (obtained with the \pkg{xcms} package)
#'   to a new \code{\link{features}} object.
#'
#' @param xdata An \code{\link{XCMSnExp}} object.
#'
#' @rdname feature-finding
#' @export
importFeaturesXCMS3 <- function(xdata, analysisInfo)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(xdata, "XCMSnExp", add = ac)
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, c("mzXML", "mzML"), add = ac)
    checkmate::reportAssertions(ac)

    feat <- importXCMSPeaks(xcms::chromPeaks(xdata), analysisInfo)

    return(featuresXCMS3(xdata = xdata, features = feat, analysisInfo = analysisInfo))
}
