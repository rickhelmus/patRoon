#' @include features.R
NULL

updateXData <- function(obj)
{
    # UNDONE: see if this could be done more efficient/selective...
    
    cat("Updating XCMS object...\n")
    # NOTE: use base method to force update as overloaded method simply returns @xdata slot
    obj@xdata <- selectMethod(getXCMSnExp, "features")(obj, TRUE)
    return(obj)
}

#' @rdname features-class
#' @export
featuresXCMS3 <- setClass("featuresXCMS3", slots = list(xdata = "ANY"), contains = "features")

setMethod("initialize", "featuresXCMS3",
          function(.Object, ...) callNextMethod(.Object, algorithm = "xcms3", ...))

#' @rdname features-class
#' @export
setReplaceMethod("featureTable", "featuresXCMS3", function(obj, value)
{
    ret <- callNextMethod()
    if (!all(mapply(featureTable(obj), featureTable(ret), FUN = XCMSFeatTblEqual)))
        ret <- updateXData(ret)
    return(ret)
})

#' @rdname features-class
#' @export
setMethod("[", c("featuresXCMS3", "ANY", "missing", "missing"), function(x, i, j, ..., drop = TRUE)
{
    x <- callNextMethod(x, i, j, ..., drop = drop)
    x@xdata <- xcms::filterFile(x@xdata, analyses(x))
    return(x)
})

#' @rdname features-class
#' @export
setReplaceMethod("[", c("featuresXCMS3", "ANY", "missing"), function(x, i, j, value)
{
    ret <- callNextMethod()
    if (!all(mapply(featureTable(x), featureTable(ret), FUN = XCMSFeatTblEqual)))
        ret <- updateXData(ret)
    return(ret)
})

#' @rdname features-class
#' @export
setReplaceMethod("[[", c("featuresXCMS3", "ANY", "missing"), function(x, i, j, value)
{
    ret <- callNextMethod()
    if (!XCMSFeatTblEqual(x[[i]], value))
        ret <- updateXData(ret)
    return(ret)
})

#' @rdname features-class
#' @export
setReplaceMethod("$", "featuresXCMS3", function(x, name, value)
{
    ret <- callNextMethod()
    if (!XCMSFeatTblEqual(`$`(ret, name), value))
        ret <- updateXData(ret)
    return(ret)
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
findFeaturesXCMS3 <- function(analysisInfo, param = xcms::CentWaveParam(), verbose = TRUE)
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
        xdata <- xcms::findChromPeaks(rawData, param = param)
    else
        suppressMessages(xdata <- xcms::findChromPeaks(rawData, param = param))

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
