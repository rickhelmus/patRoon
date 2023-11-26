#' @include features.R
NULL

#' @rdname features-class
#' @export
featuresXCMS <- setClass("featuresXCMS", slots = list(xs = "ANY"), contains = "features")

setMethod("initialize", "featuresXCMS",
          function(.Object, ...) callNextMethod(.Object, algorithm = "xcms", ...))

setMethod("reorderAnalyses", "featuresXCMS", function(obj, anas)
{
    obj <- callNextMethod()
    obj@xs <- obj@xs[, anas]
    return(obj)
})

#' @rdname features-class
#' @export
setReplaceMethod("analysisInfo", "featuresXCMS", function(obj, value)
{
    obj <- callNextMethod()
    xcms::sampclass(obj@xs) <- analysisInfo(obj)$group # sync
    return(obj)
})

#' @rdname features-class
#' @export
setMethod("delete", "featuresXCMS", function(obj, i = NULL, j = NULL, ...)
{
    old <- obj
    obj <- callNextMethod()
    
    # simple ana subset
    if (is.null(j) && !setequal(analyses(old), analyses(obj)))
        obj@xs <- obj@xs[, analyses(old) %in% analyses(obj)]
    else if (!is.null(j)) # sync features
        xcms::peaks(obj@xs) <- xcms::peaks(obj@xs)[getKeptXCMSPeakInds(old, obj), , drop = FALSE]
    
    return(obj)
})

#' Find features using XCMS (old interface)
#'
#' Uses the legacy \code{\link[xcms]{xcmsSet}} function from the \pkg{xcms} package to find features.
#'
#' @templateVar algo XCMS
#' @templateVar do automatically find features
#' @templateVar generic findFeatures
#' @templateVar algoParam xcms
#' @template algo_generator
#'
#' @details This function uses the legacy interface of \pkg{xcms}. It is recommended to use
#'   \code{\link{findFeaturesXCMS3}} instead.
#'
#'   The file format of analyses must be \code{mzML} or \code{mzXML}.
#'
#' @template centroid_note_mandatory
#'
#' @inheritParams findFeatures
#'
#' @param method The method setting used by XCMS peak finding, see \code{\link[xcms:findPeaks-methods]{xcms::findPeaks}}
#' @param \dots Further parameters passed to \code{\link[xcms]{xcmsSet}}.
#'
#' @references \addCitations{xcms}{1} \cr\cr \addCitations{xcms}{2} \cr\cr \addCitations{xcms}{3}
#'
#' @inherit findFeatures return
#'
#' @seealso \code{\link{findFeaturesXCMS3}}
#'
#' @export
findFeaturesXCMS <- function(analysisInfo, method = "centWave", ..., verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, verifyCentroided = TRUE, c("mzXML", "mzML"), add = ac)
    checkmate::assertString(method, min.chars = 1, add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)

    files <- sapply(seq_len(nrow(analysisInfo)),
                    function(i) getMzMLOrMzXMLAnalysisPath(analysisInfo$analysis[i], analysisInfo$path[i]),
                    USE.NAMES = FALSE)

    hash <- makeHash(analysisInfo[, c("analysis", "path", "group"), with = FALSE], do.call(makeFileHash, as.list(files)),
                     method, list(...))
    cachef <- loadCacheData("featuresXCMS", hash)
    if (!is.null(cachef))
        return(cachef)

    if (verbose)
        printf("Finding features with XCMS for %d analyses ...\n", nrow(analysisInfo))

    if (verbose)
        xs <- xcms::xcmsSet(files, analysisInfo$analysis, analysisInfo$group, method = method, ...)
    else
        suppressMessages(xs <- xcms::xcmsSet(files, analysisInfo$analysis, analysisInfo$group, method = method, ...))

    ret <- importFeaturesXCMS(xs, analysisInfo)

    saveCacheData("featuresXCMS", ret, hash)

    if (verbose)
    {
        printf("Done!\n")
        printFeatStats(ret@features)
    }

    return(ret)
}

#' Imports features from XCMS (old interface)
#'
#' Imports feature data generated with the legacy \code{\link[xcms]{xcmsSet}} function from the \pkg{xcms} package.
#'
#' @templateVar algo XCMS
#' @templateVar generic importFeatures
#' @templateVar algoParam xcms
#' @template algo_importer
#'
#' @inheritParams importFeatures
#' 
#' @param xs An \code{\link{xcmsSet}} object.
#' 
#' @inherit findFeaturesXCMS references
#' @inherit importFeatures return
#'
#' @seealso \code{\link{importFeaturesXCMS3}}
#'
#' @export
importFeaturesXCMS <- function(xs, analysisInfo)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(xs, "xcmsSet", add = ac)
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, c("mzXML", "mzML"), verifyCentroided = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    feat <- importXCMSPeaks(xcms::peaks(xs), analysisInfo)

    return(featuresXCMS(xs = xs, features = feat, analysisInfo = analysisInfo))
}
