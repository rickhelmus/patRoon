#' @include main.R
#' @include feature_groups.R
NULL

#' @rdname featureGroups-class
#' @export
featureGroupsXCMS3 <- setClass("featureGroupsXCMS3", slots = c(xdata = "ANY"), contains = "featureGroups")

setMethod("initialize", "featureGroupsXCMS3",
          function(.Object, ...) callNextMethod(.Object, algorithm = "xcms3", ...))


#' Group features using XCMS (new interface)
#'
#' Uses the new \code{xcms3} interface from the \pkg{xcms} package to find features.
#'
#' @templateVar algo XCMS3
#' @templateVar do group features
#' @templateVar generic groupFeatures
#' @templateVar algoParam xcms3
#' @template algo_generator
#' 
#' @details Grouping of features and alignment of their retention times are performed with the
#'   \code{\link[xcms:groupChromPeaks]{xcms::groupChromPeaks}} and \code{\link[xcms:adjustRtime]{xcms::adjustRtime}}
#'   functions, respectively. Both of these functions support an extensive amount of parameters that modify their
#'   behavior and may therefore require optimization.
#'
#' @template feat-arg
#' @template rtalign-arg
#' @template loadrawdata-arg
#'
#' @param groupParam,retAlignParam parameter object that is directly passed to
#'   \code{\link[xcms:groupChromPeaks]{xcms::groupChromPeaks}} and \code{\link[xcms:adjustRtime]{xcms::adjustRtime}},
#'   respectively.
#' @param preGroupParam grouping parameters applied when features are grouped \emph{prior} to alignment (only with peak
#'   groups alignment).
#'
#' @inheritParams groupFeatures
#'
#' @inherit groupFeatures return
#'
#' @template sets-loadrawdata-RTalign-note
#'
#' @references \addCitations{xcms}{1} \cr\cr \addCitations{xcms}{2} \cr\cr \addCitations{xcms}{3}
#'
#' @templateVar what groupFeaturesXCMS3
#' @templateVar cl features
#' @template main-rd-method
#' @export
setMethod("groupFeaturesXCMS3", "features", function(feat, rtalign = TRUE, loadRawData = TRUE,
                                                     groupParam = xcms::PeakDensityParam(sampleGroups = analysisInfo(feat)$group),
                                                     preGroupParam = groupParam,
                                                     retAlignParam = xcms::ObiwarpParam(), verbose = TRUE)
{
    # UNDONE: aligning gives XCMS errors when multithreading

    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(feat, "features", add = ac)
    aapply(checkmate::assertFlag, . ~ rtalign + loadRawData + verbose, fixed = list(add = ac))
    aapply(assertS4, . ~ groupParam + preGroupParam + retAlignParam, fixed = list(add = ac))
    checkmate::reportAssertions(ac)

    xdata <- getXCMSnExp(feat, verbose = verbose, loadRawData = loadRawData)
    return(doGroupFeaturesXCMS3(xdata, feat, rtalign, loadRawData, groupParam, preGroupParam, retAlignParam, verbose))
})

#' @rdname groupFeaturesXCMS3
#' @export
setMethod("groupFeaturesXCMS3", "featuresSet", function(feat,
                                                        groupParam = xcms::PeakDensityParam(sampleGroups = analysisInfo(feat)$group),
                                                        verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFlag(verbose, add = ac)
    assertS4(groupParam, add = ac)
    checkmate::reportAssertions(ac)
    
    # HACK: force non-set features method to allow grouping of neutralized features
    # UNDONE: or simply export this functionality with a flag?
    xdata <- selectMethod("getXCMSnExp", "features")(feat, verbose = verbose, loadRawData = FALSE)
    
    return(doGroupFeaturesXCMS3(xdata, feat, rtalign = FALSE, loadRawData = FALSE, groupParam, groupParam,
                                xcms::ObiwarpParam(), verbose))
})

doGroupFeaturesXCMS3 <- function(xdata, feat, rtalign, loadRawData, groupParam, preGroupParam, retAlignParam, verbose)
{
    anaInfo <- analysisInfo(feat)
    
    if (length(feat) == 0)
        return(featureGroupsXCMS(features = feat))
    
    hash <- makeHash(feat, rtalign, loadRawData, groupParam, preGroupParam, retAlignParam)
    cachefg <- loadCacheData("featureGroupsXCMS3", hash)
    if (!is.null(cachefg))
        return(cachefg)
    
    if (verbose)
        cat("Grouping features with XCMS...\n===========\n")
    
    if (!loadRawData && rtalign)
    {
        if (verbose)
            cat("Skipping RT alignment: no raw data\n")
    }
    else if (rtalign)
    {
        # group prior to alignment when using the peaks group algorithm
        if (isXCMSClass(retAlignParam, "PeakGroupsParam"))
        {
            if (verbose)
                cat("Performing grouping prior to retention time alignment...\n")
            xdata <- verboseCall(xcms::groupChromPeaks, list(xdata, preGroupParam), verbose)
        }
        
        if (verbose)
            cat("Performing retention time alignment...\n")
        xdata <- verboseCall(xcms::adjustRtime, list(xdata, retAlignParam), verbose)
    }
    
    if (verbose)
        cat("Performing grouping...\n")
    xdata <- verboseCall(xcms::groupChromPeaks, list(xdata, groupParam), verbose)
    
    ret <- importFeatureGroupsXCMS3FromFeat(xdata, anaInfo, feat)
    saveCacheData("featureGroupsXCMS3", ret, hash)
    
    if (verbose)
        cat("\n===========\nDone!\n")
    return(ret)
}

getFeatIndicesFromXCMSnExp <- function(xdata)
{
    plist <- as.data.table(xcms::chromPeaks(xdata))
    plist[, subind := seq_len(.N), by = "sample"]

    xdftidx <- xcms::featureValues(xdata, value = "index")
    ret <- as.data.table(t(xdftidx))

    for (grp in seq_along(ret))
    {
        idx <- plist[ret[, grp, with = FALSE][[1]], subind]
        set(ret, j = grp, value = idx)
    }

    ret[is.na(ret)] <- 0

    return(ret)
}

importFeatureGroupsXCMS3FromFeat <- function(xdata, analysisInfo, feat)
{
    xcginfo <- xcms::featureDefinitions(xdata)
    gNames <- makeFGroupName(seq_len(nrow(xcginfo)), xcginfo[, "rtmed"], xcginfo[, "mzmed"])
    gInfo <- data.frame(rts = xcginfo[, "rtmed"], mzs = xcginfo[, "mzmed"], row.names = gNames, stringsAsFactors = FALSE)

    groups <- data.table(t(xcms::featureValues(xdata, value = "maxo")))
    setnames(groups, gNames)
    groups[is.na(groups)] <- 0

    ret <- featureGroupsXCMS3(xdata = xdata, groups = groups, groupInfo = gInfo, features = feat,
                              ftindex = setnames(getFeatIndicesFromXCMSnExp(xdata), gNames))
    
    # synchronize features: any that were without group have been removed
    ret@xdata <- xcms::filterChromPeaks(ret@xdata, getKeptXCMSPeakInds(feat, ret@features))
    
    return(ret)
}

#' Imports feature groups from XCMS (new interface)
#'
#' Imports grouped features from a \code{\link{XCMSnExp}} object from the \pkg{xcms} package.
#'
#' @template analysisInfo-arg
#' @param xdata An \code{\link{XCMSnExp}} object.
#' 
#' @inherit groupFeaturesXCMS3 references
#' @inherit importFeatureGroups return
#'
#' @seealso \code{\link{groupFeatures}}
#' 
#' @export
importFeatureGroupsXCMS3 <- function(xdata, analysisInfo)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(xdata, "XCMSnExp", add = ac)
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, c("mzXML", "mzML"), add = ac)
    checkmate::reportAssertions(ac)

    if (length(xcms::hasFeatures(xdata)) == 0)
        stop("Provided XCMS data does not contain any grouped features!")

    feat <- importFeaturesXCMS3(xdata, analysisInfo)
    return(importFeatureGroupsXCMS3FromFeat(xdata, analysisInfo, feat))
}

#' @rdname featureGroups-class
#' @export
setMethod("delete", "featureGroupsXCMS3", function(obj, ...)
{
    # UNDONE: update individual features somehow too?
    
    old <- obj
    obj <- callNextMethod()

    if (length(old) > length(obj))
        obj@xdata <- xcms::filterFeatureDefinitions(obj@xdata, names(old) %chin% names(obj))
    
    # simple ana subset
    if (!setequal(analyses(old), analyses(obj)))
        obj@xdata <- xcms::filterFile(obj@xdata, which(analyses(old) %in% analyses(obj)), keepFeatures = TRUE)
    
    if (nrow(xcms::chromPeaks(obj@xdata)) != length(obj@features)) # sync features
        obj@xdata <- xcms::filterChromPeaks(obj@xdata, getKeptXCMSPeakInds(old, obj@features))

    return(obj)
})
