# SPDX-FileCopyrightText: 2016 - 2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

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
#' @references \addCitations{xcms}
#'
#' @templateVar what groupFeaturesXCMS3
#' @templateVar cl features
#' @template main-rd-method
#' @export
setMethod("groupFeaturesXCMS3", "features", function(feat, rtalign = TRUE, loadRawData = TRUE,
                                                     groupParam = xcms::PeakDensityParam(sampleGroups = analysisInfo(feat)$replicate),
                                                     preGroupParam = groupParam,
                                                     retAlignParam = xcms::ObiwarpParam(),
                                                     IMSWindow = defaultLim("mobility", "medium"), verbose = TRUE)
{
    # UNDONE: aligning gives XCMS errors when multithreading

    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(feat, "features", add = ac)
    aapply(checkmate::assertFlag, . ~ rtalign + loadRawData, fixed = list(add = ac))
    aapply(assertS4, . ~ groupParam + preGroupParam + retAlignParam, fixed = list(add = ac))
    checkmate::assertNumber(IMSWindow, lower = 0, finite = TRUE, add = ac)
    assertGroupFeatVerbose(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    doG <- function(feat, ...)
    {
        xdata <- getXCMSnExp(feat, verbose = verbose, loadRawData = loadRawData, IMS = "both")
        doGroupFeaturesXCMS3(feat, xdata, ...)
    }
    
    return(doGroupFeatures(feat, doG, "xcms3", rtalign = rtalign, loadRawData = loadRawData, groupParam = groupParam,
                           preGroupParam = preGroupParam, retAlignParam = retAlignParam, IMSWindow = IMSWindow,
                           verbose = verbose))
})

#' @rdname groupFeaturesXCMS3
#' @export
setMethod("groupFeaturesXCMS3", "featuresSet", function(feat,
                                                        groupParam = xcms::PeakDensityParam(sampleGroups = analysisInfo(feat)$replicate),
                                                        IMSWindow = defaultLim("mobility", "medium"), verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    assertS4(groupParam, add = ac)
    checkmate::assertNumber(IMSWindow, lower = 0, finite = TRUE, add = ac)
    assertGroupFeatVerbose(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    doG <- function(feat, ...)
    {
        # HACK: force non-set features method to allow grouping of neutralized features
        # UNDONE: or simply export this functionality with a flag?
        xdata <- selectMethod("getXCMSnExp", "features")(feat, verbose = verbose, loadRawData = FALSE, IMS = "both")
        doGroupFeaturesXCMS3(feat, xdata, ...)
    }
    
    return(doGroupFeatures(feat, doG, "xcms3", rtalign = FALSE, loadRawData = FALSE, groupParam = groupParam,
                           preGroupParam = groupParam, retAlignParam = xcms::ObiwarpParam(), IMSWindow = IMSWindow,
                           verbose = verbose))
})

doGroupFeaturesXCMS3 <- function(feat, xdata, rtalign, loadRawData, groupParam, preGroupParam, retAlignParam, verbose)
{
    anaInfo <- analysisInfo(feat)
    
    if (length(feat) == 0)
        return(featureGroupsXCMS3(features = feat))
    
    verbose <- !isFALSE(verbose)
    
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

    if (nrow(xcginfo) == 0)
    {
        # no results (e.g. due to minFraction>0)
        groups <- ftindex <- data.table()
        gInfo <- data.table(group = character(), ret = numeric(), mzs = numeric())
    }
    else
    {
        gNames <- makeFGroupName(seq_len(nrow(xcginfo)), xcginfo[, "rtmed"], xcginfo[, "mzmed"])
        gInfo <- data.table(group = gNames, ret = xcginfo[, "rtmed"], mz = xcginfo[, "mzmed"])
        
        groups <- data.table(t(xcms::featureValues(xdata, value = "maxo")))
        setnames(groups, gNames)
        groups[is.na(groups)] <- 0
        
        ftindex <- setnames(getFeatIndicesFromXCMSnExp(xdata), gNames)
    }
    
    ret <- featureGroupsXCMS3(xdata = xdata, groups = groups, groupInfo = gInfo, features = feat, ftindex = ftindex)
    
    # synchronize features: any that were without group have been removed
    ret@xdata <- xcms::filterChromPeaks(ret@xdata, getKeptXCMSPeakInds(feat, ret@features))
    
    return(ret)
}

#' Imports feature groups from XCMS (new interface)
#'
#' Imports grouped features from a \code{\link{XCMSnExp}} object from the \pkg{xcms} package.
#'
#' @template analysisInfo-arg
#' @param input An \code{\link{XCMSnExp}} object.
#' 
#' @inherit groupFeaturesXCMS3 references
#' @inherit importFeatureGroups return
#'
#' @seealso \code{\link{groupFeatures}}
#' 
#' @export
importFeatureGroupsXCMS3 <- function(input, analysisInfo)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(input, "XCMSnExp", add = ac)
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, fileTypes = "centroid", add = ac)
    checkmate::reportAssertions(ac)

    if (length(xcms::hasFeatures(input)) == 0)
        stop("Provided XCMS data does not contain any grouped features!")

    feat <- importFeaturesXCMS3(input, analysisInfo)
    return(importFeatureGroupsXCMS3FromFeat(input, analysisInfo, feat))
}

#' @rdname featureGroups-class
#' @export
setMethod("delete", "featureGroupsXCMS3", function(obj, ...)
{
    # UNDONE: update individual features somehow too?
    
    old <- obj
    obj <- callNextMethod()

    if (!hasMobilities(obj))
    {
        if (length(old) > length(obj))
            obj@xdata <- xcms::filterFeatureDefinitions(obj@xdata, names(old) %chin% names(obj))
        # simple ana subset
        if (!setequal(analyses(old), analyses(obj)))
            obj@xdata <- xcms::filterFile(obj@xdata, which(analyses(old) %in% analyses(obj)), keepFeatures = TRUE)
        if (nrow(xcms::chromPeaks(obj@xdata)) != length(obj@features)) # sync features
            obj@xdata <- xcms::filterChromPeaks(obj@xdata, getKeptXCMSPeakInds(old, obj@features))
    }
    
    return(obj)
})
