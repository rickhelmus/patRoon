# SPDX-FileCopyrightText: 2016 - 2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
#' @include feature_groups.R
NULL

#' @rdname featureGroups-class
#' @export
featureGroupsKPIC2 <- setClass("featureGroupsKPIC2", slots = c(picsSetGrouped = "list"), contains = "featureGroups")

setMethod("initialize", "featureGroupsKPIC2",
          function(.Object, ...) callNextMethod(.Object, algorithm = "kpic2", ...))

#' Group features using KPIC2
#'
#' Uses the the \href{https://github.com/hcji/KPIC2}{KPIC2} \R package for grouping of features.
#'
#' @templateVar algo KPIC2
#' @templateVar do group features
#' @templateVar generic groupFeatures
#' @templateVar algoParam kpic2
#' @template algo_generator
#'
#' @details Grouping of features and alignment of their retention times are performed with the
#'   \code{\link[KPIC:PICset.group]{KPIC::PICset.group}} and \code{\link[KPIC:PICset.align]{KPIC::PICset.align}}
#'   functions, respectively.
#'
#' @template feat-arg
#' @template rtalign-arg
#' @template loadrawdata-arg
#'
#' @param groupArgs,alignArgs Named \code{character} vector that may contain extra parameters to be used by
#'   \code{\link[KPIC:PICset.group]{KPIC::PICset.group}} and \code{\link[KPIC:PICset.align]{KPIC::PICset.align}},
#'   respectively.
#'
#' @inheritParams groupFeatures
#'
#' @inherit groupFeatures return
#'
#' @template sets-loadrawdata-RTalign-note
#'
#' @references \insertRef{Ji2017}{patRoon}
#'
#' @templateVar what groupFeaturesKPIC2
#' @templateVar cl features
#' @template main-rd-method
#' @export
setMethod("groupFeaturesKPIC2", "features", function(feat, rtalign = TRUE, loadRawData = TRUE,
                                                     groupArgs = list(tolerance = c(0.005, 12)),
                                                     alignArgs = list(), IMSWindow = defaultLim("mobility", "medium"),
                                                     verbose = TRUE)
{
    checkPackage("KPIC", "rickhelmus/KPIC2")
    
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertFlag, . ~ rtalign + loadRawData, fixed = list(add = ac))
    aapply(checkmate::assertList, . ~ groupArgs + alignArgs, any.missing = FALSE, names = "unique",
           fixed = list(add = ac))
    checkmate::assertNumber(IMSWindow, lower = 0, finite = TRUE, add = ac)
    assertGroupFeatVerbose(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    doG <- \(feat, ...) doGroupFeaturesKPIC2(feat, getPICSet(feat, loadRawData = loadRawData, IMS = "both"), ...)
    return(doGroupFeatures(feat, doG, "kpic2", rtalign = rtalign, loadRawData = loadRawData, groupArgs = groupArgs,
                           alignArgs = alignArgs, IMSWindow = IMSWindow, verbose = verbose))
})

#' @rdname groupFeaturesKPIC2
#' @export
setMethod("groupFeaturesKPIC2", "featuresSet", function(feat, groupArgs = list(tolerance = c(0.005, 12)),
                                                        IMSWindow = defaultLim("mobility", "medium"), verbose = TRUE)
{
    checkPackage("KPIC", "rickhelmus/KPIC2")
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertList(groupArgs, any.missing = FALSE, names = "unique", add = ac)
    checkmate::assertNumber(IMSWindow, lower = 0, finite = TRUE, add = ac)
    assertGroupFeatVerbose(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    doG <- function(feat, ...)
    {
        # HACK: force non-set features method to allow grouping of neutralized features
        # UNDONE: or simply export this functionality with a flag?
        picsSet <- selectMethod("getPICSet", "features")(feat, loadRawData = FALSE, IMS = "both")
        doGroupFeaturesKPIC2(feat, picsSet, ...)
    }
    
    return(doGroupFeatures(feat, doG, "kpic2", rtalign = FALSE, loadRawData = FALSE, groupArgs = groupArgs,
                           alignArgs = list(), IMSWindow = IMSWindow, verbose = verbose))
})

doGroupFeaturesKPIC2 <- function(feat, picsSet, rtalign, loadRawData, groupArgs, alignArgs, verbose)
{
    if (length(feat) == 0)
        return(featureGroupsKPIC2(features = feat))
    
    verbose <- !isFALSE(verbose)
    
    hash <- makeHash(feat, rtalign, loadRawData, groupArgs, alignArgs)
    cachefg <- loadCacheData("featureGroupsKPIC2", hash)
    if (!is.null(cachefg))
        return(cachefg)
    
    if (verbose)
        cat("Grouping features with KPIC2... ")
    
    picsSetGrouped <- do.call(KPIC::PICset.group, c(list(picsSet), groupArgs))
    
    if (!loadRawData && rtalign)
    {
        if (verbose)
            cat("Skipping RT alignment: no raw data\n")
    }
    else if (rtalign)
    {
        picsSetGrouped <- do.call(KPIC::PICset.align, c(list(picsSetGrouped), alignArgs))
        # UNDONE: need to group again?
    }
    
    ret <- importFeatureGroupsKPIC2FromFeat(picsSetGrouped, analysisInfo(feat), feat)
    saveCacheData("featureGroupsKPIC2", ret, hash)
    
    if (verbose)
        cat("Done!\n")
    
    return(ret)
}

importFeatureGroupsKPIC2FromFeat <- function(picsSetGrouped, analysisInfo, feat)
{
    # NOTE: the group.info from KPIC2 is not used since retention times and m/zs are rounded...
    
    if (is.null(picsSetGrouped$peakmat)) # no results (eg filtered out due to frac argument)
    {
        groups <- ftindex <- data.table()
        gInfo <- data.table(group = character(), ret = numeric(), mz = numeric())
    }
    else
    {
        peakMat <- copy(picsSetGrouped$peakmat)
        
        gInfo <- peakMat[, .(ret = mean(rt), mz = mean(mz)), by = "group"]
        setorderv(gInfo, "mz")
        gNames <- makeFGroupName(seq_len(nrow(gInfo)), gInfo$ret, gInfo$mz)
        setnames(gInfo, "group", "group_kp2")
        gInfo[, group := gNames]
        
        # NOTE: KPIC2 group ID may not be continuous
        peakMat[, groupName := gNames[match(group, gInfo$group_kp2)]]
        peakMatSplit <- split(peakMat, by = "groupName")[gNames] # sync order after split
        groups <- data.table()
        groups[, (gNames) := lapply(peakMatSplit, function(grpTab)
        {
            ints <- numeric(nrow(analysisInfo))
            ints[grpTab$sample] <- grpTab$maxo
            return(ints)
        })]
        
        ftindex <- data.table()
        ftindex[, (gNames) := lapply(peakMatSplit, function(grpTab)
        {
            inds <- integer(nrow(analysisInfo))
            inds[grpTab$sample] <- grpTab$index
            return(inds)
        })]
        
        gInfo[, group_kp2 := NULL]
        setcolorder(gInfo, c("group", "ret", "mz"))
    }
    
    return(featureGroupsKPIC2(picsSetGrouped = picsSetGrouped, groups = groups, groupInfo = gInfo, features = feat,
                              ftindex = ftindex))
}

#' Imports feature groups from KPIC2
#'
#' Imports grouped features from an \pkg{KPIC} object.
#'
#' @template analysisInfo-arg
#' @param picsSetGrouped A grouped \code{PIC set} object (\emph{e.g.} as returned by
#'   \code{\link[KPIC:PICset.group]{KPIC::PICset.group}}).
#' 
#' @inherit groupFeaturesKPIC2 references
#' @inherit importFeatureGroups return
#'
#' @seealso \code{\link{groupFeatures}}
#' 
#' @export
importFeatureGroupsKPIC2 <- function(picsSetGrouped, analysisInfo)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertList(picsSetGrouped, add = ac)
    checkmate::assertNames(names(picsSetGrouped), must.include = c("group.info", "peakmat", "picset"), add = ac)
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, fileTypes = "centroid", add = ac)
    checkmate::reportAssertions(ac)
    
    feat <- importFeaturesKPIC2(picsSetGrouped$picset, analysisInfo)
    return(importFeatureGroupsKPIC2FromFeat(picsSetGrouped, analysisInfo, feat))
}

#' @rdname featureGroups-class
#' @export
setMethod("delete", "featureGroupsKPIC2", function(obj, ...)
{
    obj <- callNextMethod()
    
    # UNDONE: sync picsSetGroups slot?
    
    return(obj)
})
