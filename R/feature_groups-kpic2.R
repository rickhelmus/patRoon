#' @include main.R
#' @include feature_groups.R
NULL

#' @rdname featureGroups-class
#' @export
featureGroupsKPIC2 <- setClass("featureGroupsKPIC2", slots = c(picsSetGrouped = "list"), contains = "featureGroups")

setMethod("initialize", "featureGroupsKPIC2",
          function(.Object, ...) callNextMethod(.Object, algorithm = "kpic2", ...))


#' @rdname feature-grouping
#' @export
groupFeaturesKPIC2 <- function(feat, rtalign = TRUE, exportedData = TRUE,
                               groupArgs = list(tolerance = c(0.005, 12)),
                               alignArgs = list(), verbose = TRUE)
{
    checkPackage("KPIC", "https://github.com/hcji/KPIC2")
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(feat, "features", add = ac)
    aapply(checkmate::assertFlag, . ~ rtalign + exportedData + verbose, fixed = list(add = ac))
    aapply(checkmate::assertList, . ~ groupArgs + alignArgs, any.missing = FALSE, names = "unique",
           fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    if (length(feat) == 0)
        return(featureGroupsKPIC2(analysisInfo = analysisInfo(feat), features = feat))
    
    hash <- makeHash(feat, rtalign, exportedData, groupArgs, alignArgs)
    cachefg <- loadCacheData("featureGroupsKPIC2", hash)
    if (!is.null(cachefg))
        return(cachefg)
    
    if (verbose)
        cat("Grouping features with KPIC2... ")
    
    picsSet <- getPICSet(feat, exportedData = exportedData)
    picsSetGrouped <- do.call(KPIC::PICset.group, c(list(picsSet), groupArgs))
    
    if (!exportedData && rtalign)
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
    # NOTE: the group.info from KPCI2 is not used since retention times and m/zs are rounded...
    
    peakMat <- copy(picsSetGrouped$peakmat)
    
    gInfo <- peakMat[, .(rts = mean(rt), mzs = mean(mz)), by = "group"]
    setorderv(gInfo, "mzs")
    gNames <- makeFGroupName(seq_len(nrow(gInfo)), gInfo$rts, gInfo$mzs)
    
    # NOTE: KPIC2 group ID may not be continuous
    peakMat[, groupName := gNames[match(group, gInfo$group)]]
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

    gInfo <- as.data.frame(gInfo[, -"group"])
    rownames(gInfo) <- gNames
    
    return(featureGroupsKPIC2(picsSetGrouped = picsSetGrouped, groups = groups, groupInfo = gInfo,
                              analysisInfo = analysisInfo, features = feat,
                              ftindex = ftindex))
}

#' @rdname feature-grouping
#' @export
importFeatureGroupsKPIC2 <- function(picsSetGrouped, analysisInfo)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertList(picsSetGrouped, add = ac)
    checkmate::assertNames(names(picsSetGrouped), must.include = c("group.info", "peakmat", "picset"), add = ac)
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, c("mzXML", "mzML"), add = ac)
    checkmate::reportAssertions(ac)
    
    feat <- importfeaturesKPIC2(picsSetGrouped$picset, analysisInfo)
    return(importFeatureGroupsKPIC2FromFeat(picsSetGrouped, analysisInfo, feat))
}

setMethod("removeGroups", "featureGroupsKPIC2", function(fGroups, indices)
{
    fGroups <- callNextMethod(fGroups, indices)
    
    # UNDONE: sync picsSetGroups slot?
    
    return(fGroups)
})
