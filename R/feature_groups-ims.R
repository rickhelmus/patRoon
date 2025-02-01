#' @include features.R
#' @include feature_groups.R
NULL

#' @rdname featureGroups-class
#' @export
featureGroupsIMS <- setClass("featureGroupsIMS", contains = "featureGroups")

setMethod("initialize", "featureGroupsIMS",
          function(.Object, groupAlgo, ...) callNextMethod(.Object, algorithm = paste0("ims-", groupAlgo), ...))

setMethod("groupFeaturesIMS", "features", function(feat, groupAlgo, IMSWindow = 0.02, ..., verbose = FALSE)
{
    # UNDONE: doc verbose bahavior
    # UNDONE: check that all features are orphans?
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertString(groupAlgo, add = ac)
    checkmate::assertNumber(IMSWindow, lower = 0, finite = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (!hasMobilities(feat))
        stop("Features have no mobilities assigned!", call. = FALSE)
    
    anaInfo <- analysisInfo(feat)
    fTable <- featureTable(feat)
    fTableAll <- as.data.table(feat)
    
    # clusters features with similar mobilities
    if (nrow(fTableAll) == 0)
        fTableAll[, cl := integer()]
    else if (nrow(fTableAll) == 1)
        fTableAll[, cl := 1L]
    else
    {
        hc <- fastcluster::hclust(dist(fTableAll$mobility))
        fTableAll[, cl := cutree(hc, h = IMSWindow)]
    }
    fgIMSClusts <- lapply(seq_len(max(fTableAll$cl)), function(clust)
    {
        featSub <- delete(feat, j = function(ft, ana, ...) !ft$ID %chin% fTableAll[analysis == ana & cl == clust]$ID)
        ret <- groupFeatures(featSub, groupAlgo, ..., verbose = verbose)
        ft <- as.data.table(getFeatures(ret))
        ret@groupInfo[, mobility := mean(ft$mobility[ft$group == group]), by = "group"]
        return(ret)
    })
    
    fgInfoAll <- rbindlist(lapply(fgIMSClusts, function(fg) copy(groupInfo(fg))), idcol = "IMSClust")
    fgInfoAll[, group_ims := makeIMSFGroupName(.I, ret, mz, mobility)]

    gTable <- data.table()
    gTable[, (fgInfoAll$group_ims) := Map(fgInfoAll$IMSClust, fgInfoAll$group, f = function(cl, g)
    {
        ints <- numeric(length(fTable))
        fg <- fgIMSClusts[[cl]][, g]
        ints[match(analyses(fg), anaInfo$analysis)] <- fg[[1]]
        return(ints)
    })]
    
    ftind <- data.table()
    ftind[, (fgInfoAll$group_ims) := Map(fgInfoAll$IMSClust, fgInfoAll$group, f = function(cl, g)
    {
        inds <- integer(length(fTable))
        fg <- removeEmptyAnalyses(fgIMSClusts[[cl]][, g])
        anaInds <- match(analyses(fg), anaInfo$analysis)
        inds[anaInds] <- mapply(featureTable(fg), anaInds, FUN = function(ft, ai) match(ft$ID, fTable[[ai]]$ID))
        return(inds)
    })]
    
    gInfo <- subsetDTColumnsIfPresent(copy(fgInfoAll), c("group_ims", "ret", "mz", "mobility"))
    gInfo[, ims_parent_group := NA_character_]
    setnames(gInfo, "group_ims", "group")
    
    return(featureGroupsIMS(groupAlgo = groupAlgo, groups = gTable, groupInfo = gInfo, features = feat,
                            ftindex = ftind))
})
