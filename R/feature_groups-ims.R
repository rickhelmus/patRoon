#' @include features.R
#' @include feature_groups.R
NULL

#' @rdname featureGroups-class
#' @export
featureGroupsIMS <- setClass("featureGroupsIMS", contains = "featureGroups")

setMethod("initialize", "featureGroupsIMS",
          function(.Object, groupAlgo, ...) callNextMethod(.Object, algorithm = paste0("ims-", groupAlgo), ...))

# UNDONE: this was an attempt to use classical grouping algorithms with IMS features. It first clusters all features by
# mobility and then does the grouping and merges the final groups. The problem is that mobility clustering cannot work
# well with reasonable amounts of features. Maybe fix this someday?
setMethod("groupFeaturesIMS", "features", function(feat, grouper, groupAlgo, ..., IMSWindow, verbose)
{
    hash <- makeHash(feat, groupAlgo, list(...), IMSWindow)
    cd <- loadCacheData("groupFeaturesIMS", hash)
    if (!is.null(cd))
        return(cd)
    
    anaInfo <- analysisInfo(feat)
    fTable <- featureTable(feat)
    fTableAll <- as.data.table(feat)
    
    if (any(!is.na(fTableAll$ims_parent_ID)))
        warning("Any links between IMS parents and mobility features will be removed!", call. = FALSE)
        
    # clusters features with similar mobilities
    fTableAll <- clusterFTableMobilities(feat, IMSWindow, byGroup = FALSE)
    
    if (!isFALSE(verbose))
        printf("Grouping features in %d IMS clusters... \n", uniqueN(fTableAll$IMSClust))
    
    fgIMSClusts <- doApply("lapply", doPar = FALSE, prog = isTRUE(verbose), unique(fTableAll$IMSClust), \(clust)
    {
        if (is.na(clust))
        {
            # special case: these are non-IMS features
            featSub <- delete(feat, j = \(ft, ...) !is.na(ft$mobility))
        }
        else
            featSub <- delete(feat, j = \(ft, ana, ...) !ft$ID %chin% fTableAll[analysis == ana & IMSClust == clust]$ID)
        ret <- grouper(featSub, ..., verbose = if (identical(verbose, "full")) verbose else FALSE)
        ft <- as.data.table(getFeatures(ret))
        ret@groupInfo[, mobility := mean(ft$mobility[ft$group == group]), by = "group"]
        if (!is.null(ft[["CCS"]]))
            ret@groupInfo[, CCS := mean(ft$CCS[ft$group == group]), by = "group"]
        doProgress()
        return(ret)
    })
    
    fgInfoAll <- rbindlist(lapply(fgIMSClusts, function(fg) copy(groupInfo(fg))), idcol = "IMSClust")
    setorderv(fgInfoAll, c("ret", "mz", "mobility"), na.last = FALSE)
    fgInfoAll[, group_ims := fifelse(!is.na(mobility), makeIMSFGroupName(.I, ret, mz, mobility), makeFGroupName(.I, ret, mz))]

    gTable <- data.table()
    gTable[, (fgInfoAll$group_ims) := Map(fgInfoAll$IMSClust, fgInfoAll$group, f = \(cl, g) fgIMSClusts[[cl]][[g]])]

    ftind <- data.table()
    ftind[, (fgInfoAll$group_ims) := Map(fgInfoAll$IMSClust, fgInfoAll$group, f = function(cl, g)
    {
        mapply(featureTable(fgIMSClusts[[cl]]), seq_len(nrow(anaInfo)), FUN = function(ft, ai)
        {
            ID <- ft$ID[ft$group == g]
            if (length(ID) == 0) 0L else match(ID, fTable[[ai]]$ID, nomatch = 0L)
        })
    })]
    
    gInfo <- subsetDTColumnsIfPresent(copy(fgInfoAll), c("group_ims", "ret", "mz", "mobility", "CCS"))
    gInfo[, ims_parent_group := NA_character_]
    setnames(gInfo, "group_ims", "group")
    setcolorder(gInfo, "CCS", after = last(names(gInfo)), skip_absent = TRUE) # put CCS at the end of the table
    
    ret <- featureGroupsIMS(groupAlgo = groupAlgo, groups = gTable, groupInfo = gInfo, features = feat,
                            ftindex = ftind)
    saveCacheData("groupFeaturesIMS", ret, hash)
    return(ret)
})
