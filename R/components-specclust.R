#' @include main.R
#' @include components-clust.R
NULL

componentsSpecClust <- setClass("componentsSpecClust", contains = "componentsClust")


#' @export
setMethod("generateComponentsSpecClust", "featureGroups", function(fGroups, MSPeakLists, method = "complete",
                                                                   maxTreeHeight = 1, deepSplit = TRUE,
                                                                   minModuleSize = 1)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertString(method, add = ac)
    assertDynamicTreeCutArgs(maxTreeHeight, deepSplit, minModuleSize, ac)
    checkmate::reportAssertions(ac)
    
    if (length(fGroups) == 0)
        return(componentsSpecClust(distm = NULL, method = method, gInfo = groupInfo(fGroups),
                                   properties = list(), maxTreeHeight = maxTreeHeight, deepSplit = deepSplit,
                                   minModuleSize = minModuleSize, algorithm = "specclust"))

    cat("Calculating distance matrix... ")
    
    MSPeakLists <- MSPeakLists[, intersect(groupNames(MSPeakLists), groupNames(fGroups))]
    allMSMS <- pruneList(sapply(averagedPeakLists(MSPeakLists), "[[", "MSMS", simplify = FALSE))
    
    if (F)
    {
        # UNDONE: parameters for spec similarity
        distm <- 1 - proxy::simil(allMSMS, method = function(x, y) specSimilarityR(x, y, "cosine"))
        class(distm) <- "dist" # has both simul and dist, which confuses S4 validity checks
    }
    else
        distm <- 1 - as.dist(specDistMatrix(allMSMS, "cosine", "none", 0, 1, 0.002))
    cat("Done!\n")
    
    gInfo <- groupInfo(fGroups)[names(allMSMS), ] # make sure to subset!
    
    # UNDONE: properties: spec similarity properties
    return(componentsSpecClust(distm = distm, method = method, gInfo = gInfo, properties = list(),
                               maxTreeHeight = maxTreeHeight, deepSplit = deepSplit,
                               minModuleSize = minModuleSize, algorithm = "specclust"))
})
