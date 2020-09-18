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
    
    # UNDONE: handle empty objects here or in componentsClust
    # if (length(fGroups) == 0)
    #     return(componentsSpecClust(components = list(), componentInfo = data.table(), clusterm = matrix(),
    #                               distm = structure(list(), class = "dissimilarity"),
    #                               clust = structure(list(), class = "hclust"), cutClusters = numeric(),
    #                               gInfo = data.frame(), properties = list()))
    
    cat("Calculating distance matrix... ")
    
    MSPeakLists <- MSPeakLists[, intersect(groupNames(MSPeakLists), groupNames(fGroups))]
    allMSMS <- pruneList(sapply(averagedPeakLists(MSPeakLists), "[[", "MSMS", simplify = FALSE))
    
    # UNDONE: parameters for spec similarity
    distm <- 1 - proxy::simil(allMSMS, method = function(x, y) specSimilarity(x, y, "cosine"))
    class(distm) <- "dist" # has both simul and dist, which confuses S4 validity checks
    cat("Done!\n")
    
    gInfo <- groupInfo(fGroups)[names(allMSMS), ] # make sure to subset!
    
    # UNDONE: properties: spec similarity properties
    return(componentsSpecClust(distm = distm, method = method, gInfo = gInfo, properties = list(),
                               maxTreeHeight = maxTreeHeight, deepSplit = deepSplit,
                               minModuleSize = minModuleSize, algorithm = "specclust"))
})
