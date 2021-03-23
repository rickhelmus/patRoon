#' @include main.R
#' @include components-clust.R
NULL

componentsSpecClust <- setClass("componentsSpecClust", contains = "componentsClust")

#' @export
setMethod("generateComponentsSpecClust", "featureGroups", function(fGroups, MSPeakLists, method = "complete",
                                                                   specSimParams = getDefSpecSimParams(),
                                                                   maxTreeHeight = 1, deepSplit = TRUE,
                                                                   minModuleSize = 1)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    assertSpecSimParams(specSimParams, add = ac)
    assertDynamicTreeCutArgs(maxTreeHeight, deepSplit, minModuleSize, ac)
    checkmate::reportAssertions(ac)
    
    if (length(fGroups) == 0)
        return(componentsSpecClust(distm = NULL, method = method, gInfo = groupInfo(fGroups),
                                   properties = list(), maxTreeHeight = maxTreeHeight, deepSplit = deepSplit,
                                   minModuleSize = minModuleSize, algorithm = "specclust"))

    gNames <- names(fGroups)
    
    cat("Calculating distance matrix... ")
    sims <- spectrumSimilarity(MSPeakLists, gNames, NULL, MSLevel = 2, specSimParams = specSimParams, NAToZero = TRUE,
                               drop = FALSE)
    
    # figure out fGroups with results: these must have non-zero columns (or rows), since there must be at least a 1.0
    # similarity with itself.
    grpsResults <- gNames[colSums(sims) > 0]
    sims <- sims[grpsResults, grpsResults]
    
    distm <- 1 - as.dist(sims)
    cat("Done!\n")
    
    gInfo <- groupInfo(fGroups)[grpsResults, ]

    return(componentsSpecClust(distm = distm, method = method, gInfo = gInfo,
                               properties = list(specSimParams = specSimParams),
                               maxTreeHeight = maxTreeHeight, deepSplit = deepSplit,
                               minModuleSize = minModuleSize, algorithm = "specclust"))
})
