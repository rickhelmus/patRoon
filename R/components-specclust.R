#' @include main.R
#' @include components-clust.R
NULL

componentsSpecClust <- setClass("componentsSpecClust", contains = "componentsClust")

#' @export
setMethod("generateComponentsSpecClust", "featureGroups", function(fGroups, MSPeakLists, method = "complete",
                                                                   simMethod, shift = "none", removePrecursor = FALSE,
                                                                   mzWeight = 0, intWeight = 1, absMzDev = 0.005,
                                                                   relMinIntensity = 0.05, minSimMSMSPeaks = 0,
                                                                   maxTreeHeight = 1, deepSplit = TRUE,
                                                                   minModuleSize = 1)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertString(method, add = ac)
    checkmate::assertChoice(simMethod, c("cosine", "jaccard"), add = ac)
    checkmate::assertChoice(shift, c("none", "precursor", "both"), add = ac)
    checkmate::assertFlag(removePrecursor, add = ac)
    aapply(checkmate::assertNumber, . ~ mzWeight + intWeight + absMzDev + relMinIntensity,
           finite = TRUE, fixed = list(add = ac))
    checkmate::assertCount(minSimMSMSPeaks, add = ac)
    assertDynamicTreeCutArgs(maxTreeHeight, deepSplit, minModuleSize, ac)
    checkmate::reportAssertions(ac)
    
    if (length(fGroups) == 0)
        return(componentsSpecClust(distm = NULL, method = method, gInfo = groupInfo(fGroups),
                                   properties = list(), maxTreeHeight = maxTreeHeight, deepSplit = deepSplit,
                                   minModuleSize = minModuleSize, algorithm = "specclust"))

    gNames <- names(fGroups)
    
    cat("Calculating distance matrix... ")
    sims <- spectrumSimilarity(MSPeakLists, gNames, NULL, MSLevel = 2, method = simMethod,
                               shift = shift, removePrecursor = removePrecursor,
                               mzWeight = mzWeight, intWeight = intWeight, absMzDev = absMzDev,
                               relMinIntensity = relMinIntensity, minPeaks = minSimMSMSPeaks,
                               NAToZero = TRUE, drop = FALSE)
    
    # figure out fGroups with results: these must have non-zero columns (or rows), since there must be at least a 1.0
    # similarity with itself.
    grpsResults <- gNames[colSums(sims) > 0]
    sims <- sims[grpsResults, grpsResults]
    
    distm <- 1 - as.dist(sims)
    cat("Done!\n")
    
    gInfo <- groupInfo(fGroups)[grpsResults, ]

    return(componentsSpecClust(distm = distm, method = method, gInfo = gInfo,
                               properties = list(simMethod = simMethod, shift = shift,
                                                 removePrecursor = removePrecursor,
                                                 mzWeight = mzWeight, intWeight = intWeight,
                                                 absMzDev = absMzDev),
                               maxTreeHeight = maxTreeHeight, deepSplit = deepSplit,
                               minModuleSize = minModuleSize, algorithm = "specclust"))
})
