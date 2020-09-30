#' @include main.R
#' @include components-clust.R
NULL

componentsSpecClust <- setClass("componentsSpecClust", contains = "componentsClust")

#' @export
setMethod("generateComponentsSpecClust", "featureGroups", function(fGroups, MSPeakLists, method = "complete",
                                                                   simMethod, shift = "none", removePrecursor = FALSE,
                                                                   mzWeight = 0, intWeight = 1, absMzDev = 0.005,
                                                                   relMinIntensity = 0.1, maxTreeHeight = 1, deepSplit = TRUE,
                                                                   minModuleSize = 1)
{
    # UNDONE: document that relative intensity filter is applied after removing precursors
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertString(method, add = ac)
    checkmate::assertChoice(simMethod, c("cosine", "jaccard"), add = ac)
    checkmate::assertChoice(shift, c("none", "precursor", "both"), add = ac)
    checkmate::assertFlag(removePrecursor, add = ac)
    aapply(checkmate::assertNumber, . ~ mzWeight + intWeight + absMzDev + relMinIntensity,
           finite = TRUE, fixed = list(add = ac))
    assertDynamicTreeCutArgs(maxTreeHeight, deepSplit, minModuleSize, ac)
    checkmate::reportAssertions(ac)
    
    if (length(fGroups) == 0)
        return(componentsSpecClust(distm = NULL, method = method, gInfo = groupInfo(fGroups),
                                   properties = list(), maxTreeHeight = maxTreeHeight, deepSplit = deepSplit,
                                   minModuleSize = minModuleSize, algorithm = "specclust"))

    MSPeakLists <- MSPeakLists[, intersect(groupNames(MSPeakLists), groupNames(fGroups))]
    allMSMS <- pruneList(sapply(averagedPeakLists(MSPeakLists), "[[", "MSMS", simplify = FALSE))
    
    if (removePrecursor || relMinIntensity > 0)
    {
        allMSMS <- pruneList(lapply(allMSMS, function(pl)
        {
            if (removePrecursor)
                pl <- pl[precursor == FALSE]
            if (relMinIntensity > 0 && nrow(pl) > 0)
            {
                thr <- relMinIntensity * max(pl$intensity)
                pl <- pl[intensity >= thr]
            }
            return(pl)
        }), checkZeroRows = TRUE)
    }
    
    gInfo <- groupInfo(fGroups)[names(allMSMS), ] # make sure to subset!
    
    precMZs <- sapply(names(allMSMS), function(g) gInfo[g, "mzs"])
    
    cat("Calculating distance matrix... ")
    
    if (F)
    {
        # UNDONE: parameters for spec similarity
        distm <- 1 - proxy::simil(allMSMS, method = function(x, y) specSimilarityR(x, y, simMethod, shift, FALSE, mzWeight,
                                                                                   intWeight, absMzDev, 0))
        class(distm) <- "dist" # has both simul and dist, which confuses S4 validity checks
    }
    else
        distm <- 1 - as.dist(specDistMatrix(allMSMS, simMethod, shift, precMZs, mzWeight, intWeight, absMzDev))
    cat("Done!\n")
    
    return(componentsSpecClust(distm = distm, method = method, gInfo = gInfo,
                               properties = list(simMethod = simMethod, shift = shift,
                                                 removePrecursor = removePrecursor,
                                                 mzWeight = mzWeight, intWeight = intWeight,
                                                 absMzDev = absMzDev),
                               maxTreeHeight = maxTreeHeight, deepSplit = deepSplit,
                               minModuleSize = minModuleSize, algorithm = "specclust"))
})

componentsSpecClustSet <- setClass("componentsSpecClustSet", contains = "componentsClustSet")

#' @export
setMethod("generateComponentsSpecClust", "featureGroupsSet", function(fGroups, ...)
{
    cset <- generateComponentsSet(fGroups, generateComponentsSpecClust, setIonization = FALSE, ...)
    browser()
    return(componentsSpecClustSet(adducts = adducts(cset), setObjects = setObjects(cset),
                                  components = componentTable(cset), componentInfo = componentInfo(cset),
                                  algorithm = "specclust-set"))
})
