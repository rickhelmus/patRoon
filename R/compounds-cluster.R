#' @include main.R
#' @include compounds.R
NULL

# This way of compound clustering largely based on metfRag's chemclust.R and the
# package vignette of rcdk

compoundsCluster <- setClass("compoundsCluster",
                             slots = c(clusters = "list", molecules = "list", cutClusters = "list"),
                             prototype = list(clusters = list(), molecules = list(), cutClusters = list()))


setMethod("makeHCluster", "compounds", function(obj, method, fpType = "extended",
                                                maxTreeHeight = 1, deepSplit = TRUE,
                                                minModuleSize = 1)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertString, . ~ method + fpType, min.chars = 1, fixed = list(add = ac))
    checkmate::assertNumber(maxTreeHeight, 0, finite = TRUE, add = ac)
    checkmate::assertFlag(deepSplit, add = ac)
    checkmate::assertCount(minModuleSize, positive = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    compTable <- compoundTable(obj)
    mols <- sapply(compTable, function(ct) rcdk::parse.smiles(ct$SMILES),
                   simplify = FALSE)
    
    clust <- lapply(seq_along(mols), function(i)
    {
        for (j in seq_along(mols[[i]]))
        {
            rcdk::do.typing(mols[[i]][[j]])
            rcdk::do.aromaticity(mols[[i]][[j]])
        }
        
        fps <- lapply(mols[[i]], rcdk::get.fingerprint, type = fpType)
        dist <- as.dist(1 - fingerprint::fp.sim.matrix(fps))
        return(hclust(dist, method))
    })
    
    cutClusters <- lapply(clust, function(d)
    {
        if (minModuleSize == 1)
        {
            # workaround adapted from RAMClustR (ramclustR.R)
            ret <- dynamicTreeCut::cutreeDynamicTree(dendro = d, maxTreeHeight = maxTreeHeight,
                                                     deepSplit = deepSplit, minModuleSize = 2)
            single <- which(ret == 0) # all unassigned have length = 1
            ret[single] <- max(ret) + seq_len(length(single))
        }
        else
            ret <- dynamicTreeCut::cutreeDynamicTree(dendro = d, maxTreeHeight = maxTreeHeight,
                                                     deepSplit = deepSplit, minModuleSize = minModuleSize)
        return(ret)
    })
    
    names(clust) <- names(compTable)
    names(cutClusters) <- names(compTable)
    
    return(compoundsCluster(clusters = clust, molecules = mols, cutClusters = cutClusters))
})

setMethod("plot", "compoundsCluster", function(x, groupName, ...)
{
    assertChoiceSilent(groupName, names(x@clusters))
    plot(x@clusters[[groupName]], ...)
})
