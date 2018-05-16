#' @include main.R
#' @include compounds.R
NULL

# already done in cluster.R
# setOldClass("hclust")

compoundsCluster <- setClass("compoundsCluster",
                             slots = c(clusters = "list", molecules = "list", cutClusters = "list"),
                             prototype = list(clusters = list(), molecules = list(), cutClusters = list()))


setMethod("makeHCluster", "compounds", function(obj, method, fpType = "extended")
{
    aapply(checkmate::assertString, . ~ method + fpType, min.chars = 1)
    
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
    
    return(compoundsCluster(clusters = clust, molecules = mols))
})
