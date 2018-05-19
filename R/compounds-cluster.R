#' @include main.R
#' @include compounds.R
NULL

# This way of compound clustering largely based on metfRag's chemclust.R and the
# package vignette of rcdk

compoundsCluster <- setClass("compoundsCluster",
                             slots = c(clusters = "list", molecules = "list", cutClusters = "list",
                                       properties = "list"),
                             prototype = list(clusters = list(), molecules = list(), cutClusters = list(),
                                              properties = list()))

setMethod("length", "compoundsCluster", function(x) sum(lengths(x)))

setMethod("lengths", "compoundsCluster", function(x, use.names = TRUE) sapply(x@cutClusters,
                                                                              function(cc) length(unique(cc)),
                                                                              USE.NAMES = use.names))

#' @describeIn compounds-clust Show summary information for this object.
#' @export
setMethod("show", "compoundsCluster", function(object)
{
    printf("A compounds cluster object (%s)\n", class(object))
    
    printf("Number of feature groups with compound clusters in this object: %d\n", length(object@clusters))
    
    ls <- lengths(object)
    printf("Number of clusters: %d (total), %.1f (mean), %d - %d (min - max)\n",
           sum(ls), mean(ls), min(ls), max(ls))
    
    printf("Clustering properties:\n")
    printf(" - hclust method: %s\n", object@properties$method)
    printf(" - fingerprint type: %s\n", object@properties$fpType)
    printf(" - fingerprint similarity method: %s\n", object@properties$fpSimMethod)
    
    showObjectSize(object)
})

setMethod("cutCluster", "compoundsCluster", function(obj, k = NULL, h = NULL, groupName)
{
    if (is.null(k) && is.null(h))
        stop("Either k or h should be specified")
    
    ac <- checkmate::makeAssertCollection()
    assertChoiceSilent(groupName, names(obj@clusters), add = ac)
    checkmate::assertCount(k, positive = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertNumber(h, lower = 0, finite = TRUE, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    obj@cutClusters[[groupName]] <- cutree(obj@clusters[[groupName]], k, h)
    return(obj)
})

setMethod("plot", "compoundsCluster", function(x, groupName, pal = "Paired", ...)
{
    assertChoiceSilent(groupName, names(x@clusters))
    
    dend <- as.dendrogram(x@clusters[[groupName]])
    ct <- x@cutClusters[[groupName]]
    ct <- ct[order.dendrogram(dend)] # re-order for dendrogram
    nclust <- length(unique(ct[ct != 0])) # without unassigned
    cols <- getBrewerPal(nclust, pal)
    dend <- dendextend::color_branches(dend, clusters = ct, col = cols)
    lcols <- dendextend::get_leaves_branches_col(dend)
    dendextend::labels_colors(dend) <- lcols
    # dendextend::labels(dend) <- ct
    
    withr::with_par(list(mai = par("mai") + c(0, 0, 0, 0.5)),
    {
        plot(dend, ylab = "Tanimoto dist", ...)
        legend("topright", legend = seq_len(nclust),
               bty = "n", cex = 1, fill = cols, inset = c(-0.1, 0), xpd = NA,
               title = "cluster")
    })
})

setMethod("getMCS", "compoundsCluster", function(obj, groupName, cluster)
{
    ac <- checkmate::makeAssertCollection()
    assertChoiceSilent(groupName, names(obj@clusters), add = ac)
    
    cc <- obj@cutClusters[[groupName]]
    nclust <- length(unique(cc))
    checkmate::assertInt(cluster, lower = 0, upper = nclust)
    checkmate::reportAssertions(ac)

    mols <- obj@molecules[[groupName]][cc == cluster]
    mcons <- mols[[1]]
    if (length(mols) > 1)
    {
        for (i in seq(2, length(mols)))
            mcons <- rcdk::get.mcs(mcons, mols[[i]])
    }

    return(mcons)    
})

setMethod("plotStructure", "compoundsCluster", function(obj, groupName, cluster,
                                                        width = 500, height = 500)
{
    aapply(checkmate::assertNumber, . ~ width + height, lower = 0, finite = TRUE)
    rcdkplot(getMCS(obj, groupName, cluster), width, height)
})

setMethod("makeHCluster", "compounds", function(obj, method, fpType = "extended",
                                                fpSimMethod = "tanimoto",
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
        dist <- as.dist(1 - fingerprint::fp.sim.matrix(fps, method = fpSimMethod))
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
    
    return(compoundsCluster(clusters = clust, molecules = mols, cutClusters = cutClusters,
                            properties = list(method = method, fpType = fpType,
                                              fpSimMethod = fpSimMethod)))
})
