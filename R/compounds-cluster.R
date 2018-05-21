#' @include main.R
#' @include compounds.R
NULL

doDynamicTreeCut <- function(dendro, maxTreeHeight, deepSplit,
                             minModuleSize)
{
    if (minModuleSize == 1)
    {
        # workaround adapted from RAMClustR (ramclustR.R)
        ret <- dynamicTreeCut::cutreeDynamicTree(dendro = dendro, maxTreeHeight = maxTreeHeight,
                                                 deepSplit = deepSplit, minModuleSize = 2)
        single <- which(ret == 0) # all unassigned have length = 1
        ret[single] <- max(ret) + seq_len(length(single))
    }
    else
        ret <- dynamicTreeCut::cutreeDynamicTree(dendro = dendro, maxTreeHeight = maxTreeHeight,
                                                 deepSplit = deepSplit, minModuleSize = minModuleSize)
    
    return(ret)
}

# This way of compound clustering largely based on metfRag's chemclust.R and the
# package vignette of rcdk

compoundsCluster <- setClass("compoundsCluster",
                             slots = c(clusters = "list", molecules = "list", cutClusters = "list",
                                       properties = "list"),
                             prototype = list(clusters = list(), molecules = list(), cutClusters = list(),
                                              properties = list()))

setMethod("clusters", "compoundsCluster", function(obj) obj@clusters)

setMethod("molecules", "compoundsCluster", function(obj) obj@molecules)

setMethod("cutClusters", "compoundsCluster", function(obj) obj@cutClusters)

setMethod("clusterProperties", "compoundsCluster", function(obj) obj@properties)

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

setMethod("treeCut", "compoundsCluster", function(obj, k = NULL, h = NULL, groupName)
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

setMethod("treeCutDynamic", "compoundsCluster", function(obj, maxTreeHeight, deepSplit,
                                                         minModuleSize, groupName)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertNumber(maxTreeHeight, 0, finite = TRUE, add = ac)
    checkmate::assertFlag(deepSplit, add = ac)
    checkmate::assertCount(minModuleSize, positive = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    obj@cutClusters[[groupName]] <- doDynamicTreeCut(obj@clusters[[groupName]], maxTreeHeight,
                                                     deepSplit, minModuleSize)
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
    dend <- dendextend::color_branches(dend, clusters = ct, col = cols[unique(ct)]) # unique: fixup colour order
    lcols <- dendextend::get_leaves_branches_col(dend)
    dendextend::labels_colors(dend) <- lcols
    
    withr::with_par(list(mar = c(1, 4, 0, 5.5)),
    {
        plot(dend, ylab = "Tanimoto dist", ...)
        legend("topright", legend = seq_len(nclust),
               bty = "n", cex = 1, fill = cols, inset = c(-0.18, 0), xpd = NA,
               ncol = 2, title = "cluster")
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
        {
            # might fail if there is no overlap...
            tryCatch(mcons <- rcdk::get.mcs(mcons, mols[[i]]), error = function(e) FALSE)
            if (mcons == FALSE)
                return(rcdk::parse.smiles("")) # return empty molecule
        }
    }

    return(mcons)    
})

setMethod("plotStructure", "compoundsCluster", function(obj, groupName, cluster,
                                                        width = 500, height = 500,
                                                        withTitle = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertNumber, . ~ width + height, lower = 0, finite = TRUE,
           fixed = list(add = ac))
    checkmate::assertFlag(withTitle, add = ac)
    checkmate::reportAssertions(ac)
    
    rcdkplot(getMCS(obj, groupName, cluster), width, height)
    
    if (withTitle)
    {
        count <- sum(obj@cutClusters[[groupName]] == cluster)
        title(sprintf("cluster %d (n=%d)", cluster, count))
    }
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
    

    cat("Performing clustering ...\n")    
    prog <- txtProgressBar(0, length(mols), style = 3)
    clust <- lapply(seq_along(mols), function(i)
    {
        for (j in seq_along(mols[[i]]))
        {
            rcdk::do.typing(mols[[i]][[j]])
            rcdk::do.aromaticity(mols[[i]][[j]])
        }
        
        fps <- lapply(mols[[i]], rcdk::get.fingerprint, type = fpType)
        dist <- as.dist(1 - fingerprint::fp.sim.matrix(fps, method = fpSimMethod))
        
        hc <- hclust(dist, method)
        setTxtProgressBar(prog, i)
        
        return(hc)
    })
    
    setTxtProgressBar(prog, length(obj))
    close(prog)
    
    cat("Performing dynamic tree cutting ...\n")
    prog <- txtProgressBar(0, length(obj), style = 3)
    cutClusters <- lapply(seq_along(clust), function(ci)
    {
        dendro <- clust[[ci]]
        ret <- doDynamicTreeCut(dendro, maxTreeHeight, deepSplit, minModuleSize)
        setTxtProgressBar(prog, i)
        return(ret)
    })

    setTxtProgressBar(prog, length(obj))
    close(prog)
    
    names(clust) <- names(compTable)
    names(cutClusters) <- names(compTable)
    
    return(compoundsCluster(clusters = clust, molecules = mols, cutClusters = cutClusters,
                            properties = list(method = method, fpType = fpType,
                                              fpSimMethod = fpSimMethod)))
})
