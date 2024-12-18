# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
#' @include compounds.R
NULL

#' Hierarchical clustering of compounds
#'
#' Perform hierarchical clustering of structure candidates based on chemical
#' similarity and obtain overall structural information based on the maximum
#' common structure (MCS).
#'
#' Often many possible chemical structure candidates are found for each feature
#' group when performing \link[=generateCompounds]{compound annotation}.
#' Therefore, it may be useful to obtain an overview of their general structural
#' properties. One strategy is to perform hierarchical clustering based on their
#' chemical (dis)similarity, for instance, using the Tanimoto score. The
#' resulting clusters can then be characterized by evaluating their
#' \emph{maximum common substructure} (MCS).
#'
#' @section Source: The methodology applied here has been largely derived from
#'   \file{chemclust.R} from the \pkg{metfRag} package and the package vignette
#'   of \CRANpkg{rcdk}.
#'
#' @name compounds-cluster
#' @seealso compoundsCluster
NULL

#' Compounds cluster class
#'
#' Objects from this class are used to store hierarchical clustering data of
#' candidate structures within \code{\link{compounds}} objects.
#'
#' Objects from this type are returned by the \code{compounds} method for
#' \code{\link[=makeHCluster,compounds-method]{makeHCluster}}.
#'
#' @slot clusters A \code{list} with \code{\link{hclust}} objects for each
#'   feature group.
#' @slot dists A \code{list} with distance matrices for each feature group.
#' @slot SMILES A \code{list} containing a vector with \code{SMILES} for all
#'   candidate structures per feature group.
#' @slot cutClusters A \code{list} with assigned clusters for all candidates per
#'   feature group (same format as what \code{\link{cutree}} returns).
#' @slot properties A list containing general properties and parameters used for
#'   clustering.
#'
#' @param obj,x,object A \code{compoundsCluster} object.
#' @param groupName A character specifying the feature group name.
#' @param cluster A numeric value specifying the cluster.
#' @param \dots Further arguments passed directly to the plotting function
#'   (\code{plot} or \code{\link{plot.dendrogram}}).
#'
#' @templateVar seli feature groups
#' @templateVar selOrderi groupNames()
#' @templateVar noextract TRUE
#' @template sub_sel_del-args
#'
#' @return \code{cutTree} and \code{cutTreeDynamic} return the modified
#'   \code{compoundsCluster} object.
#'
#' @export
compoundsCluster <- setClass("compoundsCluster",
                             slots = c(clusters = "list", dists = "list", SMILES = "list", cutClusters = "list",
                                       properties = "list"))

#' @describeIn compoundsCluster Accessor method to the \code{clusters} slot.
#'   Returns a list that contains for each feature group an object as returned
#'   by \code{\link{hclust}}.
#' @export
setMethod("clusters", "compoundsCluster", function(obj) obj@clusters)

#' @describeIn compoundsCluster Accessor method to the \code{cutClusters} slot.
#'   Returns a list that contains for each feature group a vector with cluster
#'   membership for each candidate (format as \code{\link{cutree}}).
#' @export
setMethod("cutClusters", "compoundsCluster", function(obj) obj@cutClusters)

#' @describeIn compoundsCluster Returns a list with properties on how the
#'   clustering was performed.
#' @export
setMethod("clusterProperties", "compoundsCluster", function(obj) obj@properties)

#' @templateVar class compoundsCluster
#' @templateVar what feature groups
#' @template strmethod
#' @export
setMethod("groupNames", "compoundsCluster", function(obj) names(obj@clusters))

#' @describeIn compoundsCluster Returns the total number of clusters.
#' @export
setMethod("length", "compoundsCluster", function(x) sum(lengths(x)))

#' @describeIn compoundsCluster Returns a \code{vector} with the number of
#'   clusters per feature group.
#' @param use.names A logical value specifying whether the returned vector
#'   should be named with the feature group names.
#' @export
setMethod("lengths", "compoundsCluster", function(x, use.names = TRUE)
{
    if (length(x@cutClusters) == 0)
        return(0)
    sapply(x@cutClusters, function(cc) length(unique(cc)), USE.NAMES = use.names)
})

#' @describeIn compoundsCluster Show summary information for this object.
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

#' @describeIn compoundsCluster Subset on feature groups.
#' @export
setMethod("[", c("compoundsCluster", "ANY", "missing", "missing"), function(x, i, ...)
{
    if (!missing(i))
    {
        i <- assertSubsetArgAndToChr(i, groupNames(x))
        x@clusters <- x@clusters[i]
        x@SMILES <- x@SMILES[i]
        x@cutClusters <- x@cutClusters[i]
    }
    
    return(x)
})

#' @describeIn compoundsCluster Manually (re-)cut a dendrogram that was
#'   generated for a feature group.
#' @param k,h Desired number of clusters or tree height to be used for cutting
#'   the dendrogram, respecitively. One or the other must be specified.
#'   Analogous to \code{\link{cutree}}.
#' @export
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

#' @describeIn compoundsCluster Automatically (re-)cut a dendrogram that was
#'   generated for a feature group using the \code{\link{cutreeDynamicTree}}
#'   function from \pkg{\link{dynamicTreeCut}}.
#' 
#' @template dynamictreecut
#' 
#' @export
setMethod("treeCutDynamic", "compoundsCluster", function(obj, maxTreeHeight, deepSplit,
                                                         minModuleSize, groupName)
{
    ac <- checkmate::makeAssertCollection()
    assertChoiceSilent(groupName, names(obj@clusters), add = ac)
    assertDynamicTreeCutArgs(maxTreeHeight, deepSplit, minModuleSize, ac)
    checkmate::reportAssertions(ac)
    
    obj@cutClusters[[groupName]] <- doDynamicTreeCut(obj@clusters[[groupName]], maxTreeHeight,
                                                     deepSplit, minModuleSize)
    return(obj)
})

#' @describeIn compoundsCluster Plot the dendrogram for clustered compounds of a
#'   feature group. Clusters are highlighted using \CRANpkg{dendextend}.
#' @templateVar withoutDots TRUE
#' @template plot_clust
#' @export
setMethod("plot", c(x = "compoundsCluster", y = "missing"), function(x, ..., groupName, pal = "Paired",
                                                                     colourBranches = lengths(x)[groupName] < 50,
                                                                     showLegend = lengths(x)[groupName] < 20)
{
    assertChoiceSilent(groupName, names(x@clusters))
    checkmate::assertString(pal, min.chars = 1)
    plotDendroWithClusters(as.dendrogram(x@clusters[[groupName]]), x@cutClusters[[groupName]], pal,
                           colourBranches, showLegend, ylab = "Tanimoto dist", ...)
    invisible(NULL)
})

setMethod("plotHash", "compoundsCluster", function(x, groupName,  ...)
{
    makeHash(x@clusters[[groupName]], x@cutClusters[[groupName]], ...)
})

#' @describeIn compoundsCluster Calculates the maximum common substructure (MCS)
#'   for all candidate structures within a specified cluster. This method uses
#'   the \code{\link{get.mcs}} function from \CRANpkg{rcdk}.
#' @return \code{getMCS} returns an \CRANpkg{rcdk} molecule object
#'   (\code{IAtomContainer}).
#' @export
setMethod("getMCS", "compoundsCluster", function(obj, groupName, cluster)
{
    assertChoiceSilent(groupName, names(obj@clusters))
    
    ac <- checkmate::makeAssertCollection()
    cc <- obj@cutClusters[[groupName]]
    nclust <- length(unique(cc))
    checkmate::assertInt(cluster, lower = 0, upper = nclust, add = ac)
    checkmate::reportAssertions(ac)

    mols <- getMoleculesFromSMILES(obj@SMILES[[groupName]][cc == cluster], doTyping = TRUE,
                                   emptyIfFails = TRUE)
    mcons <- mols[[1]]
    if (length(mols) > 1)
    {
        for (i in seq(2, length(mols)))
        {
            if (!isValidMol(mols[[i]]))
                return(emptyMol())
            
            # might fail if there is no overlap...
            tryCatch(mcons <- rcdk::get.mcs(mcons, mols[[i]]), error = function(e) FALSE)
            if (is.null(mcons) || mcons == FALSE)
                return(emptyMol())
        }
    }

    return(mcons)    
})

#' @describeIn compoundsCluster Plots the maximum common substructure (MCS) for
#'   all candidate structures within a specified cluster.
#'
#' @param width,height The dimensions (in pixels) of the raster image that
#'   should be plotted.
#' @param withTitle A logical value specifying whether a title should be added.
#' 
#' @export
setMethod("plotStructure", "compoundsCluster", function(obj, groupName, cluster,
                                                        width = 500, height = 500,
                                                        withTitle = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    assertChoiceSilent(groupName, names(obj@clusters), add = ac)
    aapply(checkmate::assertNumber, . ~ width + height, lower = 0, finite = TRUE,
           fixed = list(add = ac))
    checkmate::assertFlag(withTitle, add = ac)
    checkmate::reportAssertions(ac)
    
    plot(getRCDKStructurePlot(getMCS(obj, groupName, cluster), width, height))
    
    if (withTitle)
    {
        count <- sum(obj@cutClusters[[groupName]] == cluster)
        title(sprintf("cluster %d (n=%d)", cluster, count))
    }
})

setMethod("plotStructureHash", "compoundsCluster", function(obj, groupName, cluster,
                                                            width = 500, height = 500,
                                                            withTitle = TRUE)
{
    cc <- obj@cutClusters[[groupName]]
    ret <- makeHash(obj@SMILES[[groupName]][cc == cluster], width, height)
    if (withTitle)
        ret <- makeHash(ret, sum(obj@cutClusters[[groupName]] == cluster))
    return(ret)
})

#' @templateVar class compoundsCluster
#' @template plotsil
#' @export
setMethod("plotSilhouettes", "compoundsCluster", function(obj, kSeq, groupName, pch = 16, type = "b", ...)
{
    assertChoiceSilent(groupName, names(obj@clusters))
    checkmate::assertIntegerish(kSeq, lower = 2, any.missing = FALSE)
    doPlotSilhouettes(obj@clusters[[groupName]], obj@dists[[groupName]], kSeq, pch, type, ...)
    invisible(NULL)
})

#' @details \code{makeHCluster} performs hierarchical clustering of all
#'   structure candidates for each feature group within a
#'   \code{\link{compounds}} object. The resulting dendrograms are automatically
#'   cut using the \code{\link{cutreeDynamicTree}} function from the
#'   \pkg{\link{dynamicTreeCut}} package. The returned
#'   \code{\link{compoundsCluster}} object can then be used, for instance, for
#'   plotting dendrograms and MCS structures and manually re-cutting specific
#'   clusters.
#'
#' @param obj The \code{\link{compounds}} object to be clustered.
#' @param method The clustering method passed to \code{\link{hclust}}.
#'
#' @template fp-args
#' @template dynamictreecut
#'
#' @return \code{makeHCluster} returns an \code{\link{compoundsCluster}} object.
#' 
#' @references \addCitations{rcdk}{1}
#' 
#' @rdname compounds-cluster
#' @aliases makeHCluster
#' @export
setMethod("makeHCluster", "compounds", function(obj, method, fpType = "extended",
                                                fpSimMethod = "tanimoto",
                                                maxTreeHeight = 1, deepSplit = TRUE,
                                                minModuleSize = 1)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertString, . ~ method + fpType + fpSimMethod, min.chars = 1, fixed = list(add = ac))
    assertDynamicTreeCutArgs(maxTreeHeight, deepSplit, minModuleSize, ac)
    checkmate::reportAssertions(ac)

    hash <- makeHash(obj, method, fpType, fpSimMethod, maxTreeHeight, deepSplit, minModuleSize)
    cd <- loadCacheData("compoundsCluster", hash)
    if (!is.null(cd))
        return(cd)
        
    compTable <- annotations(obj)

    if (length(obj) == 0)
    {
        return(compoundsCluster(clusters = list(), SMILES = list(), cutClusters = list(),
                                properties = list(method = method, fpType = fpType,
                                                  fpSimMethod = fpSimMethod)))
    }

    mols <- sapply(compTable, function(ct) getMoleculesFromSMILES(ct$SMILES,
                                                                  doTyping = TRUE, emptyIfFails = TRUE),
                   simplify = FALSE) 
    
    cat("Calculating distances ...\n")
    prog <- openProgBar(0, length(mols))
    
    dists <- lapply(seq_along(mols), function(i)
    {
        if (length(mols[[i]]) < 2)
            return(NULL) # need multiple candidates to cluster
        
        fps <- lapply(mols[[i]], rcdk::get.fingerprint, type = fpType)
        dist <- as.dist(1 - fingerprint::fp.sim.matrix(fps, method = fpSimMethod))
        setTxtProgressBar(prog, i)
        return(dist)
    })
    hasDist <- !sapply(dists, is.null)
    dists <- dists[hasDist]
    setTxtProgressBar(prog, length(obj))
    close(prog)
    
    cat("Performing clustering ...\n")
    prog <- openProgBar(0, length(dists))
    
    clusts <- lapply(seq_along(dists), function(di)
    {
        hc <- hclust(dists[[di]], method)
        setTxtProgressBar(prog, di)
        return(hc)
    })
    
    setTxtProgressBar(prog, length(obj))
    close(prog)
    
    cat("Performing dynamic tree cutting ...\n")
    prog <- openProgBar(0, length(clusts))
    cutClusters <- lapply(seq_along(clusts), function(ci)
    {
        dendro <- clusts[[ci]]
        ret <- doDynamicTreeCut(dendro, maxTreeHeight, deepSplit, minModuleSize)
        setTxtProgressBar(prog, ci)
        return(ret)
    })

    setTxtProgressBar(prog, length(clusts))
    close(prog)
    
    gNames <- names(compTable)[hasDist]
    names(clusts) <- gNames
    names(dists) <- gNames
    names(cutClusters) <- gNames
    
    ret <- compoundsCluster(clusters = clusts, dists = dists,
                            SMILES = sapply(compTable, "[[", "SMILES", simplify = FALSE),
                            cutClusters = cutClusters,
                            properties = list(method = method, fpType = fpType,
                                              fpSimMethod = fpSimMethod))
    saveCacheData("compoundsCluster", ret, hash)
    return(ret)
})
