#' @include main.R
NULL

# UNDONE: add eawag references

# So it can be used as a S4 slot
setOldClass("hclust")
setOldClass("dissimilarity")

#' Hierarchical clustering information
#'
#' Objects from this class hold information required for
#' \link[=h-cluster]{hierarchical clustering}.
#'
#' Objects from this class are returned by \code{\link{makeHCluster}}.
#'
#' @param x,object,cInfo A \code{clusterInfo} object.
#'
#' @slot clusterm Numeric matrix with normalized feature group intensities that
#'   was used for clustering.
#' @slot distm Distance matrix that was used for clustering (obtained with
#'   \code{\link{daisy}}).
#' @slot clust Object returned by \code{\link{hclust}}.
#' @slot metric Distance metric that was used (\emph{e.g.} euclidean).
#' @slot method Clustering method that was used.
#' @slot average If \code{TRUE} then feature groups from the same replicate
#'   groups were averaged prior to clustering.
#'
#' @export
clusterInfo <- setClass("clusterInfo",
                        slots = c(clusterm = "matrix", distm = "dissimilarity", clust = "hclust", metric = "character",
                                  method = "character", average = "logical"))

#' Hierarchical clustering of feature groups
#'
#' Utilities to perform and process hierarchical clustering with feature groups.
#'
#' These functions provide an interface to common hierarchical clustering tools
#' such as \code{\link{hclust}} that can be used for clustering feature groups.
#'
#' @param fGroups A \code{\link{featureGroups}} object.
#' @param cInfo A \code{\link{clusterInfo}} object with relevant clustering
#'   information.
#' @param col Colours used for drawing the heatmap. See \code{\link{heatmap}} /
#'   \code{\link{d3heatmap}}.
#' @param \dots Further options passed to \code{\link{heatmap}} /
#'   \code{\link{d3heatmap}} (\code{drawHeatMap}), \code{\link{barplot}}
#'   (\code{plot}, \code{\link{silhouetteInfo}} signature) or
#'   \code{\link[graphics]{plot}} (\code{plot}, \code{\link{clusterInfo}} signature and
#'   \code{plotInt}).
#' @param k The desired number of branches from a given cluster object.
#' @param c The numeric identifier of the desired branch from a cut cluster.
#' @param x,obj A \code{\link{clusterInfo}}/\code{\link{silhouetteInfo}} object.
#' @name h-cluster
NULL

#' @details \code{makeHCluster} generates hierarchical clustering information
#'   for feature groups. Intensity data will optionally be averaged and then
#'   normalized. A distance matrix is calculated with \code{\link{daisy}} and
#'   clustering is performed with \code{\link{hclust}}.
#'
#' @param normFunc Function that should be used for normalization of data.
#'   Intensitity values of a feature group will be divided by the result of this
#'   function when it is called with all intensity values of that feature group.
#'   For example, when \code{\link{max}} is used normalized intensities will be
#'   between zero and one.
#' @param metric Distance metric used to calculate the distance matrix (passed
#'   to \code{\link{daisy}}).
#' @param method Clustering method that should be applied (passed to
#'   \code{\link{hclust}}).
#' @param average If \code{TRUE} then all intensitity data will be averaged for
#'   each replicate group.
#'
#' @return \code{makeHCluster} Returns an object of the
#'   \code{\link{clusterInfo}} class.
#'
#' @rdname h-cluster
#' @aliases makeHCluster
#' @export
setMethod("makeHCluster", "featureGroups", function(fGroups, normFunc = max, metric, method, average)
{
    if (average)
    {
        gTable <- averageGroups(fGroups)
        analysis <- unique(analysisInfo(fGroups)$group)
    }
    else
    {
        gTable <- groups(fGroups)
        analysis <- analysisInfo(fGroups)$analysis
    }

    clusterdt <- copy(gTable)

    cat("Normalizing data... ")
    normv <- clusterdt[, lapply(.SD, normFunc)]
    for (g in seq_along(clusterdt))
        set(clusterdt, j = g, value = clusterdt[[g]] / normv[[g]])
    cat("Done!\n")

    clusterm <- as.matrix(transpose(clusterdt))
    rownames(clusterm) <- colnames(gTable)
    colnames(clusterm) <- analysis

    cat("Calculating distance matrix... ")
    distm <- daisy(clusterm, metric)
    cat("Done!\n")

    cat("Hierarchical clustering... ")
    clust <- hclust(distm, method)
    cat("Done!\n")

    return(clusterInfo(clusterm = clusterm, distm = distm, clust = clust, metric = metric, method = method, average = average))
})

#' @details \code{hClusterFilter} isolates feature groups within a
#'   \code{\link{featureGroups}} object which are present within a given branch
#'   of a cluster object.
#' @rdname h-cluster
#' @aliases hClusterFilter
#' @export
setMethod("hClusterFilter", c("featureGroups", "clusterInfo"), function(fGroups, cInfo, k, c)
{
    ct <- cutree(cInfo@clust, k)
    return(fGroups[, colnames(groups(fGroups)) %in% names(ct)[ct == c]])
})


#' @describeIn clusterInfo \code{length} returns the total number of clusters.
#' @export
setMethod("length", "clusterInfo", function(x) length(x@clust$height))

#' @describeIn clusterInfo \code{show} prints general information about this
#'   object.
#' @export
setMethod("show", "clusterInfo", function(object)
{
    printf("A clusterInfo object (%s)\n", class(object))

    printf("Number of clusters: %d\n", length(object))

    cp <- clusterProperties(object)
    printf("Properties:\n")
    printf(" - metric: %s\n", cp$metric)
    printf(" - method: %s\n", cp$method)
    printf(" - data averaged: %s\n", cp$average)

    showObjectSize(object)
})

#' @describeIn clusterInfo returns a list with data and properties relevant to
#'   the clustering that was performed (see \verb{Slots} section for a
#'   description).
#' @return \code{clusterProperties} returns a \code{list}
#' @export
#' @aliases clusterProperties
setMethod("clusterProperties", "clusterInfo", function(cInfo)
{
    return(list(clusterm = cInfo@clusterm, distm = cInfo@distm, clust = cInfo@clust,
                metric = cInfo@metric, method = cInfo@method, average = cInfo@average))
})


#' @details \code{drawHeatMap} draws a heatmap using the \code{\link{heatmap}}
#'   or \code{\link{d3heatmap}} function.
#' @param interactive If \code{TRUE} an interactive heatmap will be drawn (with
#'   \code{\link{d3heatmap}}).
#' @return \code{drawHeatMap} returns the same as \code{\link{heatmap}} or
#'   \code{\link{d3heatmap}}.
#' @rdname h-cluster
#' @aliases drawHeatMap
#' @export
setMethod("drawHeatMap", "clusterInfo", function(cInfo, col, interactive, ...)
{
    if (interactive)
        d3heatmap::d3heatmap(cInfo@clusterm, Colv = NA, distfun = function(d) dist(d, cInfo@metric), hclustfun = function(h) hclust(h, cInfo@method),
                             scale = "none", colors = col, ...)
    else
        heatmap(cInfo@clusterm, Colv = NA, distfun = function(d) dist(d, cInfo@metric), hclustfun = function(h) hclust(h, cInfo@method),
                scale = "none", col = col, ...)
})

#' @details \code{getSilhouetteInfo} is used to obtain silhouette information of a
#'   \code{\link{clusterInfo}} object which is 'cut' by a given sequence of
#'   desired number of branches.
#' @param ranges An integer vector containing a sequence of desired amounts of
#'   branches to which silhouette information should be calculated.
#' @return \code{getSilhouetteInfo} returns a \code{\link{silhouetteInfo}}
#'   object.
#' @rdname h-cluster
#' @aliases getSilhouetteInfo
#' @export
setMethod("getSilhouetteInfo", "clusterInfo", function(cInfo, ranges)
{
    silInfo <- vector("list", length(seq))
    maxmw <- maxk <- NULL

    for (i in seq_along(ranges))
    {
        k <- ranges[i]
        ct <- cutree(cInfo@clust, k)
        sil <- silhouette(ct, cInfo@distm)
        sm <- summary(sil)
        silInfo[[i]] <- list(k = k, meanw = sm$avg.width, cmeans = sm$clus.avg.widths)

        if (is.null(maxmw) || maxmw < sm$avg.width)
        {
            maxmw <- sm$avg.width
            maxk <- k
        }
    }

    return(silhouetteInfo(ranges = ranges, maxk = maxk, maxmw = maxmw, silInfo = silInfo))
})

#' @details \code{plotInt} makes a plot for all (normalized) intensity profiles
#'   of the feature groups within a given cluster.
#' @rdname h-cluster
#' @export
setMethod("plotInt", "clusterInfo", function(obj, k, c, ...)
{
    ct <- cutree(obj@clust, k)
    plotm <- obj@clusterm[rownames(obj@clusterm) %in% names(ct)[ct == c], ]
    nsamp <- ncol(plotm)

    plot(x = c(0, nsamp), y = c(0, max(plotm)), type = "n", xlab = "", ylab = "normalized intensity", xaxt = "n", ...)
    axis(1, seq_len(nsamp), colnames(plotm), las = 2)

    px <- seq_len(nsamp)
    for (i in seq_len(nrow(plotm)))
        lines(x = px, y = plotm[i, ])
})

#' @details \code{plot} (\code{\link{clusterInfo}} signature) generates a dendrogram
#'   from a given cluster object and optionally highlights resulting branches
#'   when the cluster is cut.
#' @rdname h-cluster
#' @export
setMethod("plot", "clusterInfo", function(x, k = NULL, ...)
{
    plot(x@clust, ...)
    if (!is.null(k))
        rect.hclust(x@clust, k)
})

#' Silhouette information
#'
#' Objects from this class hold silhouette information from a
#' \code{\link{clusterInfo}} object.
#'
#' Objects from this class are returned by the \code{\link{getSilhouetteInfo}}
#' method.
#'
#' @slot ranges Numeric vector of \code{k} values that were used to cut the
#'   cluster and generate silhouette information from it.
#' @slot maxk \code{k} value (from \code{ranges} parameter) that generated the
#'   largest average silhouette width.
#' @slot maxmw Largest average silhouette width from all tested \code{ranges}.
#' @slot silInfo A \code{list} containing for each test iteration the fields
#'   \code{k} (the tested \code{k} value), \code{meanw} (mean of all silhouette
#'   widths) and \code{cmeans} (mean silhouette width for each cluster).
#'
#' @export
silhouetteInfo <- setClass("silhouetteInfo",
                           slots = c(ranges = "vector", maxk = "numeric", maxmw = "numeric", silInfo = "list"))

#' @details \code{plot} (\code{\link{silhouetteInfo}} signature) plots the various mean
#'   silhouette widths stored inside a \code{\link{silhouetteInfo}} object. The
#'   resulting barplot can be used to assess the optimal number of clusters.
#' @param y Ignored.
#' @return \code{plot} (\code{\link{silhouetteInfo}} signature) returns the result of
#'   \code{\link{barplot}}.
#' @rdname h-cluster
#' @export
setMethod("plot", "silhouetteInfo", function(x, ...)
{
    meanws <- sapply(x@silInfo, function(s) s$meanw)
    barplot(meanws, names.arg = x@ranges, xlab = "cluster count", ylab = "mean silhouette width", ...)
})
