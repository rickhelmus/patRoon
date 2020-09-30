#' @include main.R
#' @include components-set.R
NULL

componentsClustSet <- setClass("componentsClustSet", contains = "componentsSet")

# UNDONE: more methods? or document clearly that setObjects should be obtained for that?

#' @describeIn componentsClustSet Manually (re-)cut the dendrogram.
#' @param k,h Desired number of clusters or tree height to be used for cutting
#'   the dendrogram, respectively. One or the other must be specified.
#'   Analogous to \code{\link{cutree}}.
#' @export
setMethod("treeCut", "componentsClustSet", function(obj, k = NULL, h = NULL)
{
    obj@setObjects <- lapply(setObjects(obj), treeCut, k = k, h = h)
    return(syncComponentsSetObjects(obj))
})

#' @describeIn componentsClustSet Automatically (re-)cut the dendrogram using
#'   the \code{\link{cutreeDynamicTree}} function from
#'   \pkg{\link{dynamicTreeCut}}.
#'
#' @template dynamictreecut
#'
#' @export
setMethod("treeCutDynamic", "componentsClustSet", function(obj, maxTreeHeight, deepSplit,
                                                           minModuleSize)
{
    obj@setObjects <- lapply(setObjects(obj), treeCutDynamic, maxTreeHeight = maxTreeHeight,
                             deepSplit = deepSplit)
    return(syncComponentsSetObjects(obj))
})

#' @describeIn componentsClustSet generates a dendrogram from a given cluster
#'   object and optionally highlights resulting branches when the cluster is
#'   cut.
#' @param numericLabels Set to \code{TRUE} to label with numeric indices instead
#'   of (long) feature group names.
#' @templateVar withoutDots TRUE
#' @template plot_clust
#' @export
setMethod("plot", "componentsClustSet", function(x, ..., set)
{
    checkmate::assertChoice(set, sets(x))
    plot(setObjects(x)[[set]], ...)
})

#' @templateVar class componentsClustSet
#' @template plotsil
#' @export
setMethod("plotSilhouettes", "componentsClustSet", function(obj, ..., set)
{
    checkmate::assertChoice(set, sets(x))
    plotSilhouettes(setObjects(x)[[set]], ...)
})
