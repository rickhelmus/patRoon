#' @include main.R
#' @include TP.R
NULL

#' Base transformation products (TP) class with formula information
#'
#' Holds information for all TPs for a set of parents, including chemical formulae.
#'
#' This (virtual) class is derived from the \code{\link{transformationProducts}} base class, please see its
#' documentation for more details. Objects from this class are returned by \link[=generateTPs]{TP generators}. More
#' specifically, algorithms that works with chemical structures (\emph{e.g.} \code{biotransformer}), uses this class to
#' store their results. The methods defined for this class extend the functionality for the base
#' \code{\link{transformationProducts}} class.
#'
#' @param obj,TPs \code{transformationProductsFormula} derived object to be accessed
#' @param commonParents Only consider TPs from parents that are common to all compared objects.
#' @param \dots For \code{filter}: Further argument passed to the base
#'   \code{\link[=filter,transformationProducts-method]{filter method}}.
#'
#'   For \code{plotVenn}, \code{plotUpSet} and \code{consensus}: further (unique) \code{transformationProductsFormula}
#'   objects.
#'
#' @section Comparison between objects: The methods that compare different objects (\emph{e.g.} \code{plotVenn} and
#'   \code{consensus}) use the \acronym{InChIKey} to match TPs between objects. Moreover, the parents between objects
#'   are matched by their name. Hence, it is \emph{crucial} that the input parents to \code{\link{generateTPs}}
#'   (\emph{i.e.} the \code{parents} argument) are named equally.
#'
#' @seealso The base class \code{\link{transformationProducts}} for more relevant methods and \code{\link{generateTPs}}
#'
#' @templateVar class transformationProductsFormula
#' @template class-hierarchy
#'
#' @export
transformationProductsFormula <- setClass("transformationProductsFormula",
                                          contains = c("VIRTUAL", "transformationProducts"))

setMethod("linkParentsToFGroups", "transformationProductsFormula", function(TPs, fGroups)
{
    return(screenInfo(fGroups)[name %in% names(TPs), c("name", "group"), with = FALSE])
})

#' @describeIn transformationProductsFormula Plots an interactive hierarchy graph of the transformation products. The
#'   resulting graph can be browsed interactively and allows exploration of the different TP formation pathways.
#'   Furthermore, results from \link[=generateComponentsTPs]{TP componentization} can be used to match the hierarchy
#'   with screening results. The graph is rendered with \pkg{\link{visNetwork}}.
#'
#' @param which Either a \code{character} or \code{integer} vector with one or more names/indices of the parents to
#'   plot.
#' @param components If specified (\emph{i.e.} not \code{NULL}), a \code{\link{componentsTPs}} object that is used for
#'   matching the graph with screening results. The TPs that were found will be marked. See also the \code{prune} and
#'   \code{onlyCompletePaths} arguments.
#' @param structuresMax An \code{integer} with the maximum number of structures to plot. Setting a maximum is mainly
#'   done to avoid long times needed to construct the graph.
#' @param prune If \code{TRUE} and \code{components} is set, then pathways without \emph{any} detected TPs are not shown
#'   (pruned). See also the \code{onlyCompletePaths} and \code{components} arguments.
#' @param onlyCompletePaths If \code{TRUE} and \code{components} is set, then only pathways are shown for which
#'   \emph{all} TPs were detected. See also the \code{prune} and \code{components} arguments.
#'
#' @template plotGraph
#'
#' @export
setMethod("plotGraph", "transformationProductsFormula", function(obj, which, components = NULL, prune = TRUE,
                                                                 onlyCompletePaths = FALSE)
{
    checkmate::assert(
        checkmate::checkSubset(which, names(obj), empty.ok = FALSE),
        checkmate::checkIntegerish(which, lower = 1, upper = nrow(parents(obj)), any.missing = FALSE, min.len = 1),
        .var.name = "which"
    )
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(components, "componentsTPs", null.ok = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ prune + onlyCompletePaths, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    doPlotTPGraph(as.data.table(obj[which]), parents(obj),
                  cmpTab = if (!is.null(components)) as.data.table(components) else NULL, structuresMax = 0,
                  prune = prune, onlyCompletePaths = onlyCompletePaths)
})
