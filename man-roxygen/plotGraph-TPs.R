#' @describeIn <%=class%> Plots an interactive hierarchy graph of the transformation products. The
#'   resulting graph can be browsed interactively and allows exploration of the different TP formation pathways.
#'   Furthermore, results from \link[=generateComponentsTPs]{TP componentization} can be used to match the hierarchy
#'   with screening results. The graph is rendered with \pkg{\link{visNetwork}}.
#'
#' @param which Either a \code{character} or \code{integer} vector with one or more names/indices of the parents to
#'   plot.
#' @param components If specified (\emph{i.e.} not \code{NULL}), a \code{\link{componentsTPs}} object that is used for
#'   matching the graph with screening results. The TPs that were found will be marked. See also the \code{prune} and
#'   \code{onlyCompletePaths} arguments.
#' @param prune If \code{TRUE} and \code{components} is set, then pathways without \emph{any} detected TPs are not shown
#'   (pruned). See also the \code{onlyCompletePaths} and \code{components} arguments.
#' @param onlyCompletePaths If \code{TRUE} and \code{components} is set, then only pathways are shown for which
#'   \emph{all} TPs were detected. See also the \code{prune} and \code{components} arguments.
#' @param width,height Passed to \code{\link{visNetwork}}.