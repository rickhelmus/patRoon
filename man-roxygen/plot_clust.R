#' @param pal Colour palette to be used from \pkg{\link{RColorBrewer}}.
#' @param colourBranches Whether branches from cut clusters (and their labels)
#'   should be coloured. Might be slow with large numbers of clusters, hence,
#'   the default is only \code{TRUE} when this is not the case.
#' @param showLegend If \code{TRUE} and \code{colourBranches} is also
#'   \code{TRUE} then a legend will be shown which outlines cluster numbers and
#'   their colours. By default \code{TRUE} for small amount of clusters to avoid
#'   overflowing the plot.
#' @param \dots Any arguments directly given to \code{\link{plot.dendrogram}}.
