#' @return \code{plotVenn} (invisibly) returns a list with the following fields: \itemize{
#' \item \code{gList} the \code{gList} object that was returned by
#'   the utilized \pkg{\link{VennDiagram}} plotting function.
#' \item \code{areas} The total area for each plotted group.
#' \item \code{intersectionCounts} The number of intersections between groups.
#' }
#'
#' @return The order for the \code{areas} and \code{intersectionCounts} fields is the same as the parameter order
#' from the used plotting function (see \emph{e.g.} \code{\link{draw.pairwise.venn}} and
#' \code{\link{draw.triple.venn}}).
