#' @param ni Optional argument. An expression used for subsetting the analyses. The
#'   \link[=analysis-information]{analysis information} is first subset and the remaining rows are used to determine for
#'   which analyses the results should be kept. The unevaluated \code{ni} expression is used to set the \code{i}
#'   argument of the subset operator of \CRANpkg{data.table}, which therefore brings the advanced subsetting
#'   capabilities of \CRANpkg{data.table} (see the
#'   \code{\link[data.table]{data.table}} documentation for more details). For instance, \code{<%=ex%>[replicate ==
#'   "standard"]} would subset all analyses assigned with the replicate \code{"standard"}.
#' @param reorder If \code{TRUE} then the order of the analyses is changed to match the order of the \code{i} argument.
#'
#'   \setsWF If the \code{sets} argument is specified (and \code{i} is not) then the order of sets is changed instead.
#' @param drop Ignored.
