#' @details \code{makeSet} is used to initiate a \link[=sets-workflow]{sets workflow}. See the \verb{Sets workflows}
#'   section for more details.
#'
#' @param adducts The adduct assignments to each set. Should either be a \code{list} with \code{\link{adduct}} objects
#'   or a \code{character} vector (\emph{e.g.} \code{"[M+H]+"}). The order should follow that of the objects given to
#'   the \code{obj} and \code{\dots} arguments. <%= if (exists("adductNULL")) "if \\code{{NULL}} then adduct annotations
#'   are used." %>
#' @param labels The labels, or \emph{set names}, for each set to be created. The order should follow that of the
#'   objects given to the \code{obj} and \code{\dots} arguments. If \code{NULL}, then labels are automatically generated
#'   from the polarity of the specified \code{adducts} argument (\emph{e.g.} \code{"positive"}, \code{"negative"}).
