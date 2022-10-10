#' @param transformations A \code{data.frame} with transformation reactions to be used for calculating the TPs (see
#'   details below). If \code{NULL}, a default table from Schollee \emph{et al.} is used (see references).
#'
#' @section Transformation reactions: The \code{transformations} argument specifies custom rules to calculate
#'   transformation products. This should be a \code{data.frame} with the following columns: \itemize{
#'
#'   \item \code{transformation} The name of the chemical transformation
#'
#'   \item \code{add} The elements that are added by this reaction (\emph{e.g.} \code{"O"}).
#'
#'   \item \code{sub} The elements that are removed by this reaction (\emph{e.g.} \code{"H2O"}).
#'
#'   \item \code{retDir} The expected retention time direction relative to the parent (assuming a reversed phase like LC
#'   separation). Valid values are: \samp{-1} (elutes before the parent), \samp{1} (elutes after the parent) or \samp{0}
#'   (no significant change or unknown).
#'
#'   }
#'
#' @section Source: The algorithms using transformation reactions are directly based on the work done by Schollee
#'   \emph{et al.} (see references).
#'
#' @references \insertRef{Scholle2015}{patRoon}
