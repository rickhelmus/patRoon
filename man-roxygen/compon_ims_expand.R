#' @section IMS workflows: In IMS workflows with post mobility assignment (see
#'   \code{\link[=assignMobilities_feat]{assignMobilities}}) it may be necessary to call \code{\link{expandForIMS}} when
#'   componentization was performed \emph{prior} to mobility assignments, see its documentation for more details.
#'
#'   If mobilities were already assigned prior to componentization, then the \code{IMS} argument selects which feature
#'   groups are subjected to componentization. Data for mobility feature groups that were not considered (\emph{i.e.}
#'   when \code{IMS} is \code{FALSE} or \code{"maybe"}), will be expanded similarly as is done by
#'   \code{\link{expandForIMS}}.
