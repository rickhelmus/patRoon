#' @param parents The parents for which transformation products should be obtained. This can be (1) a suspect list (see
#'   \link[=suspect-screening]{suspect screening} for more information), (2) the resulting output of
#'   \code{\link{screenSuspects}} or (3) a \code{\link{compounds}} annotation object. In the former two cases, the
#'   suspect (hits) are used as parents, whereas in the latter case all candidates are used as parents.
#'   <%= if (exists("parNULL")) "If \\code{NULL} then TPs for all parents in the library are obtained." %>
#' @param skipInvalid If set to \code{TRUE} then the parents will be skipped (with a warning) for which insufficient
#'   information (\emph{e.g.} SMILES) is available.
#'
#' @details An important advantage of this algorithm is that it provides structural information for generated TPs.
#'   However, this also means that if the input is from a parent suspect list or screening then either \acronym{SMILES}
#'   or \acronym{InChI} information must be available for the parents.
#'
#' @note When the \code{parents} argument is a \code{\link{compounds}} object, the candidate library \code{identifier}
#'   is used in case the candidate has no defined \code{compoundName}.
