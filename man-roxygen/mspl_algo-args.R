#' @param maxMSRtWindow Maximum chromatographic peak window used for spectrum averaging (in seconds, +/- retention
#'   time). If \code{NULL} all spectra from a feature will be taken into account. Lower to decrease processing time.
#' @param avgFGroupParams A \code{list} with parameters used for averaging of peak lists for feature groups. See
#'   \code{\link{getDefAvgPListParams}} for more details.
#' @param \dots For sets workflows: further arguments passed to the non-sets workflow method.
