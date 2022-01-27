#' @param errorRetries Maximum number of retries after an error occurred. This may be useful to handle e.g. connection
#'   errors.
#' @param topMost Only keep this number of candidates (per feature group) with highest score. Set to \code{NULL} to
#'   always keep all candidates, however, please note that this may result in significant usage of CPU/RAM resources for
#'   large numbers of candidates.
#' @param \dots \setsWF Further arguments passed to the non-sets workflow method.
