#' @param adduct An \code{\link{adduct}} object (or something that can be converted to it with \code{\link{as.adduct}}).
#'   Examples: \code{"[M-H]-"}, \code{"[M+Na]+"}. <%= if (!exists("plain")) "If the \\code{featureGroups} object has
#'   adduct annotations then these are used if \\code{adducts=NULL}." %>
#'
#'   <%= if (!exists("plain")) "\\setsWF The \\code{adduct} argument is not supported for sets workflows, since the
#'   adduct annotations will then always be used." %>
