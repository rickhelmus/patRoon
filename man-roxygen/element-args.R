#' @param elements Only retain candidate formulae (neutral form) that match a
#'   given elemental restriction. The format of \code{elements} is a
#'   \code{character} string with elements that should be present where each
#'   element is followed by a valid amount or a range thereof. If no number is
#'   specified then \samp{1} is assumed. For instance,
#'   \code{elements="C1-10H2-20O0-2P"}, specifies that \samp{1-10}, \samp{2-20},
#'   \samp{0-2} and \samp{1} carbon, hydrogen, oxygen and phosphorus atoms
#'   should be present, respectively. When \code{length(elements)>1} formulas
#'   are tested to follow at least one of the given elemental restrictions. For
#'   instance, \code{elements=c("P", "S")} specifies that either one phosphorus
#'   or one sulphur atom should be present. Set to \code{NULL} to ignore this
#'   filter.
#' @param fragElements,lossElements Specifies elemental restrictions for
#'   fragment or neutral loss formulae (charged form). Candidates are retained
#'   if at least one of the fragment formulae follow the given restrictions. See
#'   \code{elements} for the used format.
