#' @param absMinAbundance,relMinAbundance Minimum absolute or relative
#'   (\samp{0-1}) abundance across objects for a result to be kept. For
#'   instance, \code{relMinAbundance=0.5} means that a result should be present
#'   in at least half of the number of compared objects. Set to \samp{NULL} to
#'   ignore and keep all results. Limits cannot be set when \code{uniqueFrom} is
#'   not \code{NULL}.
#' @param uniqueFrom Set this argument to only retain <%=what%> that are unique
#'   within one or more of the objects for which the consensus is made.
#'   Selection is done by setting the value of \code{uniqueFrom} to a
#'   \code{logical} (values are recycled), \code{numeric} (select by index) or a
#'   \code{character} (as obtained with \code{algorithm(obj)}). For
#'   \code{logical} and \code{numeric} values the order corresponds to the order
#'   of the objects given for the consensus. Set to \code{NULL} to ignore.
#' @param uniqueOuter If \code{uniqueFrom} is not \code{NULL} and if
#'   \code{uniqueOuter=TRUE}: only retain data that are also unique between
#'   objects specified in \code{uniqueFrom}.
