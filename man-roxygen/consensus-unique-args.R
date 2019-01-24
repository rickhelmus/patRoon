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
