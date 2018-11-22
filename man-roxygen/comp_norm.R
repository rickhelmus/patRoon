#' @param <%=normParam%> A \code{character} that specifies how normalization of
#'   compound scorings occurs. Either \code{"none"} (no normalization),
#'   \code{"max"} (normalize to max value) or \code{"minmax"} (perform min-max
#'   normalization). Note that normalization of negative scores (e.g. output by
#'   SIRIUS) is always performed as min-max.
#' @param <%=excludeParam%> A \code{character} vector specifying any compound
#'   scoring names that should \emph{not} be normalized. Set to \code{NULL} to
#'   normalize all scorings. Note that whether any normalization occurs is set
#'   by the \code{<%=excludeParam%>} argument.
