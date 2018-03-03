#' @details \code{<%=func%>} is a generic function that will <%=what%>
#'   using one of the supported algorithms. The actual functionality is provided
#'   by algorithm specific functions such as \code{<%=ex1%>} and
#'   \code{<%=ex2%>}. While these functions may be called directly,
#'   \code{<%=func%>} provides a generic interface and is therefore usually
#'   preferred.
#'
#' @param algorithm A character string describing the algorithm that should be
#'   used: <%=paste0("\\code{\"", unlist(strsplit(algos, ",")), "\"}", collapse = ", ")%>
