# automatically retrieve defined methods for a generic and create document
# links. This only works if the arguments of the method are named obj or objX.

#' @details \code{<%=func%>} <%=desc%>
#'
#' <% cl <- showMethods(func, where = "package:patRoon", printTo = FALSE) %>
#' <% cl <- cl[grepl("obj.?=", cl)] %>
#' <% cl <- gsub("obj.?=\"([[:alpha:]]+)\"", "\\1", cl) %>
#' <% cl <- gsub(" ", "", cl) %>
#'
#' \itemize{
#'     \item Methods are defined for: <%= { paste0(sprintf("\\code{\\link[=%s,%s-method]{%s}}", func, cl, cl), collapse = ", ") } %>.
#' }
#' @rdname generics
#' @aliases <%=func%>

