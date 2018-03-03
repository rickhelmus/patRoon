#' @details \code{<%=func%>} <%=desc%>
#'
#' <% cl <- showMethods(func, where = "package:patRoon", printTo = FALSE) %>
#' <% cl <- cl[grepl("obj=", cl)] %>
#' <% cl <- sub("obj=\"([[:alpha:]]+)\"", "\\1", cl) %>
#'
#' \itemize{
#'     \item Methods are defined for: <%= { paste0(sprintf("\\code{\\link[=%s,%s-method]{%s}}", func, cl, cl), collapse = ", ") } %>.
#' }
#' @rdname generics
#' @aliases <%=func%>

