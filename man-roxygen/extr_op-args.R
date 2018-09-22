# BUG: need another block (in this case @param) before brew code block to avoid useless warnings about missing title/description

#' @describeIn <%=class%> Operator to subset on <%=whati%> <%=if (exists("whatj")) paste("and/or", whatj) %>.
#' 
#' <% fmt <- "Either a numeric, character or logical \\code{{vector}} that is used to select { what } by their index, name and logical selection, respectively (for the order/names see \\code{{{ order }}}). If missing all { what } are selected." %>
#' @param i <%=glue::glue(fmt, what = whati, order = orderi)%>
#' @param j <%=if (exists("whatj")) glue::glue(fmt, what = whatj, order = orderj) else "ignored." %>
#' @param drop ignored.

