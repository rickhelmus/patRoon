# BUG: need another block (in this case @param) before brew code block to avoid useless warnings about missing title/description

#' @param <%=if (!exists("selj")) "drop,j" else "drop" %> ignored.
#' 
#' <% fmt <- "@param { vars } A numeric or character value which is used to select { sel } by
#' their index or name, respectively (for the order/names see \\code{{{ selOrder }}})." %>
#' 
#' <% if (!exists("noextract")) fmt <- paste0(fmt, "\\cr\\cr For \\code{{[}}:") %>
#' <% fmt <- paste(fmt, "Can also be logical to perform logical selection
#' (similar to regular vectors). If missing all { sel } are selected.") %>
#' <% if (!exists("noextract")) fmt <- paste0(fmt, "\\cr\\cr For \\code{{[[}}: should be a scalar value.") %>
#' 
#' <% optFmt <- "If \\code{{j}} is not specified, \\code{{i}} selects by { selj } instead." %>
#' <% if (exists("optionalji")) fmt <- paste(fmt, optFmt) %>
#' <% if (exists("optionalj")) fmt <- paste(fmt, "\\code{{j}} is optional.") %>
#' <% vars <- if (exists("selj")) "i,j" else "i" %>
#' <% allsel <- if (exists("selj")) paste(seli, selj, sep = "/") else seli %>
#' <% allSelOrder <- if (exists("selOrderj")) paste(selOrderi, selOrderj, sep = "/") else selOrderi %>
#' 
#' <%=glue::glue(fmt, vars = vars, selj = selj, sel = allsel, selOrder = allSelOrder) %>
#' 
#' <%=if (exists("dollarOpName")) glue::glue("@param name The { name } name (partially matched).", name = dollarOpName) %>
