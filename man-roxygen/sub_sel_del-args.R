# BUG: need another block (in this case @param) before brew code block to avoid useless warnings about missing title/description

#' @param <%=if (!exists("selj") && !exists("del")) "drop,j" else "drop" %> ignored.
#' 
#' <% fmt <- "@param { vars } For { ops }: A numeric or character value which is used to select { sel } by
#' their index or name, respectively (for the order/names see \\code{{{ selOrder }}})." %>
#' <% ops <- if (!exists("noextract")) "\\code{[}/\\code{[[}" else "\\code{[}" %>
#' <% if (!exists("noextract") || exists("del")) fmt <- paste0(fmt, "\\cr\\cr For \\code{{[}}:") %>
#' <% fmt <- paste(fmt, "Can also be logical to perform logical selection
#' (similar to regular vectors). If missing all { sel } are selected.") %>
#' <% if (!exists("noextract")) fmt <- paste0(fmt, "\\cr\\cr For \\code{{[[}}: should be a scalar value.") %>
#' 
#' <% optFmt <- "If \\code{{j}} is not specified, \\code{{i}} selects by { selj } instead." %>
#' <% if (exists("optionalji")) fmt <- paste(fmt, optFmt) %>
#' <% if (exists("optionalj")) fmt <- paste(fmt, "\\code{{j}} is optional.") %>
#' 
#' <% if (exists("del")) fmt <- paste0(fmt, "\\cr\\cr For \\code{{delete}}: The data to remove from. \\code{{i}} are the
#' { deli } as numeric index, logical or character, \\code{{j}} the { delj } as { deljtype }. If either is
#' \\code{{NULL}} then data for all is removed. \\code{{j}} may also be a function: it will be called for each 
#' { delfwhat }, with { delfa1 } as first argument, { delfa2 } as second argument, and any other arguments passed as
#' \\code{{\\dots}} to \\code{{delete}}. The return value of this function specifies { delfr }.") %>
#'
#' <% vars <- if (exists("selj") || exists("del")) "i,j" else "i" %>
#' <% allsel <- if (exists("selj")) paste(seli, selj, sep = "/") else seli %>
#' <% allSelOrder <- if (exists("selOrderj")) paste(selOrderi, selOrderj, sep = "/") else selOrderi %>
#' 
#' <%=glue::glue(fmt, vars = vars, ops = ops, selj = selj, sel = allsel, selOrder = allSelOrder, delj = delj, delfa1 = delfa1, delfa2 = delfa2, delfr = delfr) %>
#' 
#' <%=if (exists("dollarOpName")) glue::glue("@param name The { name } name (partially matched).", name = dollarOpName) %>
