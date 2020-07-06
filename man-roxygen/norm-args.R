#' @param <%=normParam%> A \code{character} that specifies how normalization of
#'   annotation scorings occurs. Either
#'   <%=if (!exists("noNone")) '\\code{"none"} (no normalization),' %>
#'   \code{"max"} (normalize to max value) or \code{"minmax"} (perform min-max
#'   normalization). Note that normalization of negative scores (e.g. output by
#'   \command{SIRIUS}) is always performed as min-max. Furthermore, currently
#'   normalization for \code{compounds} takes the original min/max scoring
#'   values into account when candidates were generated. Thus, for
#'   \code{compounds} scoring, normalization is not affected when candidate
#'   results were removed after they were generated (\emph{e.g.} by use of
#'   \code{filter}).
#'
#' <%=if (exists("excludeParam")) glue::glue("@param { excludeParam } A
#'   \\code{{character}} vector specifying any compound scoring names that
#'   should \\emph{{not}} be normalized. Set to \\code{{NULL}} to normalize all
#'   scorings. Note that whether any normalization occurs is set by the
#'   \\code{{{ excludeParam }}} argument.
#'
#'   For \\code{{compounds}}: By default \\code{{score}} and
#'   \\code{{individualMoNAScore}} are set to mimic the behavior of the
#'   \\command{{MetFrag}} web interface.", excludeParam = excludeParam) %>
