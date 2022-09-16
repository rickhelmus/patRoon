#' @param generations An \code{integer} that specifies the number of transformation generations. TPs for subsequent
#'   iterations obtained by repeating the library search where the TPs from the previous generation are considered
#'   parents.
#' @param matchParentsBy A \code{character} that specifies how the input parents are matched with the data from the TP
#'   library. Valid options are: \code{"InChIKey"}, \code{"InChIKey1"}, \code{"InChI"}, \code{"SMILES"},
#'   \code{"formula"}, \code{"name"}. If the parent from the TP library is matched with multiple input parents then only
#'   the first is considered.
#' @param matchGenerationsBy Similar to \code{matchParentsBy}, but specifies how parents/TPs are matched when
#'   \code{generations>1}.
#'
#' @section TP libraries: The \code{TPLibrary} argument is used to specify a custom TP library. This should be a
#'   \code{data.frame} where each row specifies a TP for a parent, with the following columns: \itemize{
#'
#'   \item \code{parent_name} and \code{TP_name}: The name of the parent/TP.
#'
#'   \item \code{parent_<%=id%>} and \code{TP_<%=id%>} The <%=id%> of the parent/TP structure.
#'
#'   \item \code{retDir} The retention direction of the TP compared to its parent: \samp{-1} (elutes before), \samp{1}
#'   (elutes after) or \samp{0} (elutes similarly or unknown). If not specified then the \code{log P} values below may
#'   be used to calculate retention time directions. (\strong{optional})
#'
#'   \item \code{parent_LogP} and \code{TP_LogP} The \code{log P} values for the parent/TP. (\strong{optional})
#'
#'   \item \code{LogPDiff} The difference between parent and TP \code{Log P} values. Ignored if \emph{both}
#'   \code{parent_LogP} and \code{TP_LogP} are specified. (\strong{optional})
#'
#'   }
#'
#'   Other columns are allowed, and will be included in the final object. Multiple TPs for a single parent are specified
#'   by repeating the value within \code{parent_} columns.
