#' @param IMS Specifies which feature groups are considered <%=consider%> in IMS workflows. The following options are
#'   valid: \itemize{
#'
#'   \item \code{"both"}: Selects IMS and non-IMS features.
#'
#'   \item \code{"maybe"}: Selects non-IMS features and IMS features without assigned IMS parent.
#'
#'   \item \code{FALSE}: Selects only non-IMS features.
#'
#'   \item \code{TRUE}: Selects only IMS features.
#'
#'   }
#'
#'   <%= if (is.character(append)) append else "" %>
