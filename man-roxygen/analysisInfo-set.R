#' @describeIn <%=class%> Modifies the \link[=analysis-information]{analysis information} of this \code{features} object.
#'   This is primarily intended to change or add analysis metadata columns or can be used to re-order analysis. The
#'   removal or addition of analyses and changes to the \code{"analysis"} column are not supported. This function
#'   performs several internal updates after analysis information modifications. Hence, never attempt to change the
#'   \code{analysisInfo} slot directly.
