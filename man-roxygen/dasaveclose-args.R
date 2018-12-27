#' @param close,save If \code{TRUE} then Bruker files are closed and saved after
#'   processing with DataAnalysis, respectively. Setting \code{close=TRUE}
#'   prevents that many analyses might be opened simultaneously in DataAnalysis,
#'   which otherwise may use excessive memory or become slow. By default
#'   \code{save} is \code{TRUE} when \code{close} is \code{TRUE}, which is
#'   likely what you want as otherwise any processed data is lost.
