#' @describeIn <%=what%> Generates a consensus of results from multiple
#'   objects. In order to rank the consensus candidates, first
#'   each of the candidates are scored based on their original ranking
#'   (the scores are normalized and the highest ranked candidate gets value
#'   \samp{1}). The (weighted) mean is then calculated for all scorings of each
#'   candidate to derive the final ranking (if an object lacks the candidate its
#'   score will be \samp{0}). The original rankings for each object is stored in
#'   the \code{rank} columns.
#'
#' @param rankWeights A numeric vector with weights of to calculate the mean
#'   ranking score for each candidate. The value will be re-cycled if necessary,
#'   hence, the default value of \samp{1} means equal weights for all considered
#'   objects.
