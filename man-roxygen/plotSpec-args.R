#' @param groupName The name of the feature group for which a plot should be made. To compare spectra, two group names
#'   can be specified.
#' @param title The title of the plot. If \code{NULL} a title will be automatically made.
#' @param normalized Controls intensity normalization. Should be \code{FALSE} (don't normalize), \code{TRUE} (normalize)
#'   or \code{"multiple"} (only normalizes if multiple spectra are plotted).
#' @param showLegend Set to \code{TRUE} to show a legend.
#'
#' <%=if (exists("withSpecAna")) "@param analysis The name of the analysis for which a plot should be made. If
#' \\code{NULL} then data from the feature group averaged peak list is used. When comparing spectra, either
#' \\code{NULL} or the analyses for both spectra should be specified." %>
#'
#' <%=if (exists("featAnnArgs")) "@param mincex The formula annotation labels are automatically scaled. The
#' \\code{mincex} argument forces a minimum \\code{cex} value for readability." %>

