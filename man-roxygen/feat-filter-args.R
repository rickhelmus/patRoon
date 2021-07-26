#' @param absMinIntensity,relMinIntensity Minimum absolute/relative intensity for features to be kept. The relative
#'   intensity is determined from the feature with highest intensity (<%=if (feat) "within the same analysis" else "of
#'   all features from all groups" %>). Set to \samp{0} or \code{NULL} to skip this step.
#' @param retentionRange,mzRange,mzDefectRange,chromWidthRange Range of retention time (in seconds), \emph{m/z}, mass
#'   defect (defined as the decimal part of \emph{m/z} values) or chromatographic peak width (in seconds), respectively.
#'   Features outside this range will be removed. Should be a numeric vector with length of two containing the min/max
#'   values. The maximum can be \code{Inf} to specify no maximum range. Set to \code{NULL} to skip this step.
#' @param <%=if (feat) "qualityRange" else "featQualityRange" %> Used to filter features by their peak qualities/scores
#'   (see \code{\link{calculatePeakQualities}}). Should be a named \code{list} with min/max ranges for each
#'   quality/score to be filtered (the \code{\link{featureQualityNames}} function can be used to obtain valid names).
#'   Example: \code{<%=if (feat) "qualityRange" else "featQualityRange" %>=list(ModalityScore=c(0.3, Inf),
#'   SymmetryScore=c(0.5, Inf))}. Set to \code{NULL} to ignore.
#' @param negate If set to \code{TRUE} then filtering operations are performed in opposite manner.
