#' <%= if (algo != "bruker") "@param calculateFeatures If \\code{TRUE} fomulae are first calculated for all features
#' prior to feature group assignment (see \\verb{Candidate assignment} in \\code{\\link{generateFormulas}})." %>
#' @param featThreshold If \code{calculateFeatures=TRUE}: minimum presence (\samp{0-1}) of a formula in all features
#'   before it is considered as a candidate for a feature group. For instance, \code{featThreshold=0.75} dictates that a
#'   formula should be present in at least 75\% of the features inside a feature group.
#' @param featThresholdAnn As \code{featThreshold}, but only considers features with annotations. For instance,
#'   \code{featThresholdAnn=0.75} dictates that a formula should be present in at least 75\% of the features with
#'   annotations inside a feature group. <%= if (algo != "bruker") "@param topMost Only keep this number of candidates
#'   (per feature group) with highest score." %> <%= if (algo == "sirius") "Sets the \\option{--candidates} command line
#'   option." %>
#' @param absAlignMzDev When the group formula annotation consensus is made from feature annotations, the \emph{m/z}
#'   values of annotated MS/MS fragments may slightly deviate from those of the corresponding group MS/MS peak list. The
#'   \code{absAlignMzDev} argument specifies the maximum \emph{m/z} window used to re-align the mass peaks.
#' <%= if (algo != "bruker") "@param topMost Only keep this number of candidates (per feature group) with highest
#'   score." %> <%= if (algo == "sirius") "Sets the \\option{--candidates} command line option." %>
#' @param \dots \setsWF Further arguments passed to the non-sets workflow method.
