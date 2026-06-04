# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

getFeaturesOpenMSParamDefs <- memoise(\() list(
    noiseThrInt = list(
        default = 1000,
        description = "Intensity threshold",
        type = "number",
        positive = TRUE,
        finite = TRUE
    ),
    chromSNR = list(
        default = 3,
        description = "Chromatographic signal-to-noise ratio threshold",
        type = "number",
        positive = TRUE,
        finite = TRUE
    ),
    chromFWHM = list(
        default = 5,
        description = "Chromatographic full width at half maximum",
        type = "number",
        positive = TRUE,
        finite = TRUE
    ),
    mzPPM = list(
        default = defaultLim("mz", "medium_rel"),
        description = "Mass-to-charge ratio tolerance in parts per million",
        type = "number",
        positive = TRUE,
        finite = TRUE
    ),
    reEstimateMTSD = list(
        default = TRUE,
        description = "Re-estimate mass trace signal-to-noise ratio",
        type = "flag"
    ),
    traceTermCriterion = list(
        default = "sample_rate",
        description = "Criterion for trace termination",
        type = "choice",
        choices = c("sample_rate", "outliers")
    ),
    traceTermOutliers = list(
        default = 5,
        description = "Number of outliers for trace termination",
        type = "number",
        positive = TRUE,
        finite = TRUE
    ),
    minSampleRate = list(
        default = 0.5,
        description = "Minimum sample rate for feature detection",
        type = "number",
        lower = 0,
        upper = 1,
        finite = TRUE
    ),
    minTraceLength = list(
        default = 3,
        description = "Minimum trace length",
        type = "number",
        positive = TRUE,
        finite = TRUE
    ),
    maxTraceLength = list(
        default = -1,
        description = "Maximum trace length (-1 for unlimited)",
        type = "number",
        lower = -1,
        finite = TRUE
    ),
    widthFiltering = list(
        default = "fixed",
        description = "Width filtering method",
        type = "choice",
        choices = c("fixed", "off", "auto")
    ),
    minFWHM = list(
        default = 1,
        description = "Minimum full width at half maximum",
        type = "number",
        positive = TRUE,
        finite = TRUE
    ),
    maxFWHM = list(
        default = 30,
        description = "Maximum full width at half maximum",
        type = "number",
        positive = TRUE,
        finite = TRUE
    ),
    traceSNRFiltering = list(
        default = FALSE,
        description = "Enable trace signal-to-noise ratio filtering",
        type = "flag"
    ),
    localRTRange = list(
        default = 10,
        description = "Local retention time range for feature detection",
        type = "number",
        positive = TRUE,
        finite = TRUE
    ),
    localMZRange = list(
        default = 6.5,
        description = "Local mass-to-charge range for feature detection",
        type = "number",
        positive = TRUE,
        finite = TRUE
    ),
    isotopeFilteringModel = list(
        default = "metabolites (5% RMS)",
        description = "Isotope filtering model",
        type = "choice",
        choices = c("metabolites (5% RMS)", "metabolites (10% RMS)", "lipids (5% RMS)", "lipids (10% RMS)")
    ),
    MZScoring13C = list(
        default = FALSE,
        description = "Enable 13C scoring for mass-to-charge",
        type = "flag"
    ),
    useSmoothedInts = list(
        default = TRUE,
        description = "Use smoothed intensities",
        type = "flag"
    ),
    extraOpts = list(
        default = NULL,
        description = "Extra options as character list",
        type = "list",
        types = "character",
        names = "unique",
        null.ok = TRUE
    ),
    useFFMIntensities = list(
        default = FALSE,
        description = "Use feature finding MS intensities",
        type = "flag"
    )
))

FeaturesOpenMSParam <- setClass("FeaturesOpenMSParam", contains = "param")
setMethod("initialize", "FeaturesOpenMSParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "FeaturesOpenMSParam", baseName = "FeaturesOpenMSParam",
                   description = "Parameters for OpenMS feature detection", version = "1.0",
                   definitions = getFeaturesOpenMSParamDefs(), ...)
})
