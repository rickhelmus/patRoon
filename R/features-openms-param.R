# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

getFeaturesOpenMSParamDefs <- paramConfigDefsFact(list(
    noiseThrInt = list(
        advanced = FALSE,
        default = 1000,
        description = "Intensity threshold",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    chromSNR = list(
        advanced = FALSE,
        default = 3,
        description = "Chromatographic signal-to-noise ratio threshold",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    chromFWHM = list(
        advanced = FALSE,
        default = 5,
        description = "Chromatographic full width at half maximum",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    mzPPM = list(
        advanced = FALSE,
        default = defaultLim("mz", "medium_rel"),
        description = "Mass-to-charge ratio tolerance in parts per million",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    reEstimateMTSD = list(
        advanced = FALSE,
        default = TRUE,
        description = "Re-estimate mass trace signal-to-noise ratio",
        type = "flag"
    ),
    traceTermCriterion = list(
        advanced = FALSE,
        default = "sample_rate",
        description = "Criterion for trace termination",
        type = "choice",
        typeCheckArgs = list(choices = c("sample_rate", "outliers"))
    ),
    traceTermOutliers = list(
        advanced = FALSE,
        default = 5,
        description = "Number of outliers for trace termination",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    minSampleRate = list(
        advanced = FALSE,
        default = 0.5,
        description = "Minimum sample rate for feature detection",
        type = "number",
        typeCheckArgs = list(lower = 0, upper = 1, finite = TRUE)
    ),
    minTraceLength = list(
        advanced = FALSE,
        default = 3,
        description = "Minimum trace length",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    maxTraceLength = list(
        advanced = FALSE,
        default = -1,
        description = "Maximum trace length (-1 for unlimited)",
        type = "number",
        typeCheckArgs = list(lower = -1, finite = TRUE)
    ),
    widthFiltering = list(
        advanced = FALSE,
        default = "fixed",
        description = "Width filtering method",
        type = "choice",
        typeCheckArgs = list(choices = c("fixed", "off", "auto"))
    ),
    minFWHM = list(
        advanced = FALSE,
        default = 1,
        description = "Minimum full width at half maximum",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    maxFWHM = list(
        advanced = FALSE,
        default = 30,
        description = "Maximum full width at half maximum",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    traceSNRFiltering = list(
        advanced = FALSE,
        default = FALSE,
        description = "Enable trace signal-to-noise ratio filtering",
        type = "flag"
    ),
    localRTRange = list(
        advanced = TRUE,
        default = 10,
        description = "Local retention time range for feature detection",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    localMZRange = list(
        advanced = TRUE,
        default = 6.5,
        description = "Local mass-to-charge range for feature detection",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    isotopeFilteringModel = list(
        advanced = FALSE,
        default = "metabolites (5% RMS)",
        description = "Isotope filtering model",
        type = "choice",
        typeCheckArgs = list(choices = c("metabolites (5% RMS)", "metabolites (10% RMS)", "lipids (5% RMS)", "lipids (10% RMS)"))
    ),
    MZScoring13C = list(
        advanced = TRUE,
        default = FALSE,
        description = "Enable 13C scoring for mass-to-charge",
        type = "flag"
    ),
    useSmoothedInts = list(
        advanced = TRUE,
        default = TRUE,
        description = "Use smoothed intensities",
        type = "flag"
    ),
    extraOpts = list(
        advanced = TRUE,
        default = NULL,
        description = "Extra options as character list",
        type = "list",
        typeCheckArgs = list(types = "character", names = "unique", null.ok = TRUE)
    ),
    useFFMIntensities = list(
        advanced = TRUE,
        default = FALSE,
        description = "Use feature finding MS intensities",
        type = "flag"
    ),
    verbose = list(
        advanced = TRUE,
        default = TRUE,
        description = "Verbose output",
        type = "flag"
    )
))

FeaturesOpenMSParam <- setClass("FeaturesOpenMSParam", contains = "param")
setMethod("initialize", "FeaturesOpenMSParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "FeaturesOpenMSParam", baseName = "openms", subDir = "features",
                   description = "Parameters for OpenMS feature detection", version = "1.0",
                   definitions = getFeaturesOpenMSParamDefs(), ...)
})
