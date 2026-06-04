# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
NULL

getFeaturesPiekParamDefs <- paramConfigDefsFact(list(
    genEICParams = list(
        default = getPiekEICParams(),
        description = "Parameters for EIC generation (use getPiekEICParams())",
        type = "list",
        null.ok = FALSE
    ),
    peakParams = list(
        default = getDefPeakParams("chrom", "piek"),
        description = "Parameters for peak detection (use getDefPeakParams())",
        type = "list",
        null.ok = FALSE
    ),
    IMS = list(
        default = FALSE,
        description = "Use ion mobility separation (IMS)",
        type = "flag"
    ),
    assignMethod = list(
        default = "basepeak",
        description = "Method to assign m/z/mobility across EIC datapoints",
        type = "choice",
        choices = c("basepeak", "weighted.mean")
    ),
    assignRTWindow = list(
        default = defaultLim("retention", "very_narrow"),
        description = "Retention time window (+/- s) for assignment",
        type = "number",
        lower = 0,
        finite = TRUE
    ),
    rtWindowDup = list(
        default = defaultLim("retention", "narrow"),
        description = "RT window for duplicate feature detection",
        type = "number",
        lower = 0,
        finite = TRUE
    ),
    mzWindowDup = list(
        default = defaultLim("mz", "medium"),
        description = "m/z window for duplicate feature detection",
        type = "number",
        lower = 0,
        finite = TRUE
    ),
    mobWindowDup = list(
        default = defaultLim("mobility", "medium"),
        description = "Mobility window for duplicate feature detection",
        type = "number",
        lower = 0,
        finite = TRUE
    ),
    minPeakOverlapDup = list(
        default = 0.25,
        description = "Minimum retention time overlap (fraction) to consider duplicates",
        type = "number",
        lower = 0,
        upper = 1,
        finite = TRUE
    ),
    minIntensityIMS = list(
        default = 25,
        description = "Minimum intensity for IMS datapoints",
        type = "number",
        positive = TRUE,
        finite = TRUE
    ),
    EICBatchSize = list(
        default = Inf,
        description = "Number of EICs processed per batch",
        type = "number",
        positive = TRUE,
        finite = FALSE
    ),
    keepDups = list(
        default = FALSE,
        description = "Keep duplicate / non-centered features",
        type = "flag"
    ),
    verbose = list(
        default = TRUE,
        description = "Verbose output",
        type = "flag"
    )
))

FeaturesPiekParam <- setClass("FeaturesPiekParam", contains = "param")
setMethod("initialize", "FeaturesPiekParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "FeaturesPiekParam", baseName = "FeaturesPiekParam",
                   description = "Parameters for piek feature detection", version = "1.0",
                   definitions = getFeaturesPiekParamDefs(), ...)
})


setValidity("FeaturesPiekParam", function(object)
{
    parsFilled <- paramListFillDefaults(object@data, object@definitions)
    
    ac <- checkmate::makeAssertCollection()
    assertPiekGenEICParams(parsFilled$genEICParams, .var.name = "genEICParams", add = ac)
    assertFindPeakParams(parsFilled$peakParams, .var.name = "peakParams", add = ac)
    
    OK <- tryCatch(checkmate::reportAssertions(ac), error = function(e) e)
    return(OK)
})
