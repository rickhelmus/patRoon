# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
NULL

#' @export
getEstimateIDConfidenceParamDefs <- paramConfigDefsFact(list(
    absMzDev = list(
        default = defaultLim("mz", "medium"),
        description = "Maximum absolute m/z deviation",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    specSimParams = list(
        default = getDefSpecSimParams(removePrecursor = TRUE),
        description = "Spectrum similarity parameters",
        type = "specSimParams"
    ),
    checkFragments = list(
        default = c("mz", "formula", "compound"),
        description = "Types of MS/MS fragments for suspect fragment matches",
        type = "subset",
        typeCheckArgs = list(choices = c("mz", "formula", "compound"))
    ),
    formulasNormalizeScores = list(
        default = "max",
        description = "Formula score normalization method",
        type = "normalizationMethod",
        typeCheckArgs = list(withNone = FALSE)
    ),
    compoundsNormalizeScores = list(
        default = "max",
        description = "Compound score normalization method",
        type = "normalizationMethod",
        typeCheckArgs = list(withNone = FALSE)
    )
))

#' @export
EstimateIDConfidenceParam <- setClass("EstimateIDConfidenceParam", contains = "param")
setMethod("initialize", "EstimateIDConfidenceParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "EstimateIDConfidenceParam", baseName = "EstimateIDConfidenceParam",
                   description = "Parameters for estimateIDConfidence", version = "1.0",
                   definitions = getEstimateIDConfidenceParamDefs(), ...)
})