# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

getFeatureGroupsGreedyParamDefs <- paramConfigDefsFact(list(
    rtWindow = list(
        default = defaultLim("retention", "medium"),
        description = "Retention time tolerance window (seconds) for grouping",
        type = "number",
        positive = TRUE,
        finite = TRUE
    ),
    mzWindow = list(
        default = defaultLim("mz", "medium"),
        description = "m/z tolerance window for grouping",
        type = "number",
        positive = TRUE,
        finite = TRUE
    ),
    mobWindow = list(
        default = defaultLim("mobility", "medium"),
        description = "Mobility tolerance window for grouping",
        type = "number",
        positive = TRUE,
        finite = TRUE
    ),
    scoreWeights = list(
        default = c(retention = 1, mz = 1, mobility = 1, intensity = 1),
        description = "Named numeric vector specifying scoring weights for retention, mz, mobility, and intensity",
        type = "number"
    ),
    verbose = list(
        default = TRUE,
        description = "Verbose output",
        type = "flag"
    )
))

FeatureGroupsGreedyParam <- setClass("FeatureGroupsGreedyParam", contains = "param")
setMethod("initialize", "FeatureGroupsGreedyParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "FeatureGroupsGreedyParam", baseName = "FeatureGroupsGreedyParam",
                   description = "Parameters for greedy feature grouping", version = "1.0",
                   definitions = getFeatureGroupsGreedyParamDefs(), ...)
})