# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

getFeatureGroupsOpenMSParamDefs <- paramConfigDefsFact(list(
    rtalign = list(
        default = TRUE,
        description = "If TRUE then retention times are aligned",
        type = "flag"
    ),
    QT = list(
        default = FALSE,
        description = "If TRUE then FeatureLinkerUnlabeledQT is used instead of FeatureLinkerUnlabeled",
        type = "flag"
    ),
    maxAlignRT = list(
        default = defaultLim("retention", "wide"),
        description = "Maximum retention time difference (seconds) for feature pairing during alignment",
        type = "number",
        finite = TRUE
    ),
    maxAlignMZ = list(
        default = defaultLim("mz", "medium"),
        description = "Maximum m/z difference for feature pairing during alignment",
        type = "number",
        finite = TRUE
    ),
    maxGroupRT = list(
        default = defaultLim("retention", "medium"),
        description = "Maximum retention time difference (seconds) for grouping of features",
        type = "number",
        finite = TRUE
    ),
    maxGroupMZ = list(
        default = defaultLim("mz", "medium"),
        description = "Maximum m/z difference for grouping of features",
        type = "number",
        finite = TRUE
    ),
    extraOptsRT = list(
        default = NULL,
        description = "Extra command line options passed to MapAlignerPoseClustering as named list",
        type = "list",
        names = "unique",
        null.ok = TRUE
    ),
    extraOptsGroup = list(
        default = NULL,
        description = "Extra command line options passed to FeatureLinkerUnlabeled(QT) as named list",
        type = "list",
        names = "unique",
        null.ok = TRUE
    ),
    verbose = list(
        default = TRUE,
        description = "Verbose output",
        type = "flag"
    )
))

FeatureGroupsOpenMSParam <- setClass("FeatureGroupsOpenMSParam", contains = "param")
setMethod("initialize", "FeatureGroupsOpenMSParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "FeatureGroupsOpenMSParam", baseName = "FeatureGroupsOpenMSParam",
                   description = "Parameters for OpenMS feature grouping", version = "1.0",
                   definitions = getFeatureGroupsOpenMSParamDefs(), ...)
})
