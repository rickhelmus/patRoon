# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

getFeatureGroupsKPIC2ParamDefs <- paramConfigDefsFact(list(
    rtalign = list(
        default = TRUE,
        description = "If TRUE then retention times are aligned",
        type = "flag"
    ),
    loadRawData = list(
        default = TRUE,
        description = "If TRUE then raw data is loaded into memory",
        type = "flag"
    ),
    groupArgs = list(
        default = list(tolerance = c(0.005, 12)),
        description = "Named list with extra parameters passed to KPIC::PICset.group",
        type = "list",
        names = "unique"
    ),
    alignArgs = list(
        default = list(),
        description = "Named list with extra parameters passed to KPIC::PICset.align",
        type = "list",
        names = "unique"
    ),
    verbose = list(
        default = TRUE,
        description = "Verbose output",
        type = "flag"
    )
))

FeatureGroupsKPIC2Param <- setClass("FeatureGroupsKPIC2Param", contains = "param")
setMethod("initialize", "FeatureGroupsKPIC2Param", function(.Object, ...)
{
    callNextMethod(.Object, name = "FeatureGroupsKPIC2Param", baseName = "FeatureGroupsKPIC2Param",
                   description = "Parameters for KPIC2 feature grouping", version = "1.0",
                   definitions = getFeatureGroupsKPIC2ParamDefs(), ...)
})