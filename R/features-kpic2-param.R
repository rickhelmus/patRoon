# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
NULL

getFeaturesKPIC2ParamDefs <- paramConfigDefsFact(list(
    kmeans = list(
        default = TRUE,
        description = "If TRUE, use kmeans clustering for feature detection",
        type = "flag"
    ),
    level = list(
        default = 1000,
        description = "Level parameter passed to KPIC::getPIC or KPIC::getPIC.kmeans",
        type = "number",
        positive = TRUE
    ),
    extraOpts = list(
        default = NULL,
        description = "List with extra options passed to KPIC::getPIC or KPIC::getPIC.kmeans",
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

FeaturesKPIC2Param <- setClass("FeaturesKPIC2Param", contains = "param")
setMethod("initialize", "FeaturesKPIC2Param", function(.Object, ...)
{
    callNextMethod(.Object, name = "FeaturesKPIC2Param", baseName = "FeaturesKPIC2Param",
                   description = "Parameters for KPIC2 feature detection", version = "1.0",
                   definitions = getFeaturesKPIC2ParamDefs(), ...)
})
