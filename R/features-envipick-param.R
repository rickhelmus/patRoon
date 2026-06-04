# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
NULL

getFeaturesEnviPickParamDefs <- paramConfigDefsFact(list(
    extraOpts = list(
        default = NULL,
        description = "List with extra options passed to enviPick::enviPickwrap()",
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

FeaturesEnviPickParam <- setClass("FeaturesEnviPickParam", contains = "param")
setMethod("initialize", "FeaturesEnviPickParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "FeaturesEnviPickParam", baseName = "FeaturesEnviPickParam",
                   description = "Parameters for enviPick feature detection", version = "1.0",
                   definitions = getFeaturesEnviPickParamDefs(), ...)
})
