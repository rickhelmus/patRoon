# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

#' @export
getComponentsCAMERAParamDefs <- paramConfigDefsFact(list(
    onlyIsotopes = list(
        default = FALSE,
        description = "If TRUE only isotopes are considered when generating components (faster).",
        type = "flag"
    ),
    minSize = list(
        default = 2,
        description = "Minimum number of feature groups to form a component",
        type = "count",
        positive = TRUE
    ),
    relMinReplicates = list(
        default = 0.5,
        description = "The features of a feature group in a component must be present in at least this relative fraction of replicates",
        type = "number",
        lower = 0,
        finite = TRUE
    ),
    extraOpts = list(
        default = NULL,
        description = "Extra options passed directly to CAMERA::annotate() as a named character list",
        type = "list",
        types = "character",
        names = "unique",
        null.ok = TRUE
    ),
    verbose = list(
        default = TRUE,
        description = "Verbose output",
        type = "flag"
    )
))

#' @export
ComponentsCAMERAParam <- setClass("ComponentsCAMERAParam", contains = "param")

setMethod("initialize", "ComponentsCAMERAParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "ComponentsCAMERAParam", baseName = "ComponentsCAMERAParam",
                   description = "Parameters for CAMERA component generation", version = "1.0",
                   definitions = getComponentsCAMERAParamDefs(), ...)
})

