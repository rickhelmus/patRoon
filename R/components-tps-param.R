# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

#' @export
getComponentsTPsParamDefs <- paramConfigDefsFact(list(
    ignoreParents = list(
        default = FALSE,
        description = "If TRUE, feature groups present in both parents and TPs are ignored during linkage.",
        type = "flag"
    ),
    minRTDiff = list(
        default = 20,
        description = "Minimum retention time difference (seconds) to consider when linking parents and TPs.",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    specSimParams = list(
        default = getDefSpecSimParams(),
        description = "MS/MS spectral similarity parameters.",
        type = "specSimParams"
    ),
    IMS = list(
        default = "maybe",
        description = "Which feature groups are considered for componentization in IMS workflows ('maybe', 'both', FALSE, TRUE)",
        type = "IMS"
    )
))

#' @export
ComponentsTPsParam <- setClass("ComponentsTPsParam", contains = "param")

setMethod("initialize", "ComponentsTPsParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "ComponentsTPsParam", baseName = "ComponentsTPsParam",
                   description = "Parameters for TP component generation", version = "1.0",
                   definitions = getComponentsTPsParamDefs(), ...)
})
