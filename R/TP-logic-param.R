# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

getTPsLogicParamDefs <- paramConfigDefsFact(list(
    minMass = list(
        default = 40,
        description = "The minimum mass for a TP to be kept",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    )
))

#' @export
TPsLogicParam <- setClass("TPsLogicParam", contains = "param")

setMethod("initialize", "TPsLogicParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "TPsLogicParam", baseName = "Logic",
                   description = "Parameters for metabolic logic TP generation", version = "1.0",
                   definitions = getTPsLogicParamDefs(), ...)
})
