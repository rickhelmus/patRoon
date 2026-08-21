## SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
##
## SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
NULL

#' @export
getSelectIonsParamDefs <- paramConfigDefsFact(list(
    onlyMonoIso = list(
        default = TRUE,
        description = "Keep only monoisotopic features when isotope information is present",
        type = "flag"
    ),
    chargeMismatch = list(
        default = "adduct",
        description = "How to resolve charge mismatches between adduct and isotope annotations",
        type = "choice",
        typeCheckArgs = list(choices = c("isotope", "adduct", "none", "ignore"))
    )
))

#' @export
SelectIonsParam <- setClass("SelectIonsParam", contains = "param")
setMethod("initialize", "SelectIonsParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "SelectIonsParam", baseName = "SelectIonsParam",
                   description = "Parameters for selectIons", version = "1.0",
                   definitions = getSelectIonsParamDefs(), ...)
})
