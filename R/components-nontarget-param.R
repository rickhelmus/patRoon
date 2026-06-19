# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

#' @export
getComponentsNontargetParamDefs <- paramConfigDefsFact(list(
    rtRange = list(
        default = c(-120, 120),
        description = "A numeric vector containing the minimum and maximum retention time (in seconds) between homologues.",
        type = "numeric",
        typeCheckArgs = list(len = 2, any.missing = FALSE, finite = TRUE)
    ),
    mzRange = list(
        default = c(5, 120),
        description = "A numeric vector specifying the minimum and maximum m/z increment of a homologous series.",
        type = "numeric",
        typeCheckArgs = list(len = 2, lower = 0, any.missing = FALSE, finite = TRUE)
    ),
    elements = list(
        default = c("C", "H", "O"),
        description = "A character vector with elements to be considered for detection of repeating units.",
        type = "character",
        typeCheckArgs = list(min.chars = 1, min.len = 1, any.missing = FALSE)
    ),
    rtDev = list(
        default = defaultLim("retention", "wide"),
        description = "Retention time tolerance.",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    absMzDev = list(
        default = defaultLim("mz", "narrow"),
        description = "Absolute m/z tolerance.",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    absMzDevLink = list(
        default = defaultLim("mz", "medium"),
        description = "Absolute m/z tolerance when linking series.",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    traceHack = list(
        default = all(R.Version()[c("major", "minor")] >= c(3, 4)),
        description = "Enable workaround for nontarget::homol.search incompatibility with R > 3.3.3.",
        type = "flag"
    )
))

#' @export
ComponentsNontargetParam <- setClass("ComponentsNontargetParam", contains = "param")

setMethod("initialize", "ComponentsNontargetParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "ComponentsNontargetParam", baseName = "ComponentsNontargetParam",
                   description = "Parameters for nontarget component generation", version = "1.0",
                   definitions = getComponentsNontargetParamDefs(), ...)
})
