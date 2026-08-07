# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
NULL

#' @export
getScreenSuspectsParamDefs <- paramConfigDefsFact(list(
    rtWindow = list(
        default = defaultLim("retention", "medium"),
        description = "Retention time window (in seconds) for matching suspects",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    mzWindow = list(
        default = defaultLim("mz", "medium"),
        description = "m/z window for matching suspects",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    IMSMatchParams = list(
        default = NULL,
        description = "IMS matching parameters",
        type = "IMSMatchParams",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    skipInvalid = list(
        default = TRUE,
        description = "If TRUE, ignore suspects with invalid data",
        type = "flag"
    ),
    prefCalcChemProps = list(
        default = TRUE,
        description = "If TRUE, prefer to calculate chemical properties",
        type = "flag"
    ),
    neutralChemProps = list(
        default = FALSE,
        description = "If TRUE, use neutral chemical properties",
        type = "flag"
    ),
    onlyHits = list(
        default = FALSE,
        description = "If TRUE, remove feature groups without suspect matches",
        type = "flag"
    ),
    amend = list(
        default = FALSE,
        description = "If TRUE, amend existing suspect matches",
        type = "flag"
    )
))

#' @export
ScreenSuspectsParam <- setClass("ScreenSuspectsParam", contains = "param")
setMethod("initialize", "ScreenSuspectsParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "ScreenSuspectsParam", baseName = "ScreenSuspectsParam",
                   description = "Parameters for screenSuspects", version = "1.0",
                   definitions = getScreenSuspectsParamDefs(), ...)
})
