# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

getFeatureGroupsSIRIUSParamDefs <- paramConfigDefsFact(list(
    noiseIntensity = list(
        default = NULL,
        description = "Noise intensity threshold",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE, null.ok = TRUE)
    ),
    alignMaxRTDev = list(
        default = NULL,
        description = "Maximum retention time deviation for alignment",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE, null.ok = TRUE)
    ),
    minSNR = list(
        default = NULL,
        description = "Minimum signal-to-noise ratio",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE, null.ok = TRUE)
    ),
    login = list(
        default = "check",
        description = "SIRIUS login credentials or login mode"
    ),
    alwaysLogin = list(
        default = FALSE,
        description = "Always log in to SIRIUS"
    ),
    verbose = list(
        default = TRUE,
        description = "Verbose output",
        type = "flag"
    )
))

FeatureGroupsSIRIUSParam <- setClass("FeatureGroupsSIRIUSParam", contains = "param")
setMethod("initialize", "FeatureGroupsSIRIUSParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "FeatureGroupsSIRIUSParam", baseName = "FeatureGroupsSIRIUSParam",
                   description = "Parameters for SIRIUS feature grouping", version = "1.0",
                   definitions = getFeatureGroupsSIRIUSParamDefs(), ...)
})

setValidity("FeatureGroupsSIRIUSParam", function(object)
{
    parsFilled <- paramListFillDefaults(object@data, object@definitions)
    ac <- checkmate::makeAssertCollection()
    assertSIRIUSLogin(parsFilled$login, parsFilled$alwaysLogin, add = ac)
    OK <- tryCatch(checkmate::reportAssertions(ac), error = function(e) e)
    return(OK)
})
