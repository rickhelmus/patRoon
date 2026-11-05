# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

getCompoundsSIRIUSParamDefs <- paramConfigDefsFact(list(
    specSimParams = list(
        default = getDefSpecSimParams(removePrecursor = TRUE),
        description = "Spectral similarity parameters used for annotation similarity calculations",
        type = "specSimParams"
    ),
    topMost = list(
        default = 100,
        description = "Maximum number of top-ranked structures to keep per feature group",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    login = list(
        default = "check",
        description = "SIRIUS login credentials or login mode"
    ),
    alwaysLogin = list(
        default = FALSE,
        description = "Always log in to SIRIUS"
    ),
    minIMSSpecSim = list(
        default = 0,
        description = "Minimum spectral similarity to copy annotations from the IMS precursor",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    verbose = list(
        default = TRUE,
        description = "Show progress output",
        type = "flag"
    ),
    setThreshold = list(
        default = 0,
        description = "Minimum abundance for a candidate among all sets",
        type = "number",
        typeCheckArgs = list(lower = 0, upper = 1, finite = TRUE)
    ),
    setThresholdAnn = list(
        default = 0,
        description = "As setThreshold, but only taking sets into account in which features of the feature group are present",
        type = "number",
        typeCheckArgs = list(lower = 0, upper = 1, finite = TRUE)
    ),
    setAvgSpecificScores = list(
        default = FALSE,
        description = "If TRUE, average set specific scores (e.g. MS/MS match)",
        type = "flag"
    )
))

#' @export
CompoundsSIRIUSParam <- setClass("CompoundsSIRIUSParam", contains = "param")

setMethod("initialize", "CompoundsSIRIUSParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "CompoundsSIRIUSParam", baseName = "sirius_compounds", subDir = "annotation",
                   description = "Parameters for SIRIUS compound generation", version = "1.0",
                   definitions = getCompoundsSIRIUSParamDefs(), ...)
})

setValidity("CompoundsSIRIUSParam", function(object)
{
    parsFilled <- paramListFillDefaults(object@data, object@definitions)
    ac <- checkmate::makeAssertCollection()
    assertSIRIUSLogin(parsFilled$login, parsFilled$alwaysLogin, add = ac)
    OK <- tryCatch(checkmate::reportAssertions(ac), error = function(e) e)
    return(OK)
})
