# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

getFormulasSIRIUSParamDefs <- paramConfigDefsFact(list(
    specSimParams = list(
        advanced = TRUE,
        default = getDefSpecSimParams(removePrecursor = TRUE),
        description = "Spectral similarity parameters used for annotation similarity calculations",
        type = "specSimParams"
    ),
    getFingerprints = list(
        advanced = FALSE,
        default = FALSE,
        description = "Load CSI:FingerID fingerprints",
        type = "flag"
    ),
    login = list(
        advanced = FALSE,
        default = "check",
        description = "SIRIUS login credentials or login mode"
    ),
    alwaysLogin = list(
        advanced = FALSE,
        default = FALSE,
        description = "Always log in to SIRIUS"
    ),
    calculateFeatures = list(
        advanced = FALSE,
        default = FALSE,
        description = "If TRUE, calculate formulas for individual features instead of directly for the feature group",
        type = "flag"
    ),
    featThreshold = list(
        advanced = FALSE,
        default = 0,
        description = "Feature threshold for consensus formula generation (if calculateFeatures==TRUE)",
        type = "number",
        typeCheckArgs = list(lower = 0, upper = 1, finite = TRUE)
    ),
    featThresholdAnn = list(
        advanced = FALSE,
        default = 0.75,
        description = "Feature annotation threshold for consensus formula generation (if calculateFeatures==TRUE)",
        type = "number",
        typeCheckArgs = list(lower = 0, upper = 1, finite = TRUE)
    ),
    absAlignMzDev = list(
        advanced = FALSE,
        default = defaultLim("mz", "narrow"),
        description = "Absolute m/z deviation for feature m/z alignment (if calculateFeatures==TRUE)",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    minIMSSpecSim = list(
        advanced = FALSE,
        default = 0,
        description = "Minimum spectral similarity to copy annotations from the IMS precursor",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    verbose = list(
        advanced = TRUE,
        default = TRUE,
        description = "Show progress output",
        type = "flag"
    ),
    setThreshold = list(
        advanced = TRUE,
        default = 0,
        description = "Minimum abundance for a candidate among all sets",
        type = "number",
        typeCheckArgs = list(lower = 0, upper = 1, finite = TRUE)
    ),
    setThresholdAnn = list(
        advanced = TRUE,
        default = 0,
        description = "As setThreshold, but only taking sets into account in which features of the feature group are present",
        type = "number",
        typeCheckArgs = list(lower = 0, upper = 1, finite = TRUE)
    ),
    setAvgSpecificScores = list(
        advanced = TRUE,
        default = FALSE,
        description = "If TRUE, average set specific scores (e.g. MS/MS match)",
        type = "flag"
    )
))

#' @export
FormulasSIRIUSParam <- setClass("FormulasSIRIUSParam", contains = "param")

setMethod("initialize", "FormulasSIRIUSParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "FormulasSIRIUSParam", baseName = "sirius_formulas", subDir = "annotation",
                   description = "Parameters for SIRIUS formula generation", version = "1.0",
                   definitions = getFormulasSIRIUSParamDefs(), ...)
})

setValidity("FormulasSIRIUSParam", function(object)
{
    parsFilled <- paramListFillDefaults(object@data, object@definitions)
    ac <- checkmate::makeAssertCollection()
    assertSIRIUSLogin(parsFilled$login, parsFilled$alwaysLogin, add = ac)
    OK <- tryCatch(checkmate::reportAssertions(ac), error = function(e) e)
    return(OK)
})
