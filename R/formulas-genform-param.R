# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

getFormulasGenFormParamDefs <- paramConfigDefsFact(list(
    specSimParams = list(
        default = getDefSpecSimParams(removePrecursor = TRUE),
        description = "Spectral similarity parameters used for annotation similarity calculations",
        type = "specSimParams"
    ),
    relMzDev = list(
        advanced = TRUE,
        default = defaultLim("mz", "narrow_rel"),
        description = "Maximum relative deviation between measured and candidate formula m/z (in ppm)",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    elements = list(
        default = "CHNOP",
        description = "Elements to be considered for formulae calculation",
        type = "string"
    ),
    hetero = list(
        advanced = TRUE,
        default = TRUE,
        description = "Only consider formulae with at least one hetero atom",
        type = "flag"
    ),
    oc = list(
        advanced = TRUE,
        default = FALSE,
        description = "Only consider organic formulae (at least one carbon atom)",
        type = "flag"
    ),
    thrMS = list(
        advanced = TRUE,
        default = NULL,
        description = "Threshold for GenForm MS score (isoScore)",
        type = "number",
        typeCheckArgs = list(finite = TRUE),
        null.ok = TRUE
    ),
    thrMSMS = list(
        advanced = TRUE,
        default = NULL,
        description = "Threshold for GenForm MS/MS score",
        type = "number",
        typeCheckArgs = list(finite = TRUE),
        null.ok = TRUE
    ),
    thrComb = list(
        advanced = TRUE,
        default = NULL,
        description = "Threshold for GenForm combined score (combMatch)",
        type = "number",
        typeCheckArgs = list(finite = TRUE),
        null.ok = TRUE
    ),
    maxCandidates = list(
        advanced = TRUE,
        default = Inf,
        description = "Maximum number of candidates until calculations are aborted. See ?generateFormulasGenform for details!"
    ),
    extraOpts = list(
        advanced = TRUE,
        default = NULL,
        description = "Extra GenForm command line options",
        type = "character",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    calculateFeatures = list(
        advanced = TRUE,
        default = FALSE,
        description = "If TRUE, calculate formulas for individual features instead of directly for the feature group",
        type = "flag"
    ),
    featThreshold = list(
        advanced = TRUE,
        default = 0,
        description = "Feature threshold for consensus formula generation (if calculateFeatures==TRUE)",
        type = "number",
        typeCheckArgs = list(lower = 0, upper = 1, finite = TRUE)
    ),
    featThresholdAnn = list(
        advanced = TRUE,
        default = 0.75,
        description = "Feature annotation threshold for consensus formula generation (if calculateFeatures==TRUE)",
        type = "number",
        typeCheckArgs = list(lower = 0, upper = 1, finite = TRUE)
    ),
    absAlignMzDev = list(
        advanced = TRUE,
        default = defaultLim("mz", "narrow"),
        description = "Absolute m/z deviation for feature m/z alignment (if calculateFeatures==TRUE)",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    MSMode = list(
        advanced = TRUE,
        default = "both",
        description = "Do calculations with data from MS only, MS/MS only, or both",
        type = "choice",
        typeCheckArgs = list(choices = c("ms", "msms", "both"))
    ),
    isolatePrec = list(
        advanced = TRUE,
        default = TRUE,
        description = "Settings for isolation of the precursor peak with its isotopes (TRUE for defaults, FALSE to disable or a list with settings)"
    ),
    minIMSSpecSim = list(
        advanced = TRUE,
        default = 0,
        description = "Minimum spectral similarity to copy annotations from the IMS precursor",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    timeout = list(
        advanced = TRUE,
        default = 120,
        description = "Maximum execution time per GenForm command (seconds)",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    topMost = list(
        advanced = TRUE,
        default = 50,
        description = "Maximum number of top ranked candidates to keep",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    batchSize = list(
        advanced = TRUE,
        default = 8,
        description = "Maximum number of GenForm commands per parallel batch",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
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
FormulasGenFormParam <- setClass("FormulasGenFormParam", contains = "param")

setMethod("initialize", "FormulasGenFormParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "FormulasGenFormParam", baseName = "genform", subDir = "annotation",
                   description = "Parameters for GenForm formula generation", version = "1.0",
                   definitions = getFormulasGenFormParamDefs(), ...)
})

setValidity("FormulasGenFormParam", function(object)
{
    parsFilled <- paramListFillDefaults(object@data, object@definitions)
    if (is.infinite(parsFilled$maxCandidates))
        parsFilled$maxCandidates <- 0L
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertCount(parsFilled$maxCandidates, .var.name = "maxCandidates", add = ac)
    if (!is.logical(parsFilled$isolatePrec))
        assertPListIsolatePrecParams(parsFilled$isolatePrec, .var.name = "isolatePrec", add = ac)
    OK <- tryCatch(checkmate::reportAssertions(ac), error = function(e) e)
    
    return(OK)
})
