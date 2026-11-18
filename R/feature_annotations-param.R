# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
NULL

getFilterFeatAnnotationsParamDefs <- paramConfigDefsFact(list(
    minExplainedPeaks = list(
        advanced = FALSE,
        default = NULL,
        description = "Minimum number of explained peaks",
        type = "count",
        typeCheckArgs = list(positive = FALSE, null.ok = TRUE)
    ),
    scoreLimits = list(
        advanced = FALSE,
        default = NULL,
        description = "Minimum and maximum score limits",
        type = "scoreRange",
        typeCheckArgs = list(scNames = NULL)
    ),
    elements = list(
        advanced = FALSE,
        default = NULL,
        description = "Allowed elements in candidate formulae",
        type = "character",
        typeCheckArgs = list(min.chars = 1, min.len = 1, null.ok = TRUE)
    ),
    fragElements = list(
        advanced = FALSE,
        default = NULL,
        description = "Elements required in fragments",
        type = "character",
        typeCheckArgs = list(min.chars = 1, min.len = 1, null.ok = TRUE)
    ),
    lossElements = list(
        advanced = FALSE,
        default = NULL,
        description = "Elements required in neutral losses",
        type = "character",
        typeCheckArgs = list(min.chars = 1, min.len = 1, null.ok = TRUE)
    ),
    fragFormulas = list(
        advanced = FALSE,
        default = NULL,
        description = "Formulae required in fragments",
        type = "character",
        typeCheckArgs = list(min.chars = 1, min.len = 1, null.ok = TRUE)
    ),
    lossFormulas = list(
        advanced = FALSE,
        default = NULL,
        description = "Formulae required in neutral losses",
        type = "character",
        typeCheckArgs = list(min.chars = 1, min.len = 1, null.ok = TRUE)
    ),
    topMost = list(
        advanced = FALSE,
        default = NULL,
        description = "Maximum number of candidates to keep",
        type = "count",
        typeCheckArgs = list(positive = TRUE, null.ok = TRUE)
    ),
    OM = list(
        advanced = TRUE,
        default = FALSE,
        description = "Apply organic matter rules",
        type = "flag"
    ),
    maxLevel = list(
        advanced = FALSE,
        default = NULL,
        description = "Maximum identification level",
        type = "count",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    negate = list(
        advanced = TRUE,
        default = FALSE,
        description = "Apply filters in the opposite manner",
        type = "flag"
    )
))

FilterFeatAnnotationsParam <- setClass("FilterFeatAnnotationsParam", contains = "param")
setMethod("initialize", "FilterFeatAnnotationsParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "FilterFeatAnnotationsParam", baseName = "filter_featannotations",
                   subDir = "annotation", description = "Parameters for feature annotation filtering", version = "1.0",
                   definitions = getFilterFeatAnnotationsParamDefs(), ...)
})


getEstimateIDConfidenceParamDefs <- paramConfigDefsFact(list(
    absMzDev = list(
        advanced = FALSE,
        default = defaultLim("mz", "medium"),
        description = "Maximum absolute m/z deviation",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    specSimParams = list(
        advanced = FALSE,
        default = getDefSpecSimParams(removePrecursor = TRUE),
        description = "Spectrum similarity parameters",
        type = "specSimParams"
    ),
    checkFragments = list(
        advanced = TRUE,
        default = c("mz", "formula", "compound"),
        description = "Types of MS/MS fragments for suspect fragment matches",
        type = "subset",
        typeCheckArgs = list(choices = c("mz", "formula", "compound"))
    ),
    formulasNormalizeScores = list(
        advanced = TRUE,
        default = "max",
        description = "Formula score normalization method",
        type = "normalizationMethod",
        typeCheckArgs = list(withNone = FALSE)
    ),
    compoundsNormalizeScores = list(
        advanced = TRUE,
        default = "max",
        description = "Compound score normalization method",
        type = "normalizationMethod",
        typeCheckArgs = list(withNone = FALSE)
    )
))

#' @export
EstimateIDConfidenceParam <- setClass("EstimateIDConfidenceParam", contains = "param")
setMethod("initialize", "EstimateIDConfidenceParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "EstimateIDConfidenceParam", baseName = "estimate_id_confidence", subDir = "annotation",
                   description = "Parameters for estimateIDConfidence", version = "1.0",
                   definitions = getEstimateIDConfidenceParamDefs(), ...)
})

getAssignMobilitiesCompoundsParamDefs <- paramConfigDefsFact(list(
    IMS = list(
        advanced = TRUE,
        default = TRUE,
        description = "Select IMS features for mobility assignment",
        type = "IMS"
    ),
    from = list(
        advanced = FALSE,
        default = NULL,
        description = "Reference data for mobility and CCS assignment"
    ),
    matchFromBy = list(
        advanced = FALSE,
        default = "InChIKey1",
        description = "Column used to match annotations to reference data",
        type = "choice",
        typeCheckArgs = list(choices = c("InChIKey", "InChIKey1", "InChI", "SMILES", "name"))
    ),
    overwrite = list(
        advanced = FALSE,
        default = FALSE,
        description = "Overwrite existing mobility and CCS values",
        type = "flag"
    ),
    CCSParams = list(
        advanced = FALSE,
        default = NULL,
        description = "Parameters for CCS assignment",
        type = "CCSParams",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    prefCalcChemProps = list(
        advanced = TRUE,
        default = TRUE,
        description = "Prefer calculated chemical properties",
        type = "flag"
    ),
    neutralChemProps = list(
        advanced = TRUE,
        default = FALSE,
        description = "Use neutral chemical properties",
        type = "flag"
    ),
    virtualenv = list(
        advanced = TRUE,
        default = "patRoon-C3SDB",
        description = "Python virtual environment used for CCS prediction",
        type = "string",
        typeCheckArgs = list(min.chars = 1, null.ok = TRUE)
    )
))

AssignMobilitiesCompoundsParam <- setClass("AssignMobilitiesCompoundsParam", contains = "param")
setMethod("initialize", "AssignMobilitiesCompoundsParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "AssignMobilitiesCompoundsParam", baseName = "assign_mobilities_compounds", subDir = "annotation",
                   description = "Parameters for assignMobilities on compounds", version = "1.0",
                   definitions = getAssignMobilitiesCompoundsParamDefs(), ...)
})

setValidity("AssignMobilitiesCompoundsParam", function(object)
{
    parsFilled <- paramListFillDefaults(object@data, object@definitions)
    ac <- checkmate::makeAssertCollection()
    checkmate::assert(
        checkmate::checkChoice(parsFilled$from, c("pubchemlite", "c3sdb")),
        checkmate::checkDataFrame(parsFilled$from, min.rows = 1),
        checkmate::checkNull(parsFilled$from),
        .var.name = "from", add = ac
    )
    OK <- tryCatch(checkmate::reportAssertions(ac), error = function(e) e)
    return(OK)
})

#' @export
getPredictToxParamDefs <- paramConfigDefsFact(list(
    type = list(
        advanced = FALSE,
        default = "FP",
        description = "Prediction input type",
        type = "choice",
        typeCheckArgs = list(choices = c("FP", "SMILES", "both"))
    ),
    LC50Mode = list(
        advanced = FALSE,
        default = "static",
        description = "LC50 prediction mode",
        type = "choice",
        typeCheckArgs = list(choices = c("static", "flow"))
    ),
    concUnit = list(
        advanced = FALSE,
        default = "ugL",
        description = "Concentration unit",
        type = "concUnit"
    ),
    updateScore = list(
        advanced = FALSE,
        default = FALSE,
        description = "Update compound scores with predicted toxicity values (unused by SIRIUS)",
        type = "flag"
    ),
    scoreWeight = list(
        advanced = FALSE,
        default = 1,
        description = "Weight for the toxicity scoring (unused by SIRIUS)",
        type = "number",
        typeCheckArgs = list(lower = 1, finite = TRUE)
    ),
    parallel = list(
        advanced = TRUE,
        default = TRUE,
        description = "Process candidates in parallel (unused by SIRIUS)",
        type = "flag"
    )
))

#' @export
PredictToxParam <- setClass("PredictToxParam", contains = "param")
setMethod("initialize", "PredictToxParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "PredictToxParam", baseName = "predict_tox", subDir = "annotation",
                   description = "Parameters for predictTox", version = "1.0",
                   definitions = getPredictToxParamDefs(), ...)
})

#' @export
getPredictRespFactorsParamDefs <- paramConfigDefsFact(list(
    type = list(
        advanced = TRUE,
        default = "FP",
        description = "Prediction input type",
        type = "choice",
        typeCheckArgs = list(choices = c("FP", "SMILES", "both"))
    ),
    eluent = list(
        advanced = TRUE,
        default = NULL,
        description = "LC gradient program",
        type = "data.frame",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    organicModifier = list(
        advanced = TRUE,
        default = NULL,
        description = "Organic modifier of the mobile phase",
        type = "choice",
        typeCheckArgs = list(choices = c("MeOH", "MeCN"), null.ok = TRUE)
    ),
    pHAq = list(
        advanced = TRUE,
        default = NULL,
        description = "pH of the aqueous part of the mobile phase",
        type = "number",
        typeCheckArgs = list(finite = TRUE, null.ok = TRUE)
    ),
    concUnit = list(
        advanced = TRUE,
        default = "ugL",
        description = "Concentration unit",
        type = "concUnit"
    ),
    calibConcUnit = list(
        advanced = TRUE,
        default = "ugL",
        description = "Concentration unit used in the calibrants table",
        type = "concUnit"
    ),
    updateScore = list(
        advanced = TRUE,
        default = FALSE,
        description = "Update compound scores with predicted response factors (unused by SIRIUS)",
        type = "flag"
    ),
    scoreWeight = list(
        advanced = TRUE,
        default = 1,
        description = "Weight for the response-factor scoring (unused by SIRIUS)",
        type = "number",
        typeCheckArgs = list(lower = 1, finite = TRUE)
    ),
    parallel = list(
        advanced = TRUE,
        default = TRUE,
        description = "Process candidates in parallel (unused by SIRIUS)",
        type = "flag"
    )
))

#' @export
PredictRespFactorsParam <- setClass("PredictRespFactorsParam", contains = "param")
setMethod("initialize", "PredictRespFactorsParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "PredictRespFactorsParam", baseName = "predict_resp_factors", subDir = "annotation",
                   description = "Parameters for predictRespFactors", version = "1.0",
                   definitions = getPredictRespFactorsParamDefs(), ...)
})

