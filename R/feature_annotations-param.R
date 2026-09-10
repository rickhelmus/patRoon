# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
NULL

#' @export
getEstimateIDConfidenceParamDefs <- paramConfigDefsFact(list(
    absMzDev = list(
        default = defaultLim("mz", "medium"),
        description = "Maximum absolute m/z deviation",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    specSimParams = list(
        default = getDefSpecSimParams(removePrecursor = TRUE),
        description = "Spectrum similarity parameters",
        type = "specSimParams"
    ),
    checkFragments = list(
        default = c("mz", "formula", "compound"),
        description = "Types of MS/MS fragments for suspect fragment matches",
        type = "subset",
        typeCheckArgs = list(choices = c("mz", "formula", "compound"))
    ),
    formulasNormalizeScores = list(
        default = "max",
        description = "Formula score normalization method",
        type = "normalizationMethod",
        typeCheckArgs = list(withNone = FALSE)
    ),
    compoundsNormalizeScores = list(
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
    callNextMethod(.Object, name = "EstimateIDConfidenceParam", baseName = "EstimateIDConfidenceParam",
                   description = "Parameters for estimateIDConfidence", version = "1.0",
                   definitions = getEstimateIDConfidenceParamDefs(), ...)
})

#' @export
getAssignMobilitiesCompoundsParamDefs <- paramConfigDefsFact(list(
    IMS = list(
        default = TRUE,
        description = "Select IMS features for mobility assignment",
        type = "IMS"
    ),
    from = list(
        default = NULL,
        description = "Reference data for mobility and CCS assignment"
    ),
    matchFromBy = list(
        default = "InChIKey1",
        description = "Column used to match annotations to reference data",
        type = "choice",
        typeCheckArgs = list(choices = c("InChIKey", "InChIKey1", "InChI", "SMILES", "name"))
    ),
    overwrite = list(
        default = FALSE,
        description = "Overwrite existing mobility and CCS values",
        type = "flag"
    ),
    CCSParams = list(
        default = NULL,
        description = "Parameters for CCS assignment",
        type = "CCSParams",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    prefCalcChemProps = list(
        default = TRUE,
        description = "Prefer calculated chemical properties",
        type = "flag"
    ),
    neutralChemProps = list(
        default = FALSE,
        description = "Use neutral chemical properties",
        type = "flag"
    ),
    virtualenv = list(
        default = "patRoon-C3SDB",
        description = "Python virtual environment used for CCS prediction",
        type = "string",
        typeCheckArgs = list(min.chars = 1, null.ok = TRUE)
    )
))

AssignMobilitiesCompoundsParam <- setClass("AssignMobilitiesCompoundsParam", contains = "param")
setMethod("initialize", "AssignMobilitiesCompoundsParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "AssignMobilitiesCompoundsParam", baseName = "AssignMobilitiesCompoundsParam",
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
        default = "FP",
        description = "Prediction input type",
        type = "choice",
        typeCheckArgs = list(choices = c("FP", "SMILES", "both"))
    ),
    LC50Mode = list(
        default = "static",
        description = "LC50 prediction mode",
        type = "choice",
        typeCheckArgs = list(choices = c("static", "flow"))
    ),
    concUnit = list(
        default = "ugL",
        description = "Concentration unit",
        type = "concUnit"
    ),
    updateScore = list(
        default = FALSE,
        description = "Update compound scores with predicted toxicity values (unused by SIRIUS)",
        type = "flag"
    ),
    scoreWeight = list(
        default = 1,
        description = "Weight for the toxicity scoring (unused by SIRIUS)",
        type = "number",
        typeCheckArgs = list(lower = 1, finite = TRUE)
    ),
    parallel = list(
        default = TRUE,
        description = "Process candidates in parallel (unused by SIRIUS)",
        type = "flag"
    )
))

#' @export
PredictToxParam <- setClass("PredictToxParam", contains = "param")
setMethod("initialize", "PredictToxParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "PredictToxParam", baseName = "PredictToxParam",
                   description = "Parameters for predictTox", version = "1.0",
                   definitions = getPredictToxParamDefs(), ...)
})

#' @export
getPredictRespFactorsParamDefs <- paramConfigDefsFact(list(
    type = list(
        default = "FP",
        description = "Prediction input type",
        type = "choice",
        typeCheckArgs = list(choices = c("FP", "SMILES", "both"))
    ),
    eluent = list(
        default = NULL,
        description = "LC gradient program",
        type = "data.frame",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    organicModifier = list(
        default = NULL,
        description = "Organic modifier of the mobile phase",
        type = "choice",
        typeCheckArgs = list(choices = c("MeOH", "MeCN"), null.ok = TRUE)
    ),
    pHAq = list(
        default = NULL,
        description = "pH of the aqueous part of the mobile phase",
        type = "number",
        typeCheckArgs = list(finite = TRUE, null.ok = TRUE)
    ),
    concUnit = list(
        default = "ugL",
        description = "Concentration unit",
        type = "concUnit"
    ),
    calibConcUnit = list(
        default = "ugL",
        description = "Concentration unit used in the calibrants table",
        type = "concUnit"
    ),
    updateScore = list(
        default = FALSE,
        description = "Update compound scores with predicted response factors (unused by SIRIUS)",
        type = "flag"
    ),
    scoreWeight = list(
        default = 1,
        description = "Weight for the response-factor scoring (unused by SIRIUS)",
        type = "number",
        typeCheckArgs = list(lower = 1, finite = TRUE)
    ),
    parallel = list(
        default = TRUE,
        description = "Process candidates in parallel (unused by SIRIUS)",
        type = "flag"
    )
))

#' @export
PredictRespFactorsParam <- setClass("PredictRespFactorsParam", contains = "param")
setMethod("initialize", "PredictRespFactorsParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "PredictRespFactorsParam", baseName = "PredictRespFactorsParam",
                   description = "Parameters for predictRespFactors", version = "1.0",
                   definitions = getPredictRespFactorsParamDefs(), ...)
})

