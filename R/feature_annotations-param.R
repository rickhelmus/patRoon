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

