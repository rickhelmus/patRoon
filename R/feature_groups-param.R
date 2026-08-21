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


#' @export
getNormIntsParamDefs <- paramConfigDefsFact(list(
    featNorm = list(
        default = "none",
        description = "Feature normalization method",
        type = "choice",
        typeCheckArgs = list(choices = c("tic", "istd", "conc", "none"))
    ),
    groupNorm = list(
        default = FALSE,
        description = "Perform normalization on feature groups",
        type = "flag"
    ),
    normFunc = list(
        default = max,
        description = "Normalization function",
        type = "function"
    ),
    ISTDRTWindow = list(
        default = 120,
        description = "RT window for ISTD matching",
        type = "range",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    ISTDMZWindow = list(
        default = 300,
        description = "m/z window for ISTD matching",
        type = "range",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    minISTDs = list(
        default = 3,
        description = "Minimum number of ISTDs to consider",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    extraOpts = list(
        default = list(),
        description = "Extra options passed to screenSuspects() as a named list",
        type = "list",
        typeCheckArgs = list(null.ok = TRUE, any.missing = FALSE, names = "unique")
    )
))

NormIntsParam <- setClass("NormIntsParam", contains = "param")
setMethod("initialize", "NormIntsParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "NormIntsParam", baseName = "NormIntsParam",
                   description = "Parameters for normInts", version = "1.0",
                   definitions = getNormIntsParamDefs(), ...)
})


#' @export
getCalculatePeakQualitiesParamDefs <- paramConfigDefsFact(list(
    weights = list(
        default = NULL,
        description = "Weighting for quality scoring",
        type = "list",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    flatnessFactor = list(
        default = 0.05,
        description = "Flatness factor used by MetaClean",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0)
    ),
    featureQualities = list(
        default = NULL,
        description = "Subset or custom feature qualities",
        type = "list",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    EICParams = list(
        default = getDefEICParams(window = 0),
        description = "EIC parameters",
        type = "list"
    ),
    featureGroupQualities = list(
        default = NULL,
        description = "Subset or custom feature group qualities",
        type = "list",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    avgFunc = list(
        default = mean,
        description = "Function used to average feature qualities into group qualities",
        type = "function"
    ),
    parallel = list(
        default = TRUE,
        description = "Parallel processing flag",
        type = "flag"
    )
))

CalculatePeakQualitiesParam <- setClass("CalculatePeakQualitiesParam", contains = "param")
setMethod("initialize", "CalculatePeakQualitiesParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "CalculatePeakQualitiesParam", baseName = "CalculatePeakQualitiesParam",
                   description = "Parameters for calculatePeakQualities", version = "1.0",
                   definitions = getCalculatePeakQualitiesParamDefs(), ...)
})
