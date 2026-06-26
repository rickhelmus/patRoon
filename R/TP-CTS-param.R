# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

getTPsCTSParamDefs <- paramConfigDefsFact(list(
    transLibrary = list(
        default = "combined_photolysis_abiotic_hydrolysis",
        description = "The transformation library to use.",
        type = "choice",
        typeCheckArgs = list(choices = c("hydrolysis", "abiotic_reduction", "photolysis_unranked", 
                                         "photolysis_ranked", "mammalian_metabolism",
                                         "combined_abioticreduction_hydrolysis",
                                         "combined_photolysis_abiotic_hydrolysis",
                                         "pfas_environmental", "pfas_metabolism"))
    ),
    generations = list(
        default = 1L,
        description = "The number of transformation generations to predict",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    errorRetries = list(
        default = 3L,
        description = "The maximum number of connection retries",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    neutralizeTPs = list(
        default = TRUE,
        description = "If TRUE, neutralize predicted TPs",
        type = "flag"
    ),
    TPStructParams = list(
        default = getDefTPStructParams(),
        description = "Parameters for the calculation of TP structure properties",
        type = "TPStructParams"
    )
))

#' @export
TPsCTSParam <- setClass("TPsCTSParam", contains = "param")

setMethod("initialize", "TPsCTSParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "TPsCTSParam", baseName = "CTS",
                   description = "Parameters for CTS TP generation", version = "1.0",
                   definitions = getTPsCTSParamDefs(), ...)
})
