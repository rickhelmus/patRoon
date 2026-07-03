# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

getTPsAnnCompParamDefs <- paramConfigDefsFact(list(
    minRTDiff = list(
        default = 20,
        description = "Minimum retention time difference between parent and TP",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    minFitFormula = list(
        default = 0.94,
        description = "Minimum fitFormula to filter out unlikely candidates",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    minFitCompound = list(
        default = 0,
        description = "Minimum fitCompound to filter out unlikely candidates",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    minSimSusp = list(
        default = 0,
        description = "Minimum similarity with suspects",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    minFitCompOrSimSusp = list(
        default = c(0.54, 0.65),
        description = "Thresholds for fitCompound OR simSusp",
        type = "numeric",
        typeCheckArgs = list(lower = 0, finite = TRUE, len = 2)
    ),
    extraOptsFMCSR = list(
        default = NULL,
        description = "Extra options passed to fmcsR",
        type = "list",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    skipInvalid = list(
        default = TRUE,
        description = "Skip parents without sufficient chemical information (e.g. SMILES)",
        type = "flag"
    ),
    prefCalcChemProps = list(
        default = TRUE,
        description = "If TRUE, prefer calculated chemical properties over already present data in the suspect list",
        type = "flag"
    ),
    neutralChemProps = list(
        default = FALSE,
        description = "If TRUE, ensure suspects are neutral",
        type = "flag"
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
    ),
    parallel = list(
        default = TRUE,
        description = "Enable parallel processing with futures",
        type = "flag"
    )
))

#' @export
TPsAnnCompParam <- setClass("TPsAnnCompParam", contains = "param")

setMethod("initialize", "TPsAnnCompParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "TPsAnnCompParam", baseName = "AnnComp",
                   description = "Parameters for annotation-based compound TP generation", version = "1.0",
                   definitions = getTPsAnnCompParamDefs(), ...)
})
