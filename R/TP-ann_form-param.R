# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

getTPsAnnFormParamDefs <- paramConfigDefsFact(list(
    minFitFormula = list(
        default = 0.94,
        description = "Minimum fitFormula to filter out unlikely candidates",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
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
    parallel = list(
        default = TRUE,
        description = "Enable parallel processing with futures",
        type = "flag"
    )
))

#' @export
TPsAnnFormParam <- setClass("TPsAnnFormParam", contains = "param")

setMethod("initialize", "TPsAnnFormParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "TPsAnnFormParam", baseName = "AnnForm",
                   description = "Parameters for annotation-based formula TP generation", version = "1.0",
                   definitions = getTPsAnnFormParamDefs(), ...)
})
