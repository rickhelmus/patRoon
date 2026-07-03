# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

getTPsLibraryFormulaParamDefs <- paramConfigDefsFact(list(
    generations = list(
        default = 1L,
        description = "The number of transformation generations to predict",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    skipInvalid = list(
        default = TRUE,
        description = "Skip invalid parents without sufficient chemical information (e.g. formulae)",
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
    matchParentsBy = list(
        default = "name",
        description = "How to match parents in the library",
        type = "choice",
        typeCheckArgs = list(choices = c("InChIKey", "InChIKey1", "InChI", "SMILES", "formula", "name"))
    ),
    matchGenerationsBy = list(
        default = "name",
        description = "How to match generations in the library",
        type = "choice",
        typeCheckArgs = list(choices = c("InChIKey", "InChIKey1", "InChI", "SMILES", "formula", "name"))
    )
))

#' @export
TPsLibraryFormulaParam <- setClass("TPsLibraryFormulaParam", contains = "param")

setMethod("initialize", "TPsLibraryFormulaParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "TPsLibraryFormulaParam", baseName = "LibraryFormula",
                   description = "Parameters for library-based formula TP generation", version = "1.0",
                   definitions = getTPsLibraryFormulaParamDefs(), ...)
})
