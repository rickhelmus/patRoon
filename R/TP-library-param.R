# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

getTPsLibraryParamDefs <- paramConfigDefsFact(list(
    generations = list(
        default = 1L,
        description = "The number of transformation generations to predict",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    skipInvalid = list(
        default = TRUE,
        description = "Skip invalid parents",
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
        default = FALSE,
        description = "If TRUE, neutralize predicted TPs",
        type = "flag"
    ),
    matchParentsBy = list(
        default = "InChIKey",
        description = "How to match parents in the library",
        type = "choice",
        typeCheckArgs = list(choices = c("InChIKey", "InChIKey1", "InChI", "SMILES", "formula", "name"))
    ),
    matchGenerationsBy = list(
        default = "InChIKey",
        description = "How to match generations in the library",
        type = "choice",
        typeCheckArgs = list(choices = c("InChIKey", "InChIKey1", "InChI", "SMILES", "formula", "name"))
    ),
    TPStructParams = list(
        default = getDefTPStructParams(),
        description = "Parameters for the calculation of TP structure properties",
        type = "TPStructParams"
    )
))

#' @export
TPsLibraryParam <- setClass("TPsLibraryParam", contains = "param")

setMethod("initialize", "TPsLibraryParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "TPsLibraryParam", baseName = "Library",
                   description = "Parameters for library-based TP generation", version = "1.0",
                   definitions = getTPsLibraryParamDefs(), ...)
})
