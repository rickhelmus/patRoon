# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

getTPsBioTransformerParamDefs <- paramConfigDefsFact(list(
    type = list(
        default = "env",
        description = "The type of prediction.",
        type = "choice",
        typeCheckArgs = list(choices = c("ecbased", "cyp450", "phaseII", "hgut", "superbio", "allHuman", "env"))
    ),
    generations = list(
        default = 2L,
        description = "The number of generations (steps) for the predictions",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    maxExpGenerations = list(
        default = 2L,
        description = "The maximum number of additional generations during hierarchy expansion",
        type = "count"
    ),
    extraOpts = list(
        default = NULL,
        description = "Extra command line options passed to the biotransformer.jar tool",
        type = "character",
        typeCheckArgs = list(min.chars = 1, null.ok = TRUE)
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
    MP = list(
        default = FALSE,
        description = "Enables multiprocessing (generally not recommended)",
        type = "flag"
    )
))

#' @export
TPsBioTransformerParam <- setClass("TPsBioTransformerParam", contains = "param")

setMethod("initialize", "TPsBioTransformerParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "TPsBioTransformerParam", baseName = "BioTransformer",
                   description = "Parameters for BioTransformer TP generation", version = "1.0",
                   definitions = getTPsBioTransformerParamDefs(), ...)
})
