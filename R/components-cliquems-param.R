# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

#' @export
getComponentsCliqueMSParamDefs <- paramConfigDefsFact(list(
    maxCharge = list(
        default = 1,
        description = "Maximum charge to consider (passed to cliqueMS::getIsotopes)",
        type = "count",
        positive = TRUE
    ),
    maxGrade = list(
        default = 2,
        description = "Maximum number of isotopes apart from the monoisotope (passed to cliqueMS::getIsotopes)",
        type = "count",
        positive = TRUE
    ),
    ppm = list(
        default = 10,
        description = "Mass accuracy tolerance (ppm) (passed to cliqueMS functions)",
        type = "number",
        lower = 0,
        finite = TRUE
    ),
    adductInfo = list(
        default = NULL,
        description = "Adduct information data.frame passed to cliqueMS::getAnnotation; NULL to use package defaults",
        type = "data.frame",
        null.ok = TRUE
    ),
    absMzDev = list(
        default = defaultLim("mz", "medium"),
        description = "Absolute m/z deviation used for componentization",
        type = "number",
        finite = TRUE,
        lower = 0
    ),
    minSize = list(
        default = 2,
        description = "Minimum number of feature groups to form a component",
        type = "count",
        positive = TRUE
    ),
    relMinAdductAbundance = list(
        default = 0.75,
        description = "Relative threshold that an adduct should be assigned to features within the same feature group",
        type = "number",
        finite = TRUE,
        lower = 0
    ),
    adductConflictsUsePref = list(
        default = TRUE,
        description = "Resolve adduct conflicts using preferred adducts when possible",
        type = "flag"
    ),
    NMConflicts = list(
        default = c("preferential", "mostAbundant", "mostIntense"),
        description = "Strategies to resolve neutral mass conflicts",
        type = "subset",
        choices = c("preferential", "mostAbundant", "mostIntense"),
        empty.ok = FALSE
    ),
    prefAdducts = list(
        default = c("[M+H]+", "[M-H]-"),
        description = "Preferred adducts used to resolve conflicts",
        type = "character",
        min.chars = 1,
        any.missing = FALSE,
        unique = TRUE
    ),
    extraOptsCli = list(
        default = NULL,
        description = "Extra options passed to cliqueMS::getCliques as a named character list",
        type = "list",
        names = "unique",
        any.missing = FALSE,
        null.ok = TRUE
    ),
    extraOptsIso = list(
        default = NULL,
        description = "Extra options passed to cliqueMS::getIsotopes as a named character list",
        type = "list",
        names = "unique",
        any.missing = FALSE,
        null.ok = TRUE
    ),
    extraOptsAnn = list(
        default = NULL,
        description = "Extra options passed to cliqueMS::getAnnotation as a named character list",
        type = "list",
        names = "unique",
        any.missing = FALSE,
        null.ok = TRUE
    ),
    parallel = list(
        default = TRUE,
        description = "Run per-analysis cliqueMS calls in parallel",
        type = "flag"
    )
))

#' @export
ComponentsCliqueMSParam <- setClass("ComponentsCliqueMSParam", contains = "param")

setMethod("initialize", "ComponentsCliqueMSParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "ComponentsCliqueMSParam", baseName = "ComponentsCliqueMSParam",
                   description = "Parameters for cliqueMS component generation", version = "1.0",
                   definitions = getComponentsCliqueMSParamDefs(), ...)
})
