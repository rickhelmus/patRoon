# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
NULL

#' @export
getComponentsOpenMSParamDefs <- paramConfigDefsFact(list(
    chargeMin = list(
        default = 1,
        description = "Minimum charge to consider",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    chargeMax = list(
        default = 1,
        description = "Maximum charge to consider",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    chargeSpan = list(
        default = 3,
        description = "Charge span considered when generating adduct candidates",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    qTry = list(
        default = "heuristic",
        description = "Strategy to try candidate adduct assignments",
        type = "choice",
        typeCheckArgs = list(choices = c("heuristic", "all"))
    ),
    # UNDONE: this needs to be updated to work for sets, see TODO.md
    potentialAdducts = list(
        default = NULL,
        description = "Adducts to consider with probabilities, eg c(\"[M+H]+\" = 0.8, \"[M+Na]+\" = 0.2)",
        type = "list",
        typeCheckArgs = list(types = "numeric", names = "unique", null.ok = TRUE)
    ),
    minRTOverlap = list(
        default = 0.66,
        description = "Minimum relative retention time overlap",
        type = "number",
        typeCheckArgs = list(lower = 0, upper = 1)
    ),
    retWindow = list(
        default = defaultLim("retention", "very_narrow"),
        description = "Retention time tolerance",
        type = "number",
        typeCheckArgs = list(finite = TRUE, lower = 0)
    ),
    absMzDev = list(
        default = defaultLim("mz", "medium"),
        description = "Absolute m/z deviation threshold used for componentization",
        type = "number",
        typeCheckArgs = list(finite = TRUE, lower = 0)
    ),
    minSize = list(
        default = 2,
        description = "Minimum number of feature groups to form a component",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    relMinAdductAbundance = list(
        default = 0.75,
        description = "Relative threshold that an adduct should be assigned to features within the same feature group",
        type = "number",
        typeCheckArgs = list(finite = TRUE, lower = 0)
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
        typeCheckArgs = list(choices = c("preferential", "mostAbundant", "mostIntense"), empty.ok = FALSE)
    ),
    prefAdducts = list(
        default = c("[M+H]+", "[M-H]-"),
        description = "Preferred adducts used to resolve conflicts",
        type = "character",
        typeCheckArgs = list(min.chars = 1, any.missing = FALSE, unique = TRUE)
    ),
    extraOpts = list(
        default = NULL,
        description = "Extra command line options passed directly to MetaboliteAdductDecharger as a named character list",
        type = "list",
        typeCheckArgs = list(types = "character", any.missing = FALSE, names = "unique", null.ok = TRUE)
    )
))

#' @export
ComponentsOpenMSParam <- setClass("ComponentsOpenMSParam", contains = "param")

setMethod("initialize", "ComponentsOpenMSParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "ComponentsOpenMSParam", baseName = "ComponentsOpenMSParam",
                   description = "Parameters for OpenMS component generation", version = "1.0",
                   definitions = getComponentsOpenMSParamDefs(), ...)
})
