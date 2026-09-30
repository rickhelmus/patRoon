# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

#' @export
getComponentsNetParamDefs <- paramConfigDefsFact(list(
    minSize = list(
        default = 2,
        description = "Minimum number of feature groups to form a component",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    mzWindow = list(
        default = defaultLim("mz", "medium"),
        description = "Absolute m/z tolerance used for annotation and MS2 matching",
        type = "number",
        typeCheckArgs = list(finite = TRUE, lower = 0)
    ),
    componSim = list(
        default = "pearson",
        description = "Similarity metric used for feature networking",
        type = "choice",
        typeCheckArgs = list(choices = c("pearson", "cosine"))
    ),
    componMinSim = list(
        default = 0.95,
        description = "Minimum feature similarity to connect nodes",
        type = "number",
        typeCheckArgs = list(finite = TRUE, lower = 0, upper = 1)
    ),
    componMaxP = list(
        default = 0.05,
        description = "Maximum p-value for Pearson similarity",
        type = "number",
        typeCheckArgs = list(finite = TRUE, lower = 0, upper = 1)
    ),
    componMethod = list(
        default = "community",
        description = "Network componentization method",
        type = "choice",
        typeCheckArgs = list(choices = c("community", "cliques", "hcs", "hclust"))
    ),
    componArgs = list(
        default = list(),
        description = "Additional arguments passed to the network componentization method",
        type = "list",
        typeCheckArgs = list(any.missing = FALSE, names = "unique")
    ),
    groupClust = list(
        default = "complete",
        description = "Clustering method for consensus components across analyses",
        type = "string"
    ),
    groupClustH = list(
        default = 0.5,
        description = "Height at which to cut the consensus clustering tree",
        type = "number",
        typeCheckArgs = list(finite = TRUE, lower = 0, upper = 1)
    ),
    annotAlgo = list(
        default = "imss",
        description = "Annotation algorithm",
        type = "choice",
        typeCheckArgs = list(choices = c("imss", "nontarget"))
    ),
    annotAdducts = list(
        default = c("[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+", "[M-H]-", "[M-H2O-H]-"),
        description = "Adducts to use for component annotation"
    ),
    annotPrefAdducts = list(
        default = c("[M+H]+", "[M-H]-"),
        description = "Preferred adducts for nontarget annotation",
        type = "character",
        typeCheckArgs = list(min.chars = 1, any.missing = FALSE, unique = TRUE)
    ),
    annotArgs = list(
        default = list(),
        description = "Additional arguments passed to the annotation function",
        type = "list",
        typeCheckArgs = list(any.missing = FALSE, names = "unique")
    ),
    mzFragBelow = list(
        default = 10,
        description = "Minimum m/z difference below a precursor for MS2 matching",
        type = "number",
        typeCheckArgs = list(finite = TRUE, lower = 0)
    )
))

#' @export
ComponentsNetParam <- setClass("ComponentsNetParam", contains = "param")

setMethod("initialize", "ComponentsNetParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "ComponentsNetParam", baseName = "ComponentsNetParam",
                   description = "Parameters for network-based component generation", version = "1.0",
                   definitions = getComponentsNetParamDefs(), ...)
})

setValidity("ComponentsNetParam", function(object)
{
    parsFilled <- paramListFillDefaults(object@data, object@definitions)
    ac <- checkmate::makeAssertCollection()
    checkmate::assert(
        checkmate::checkCharacter(parsFilled$annotAdducts, min.chars = 1, any.missing = FALSE, min.len = 2),
        checkmate::checkList(parsFilled$annotAdducts, types = "adduct", min.len = 2, any.missing = FALSE),
        .var.name = "annotAdducts", add = ac
    )
    OK <- tryCatch(checkmate::reportAssertions(ac), error = function(e) e)
    return(OK)
})
