# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

getCompoundsMetFragParamDefs <- paramConfigDefsFact(list(
    specSimParams = list(
        default = getDefSpecSimParams(removePrecursor = TRUE),
        description = "Spectral similarity parameters used for annotation similarity calculations",
        type = "specSimParams"
    ),
    method = list(
        default = "CL",
        description = "Which method should be used for MetFrag execution: 'CL' (recommended) or 'R'.",
        type = "choice",
        typeCheckArgs = list(choices = c("CL", "R"))
    ),
    timeout = list(
        default = 300,
        description = "Maximum time (in seconds) before a MetFrag query for a feature group is stopped",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    timeoutRetries = list(
        default = 2,
        description = "Maximum number of retries after reaching a timeout",
        type = "count"
    ),
    errorRetries = list(
        default = 2,
        description = "Maximum number of retries after an error occurred (e.g. connection errors).",
        type = "count"
    ),
    topMost = list(
        default = 100,
        description = "Maximum number of top ranked candidates to keep",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    maxCandidatesToStop = list(
        default = 2500,
        description = "If more than this number of candidate structures are found then processing will be aborted",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    dbRelMzDev = list(
        default = defaultLim("mz", "narrow_rel"),
        description = "Relative mass deviation (in ppm) for database search",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    fragRelMzDev = list(
        default = defaultLim("mz", "narrow_rel"),
        description = "Relative mass deviation (in ppm) for fragment matching",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    fragAbsMzDev = list(
        default = defaultLim("mz", "narrow"),
        description = "Absolute mass deviation for fragment matching",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    database = list(
        default = "pubchemlite",
        description = "Compound database to use",
        type = "choice",
        typeCheckArgs = list(choices = c("pubchem", "chemspider", "kegg", "for-ident", "sdf", "psv",
                                         "csv", "comptox", "pubchemlite"))
    ),
    extendedPubChem = list(
        default = "auto",
        description = "Whether to use the extended PubChem database"
    ),
    chemSpiderToken = list(
        default = "",
        description = "ChemSpider security token (required when database='chemspider')",
        type = "string"
    ),
    scoreTypes = list(
        default = NULL,
        description = "Character vector defining the scoring types (set to NULL for default scoring types)",
        type = "character",
        typeCheckArgs = list(any.missing = FALSE, null.ok = TRUE)
    ),
    scoreWeights = list(
        default = 1.0,
        description = "Numeric vector containing weights of the used scoring types",
        type = "numeric",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    preProcessingFilters = list(
        default = c("UnconnectedCompoundFilter", "IsotopeFilter"),
        description = "A character vector defining pre-filters applied before fragmentation and scoring",
        type = "character",
        typeCheckArgs = list(any.missing = FALSE)
    ),
    postProcessingFilters = list(
        default = "InChIKeyFilter",
        description = "A character vector defining post-filters applied after fragmentation and scoring",
        type = "character",
        typeCheckArgs = list(any.missing = FALSE)
    ),
    minIMSSpecSim = list(
        default = 0,
        description = "Minimum spectral similarity to copy annotations from the IMS precursor",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    setThreshold = list(
        default = 0,
        description = "Minimum abundance for a candidate among all sets",
        type = "number",
        typeCheckArgs = list(lower = 0, upper = 1, finite = TRUE)
    ),
    setThresholdAnn = list(
        default = 0,
        description = "As setThreshold, but only taking sets into account in which features of the feature group are present",
        type = "number",
        typeCheckArgs = list(lower = 0, upper = 1, finite = TRUE)
    ),
    setAvgSpecificScores = list(
        default = FALSE,
        description = "If TRUE, average set specific scores (e.g. MS/MS match)",
        type = "flag"
    ),
    extraOpts = list(
        default = NULL,
        description = "Extra options passed to MetFrag as a named character list",
        type = "list",
        typeCheckArgs = list(names = "unique", any.missing = FALSE, null.ok = TRUE)
    )
))

#' @export
CompoundsMetFragParam <- setClass("CompoundsMetFragParam", contains = "param")

setMethod("initialize", "CompoundsMetFragParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "CompoundsMetFragParam", baseName = "CompoundsMetFragParam",
                   description = "Parameters for MetFrag compound generation", version = "1.0",
                   definitions = getCompoundsMetFragParamDefs(), ...)
})

setValidity("CompoundsMetFragParam", function(object)
{
    parsFilled <- paramListFillDefaults(object@data, object@definitions)
    ac <- checkmate::makeAssertCollection()
    checkmate::assert(checkmate::checkFlag(parsFilled$extendedPubChem),
                      checkmate::checkChoice(parsFilled$extendedPubChem, "auto"),
                      .var.name = "extendedPubChem", add = ac)
    OK <- tryCatch(checkmate::reportAssertions(ac), error = function(e) e)
    return(OK)
})
