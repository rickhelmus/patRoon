# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
NULL

#' @export
getMSPeakListsParamDefs <- paramConfigDefsFact(list(
    maxMSRTWindow = list(
        default = defaultLim("retention", "narrow"),
        description = "Maximum retention time window for MS peak list generation",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    fixedIsolationWidth = list(
        default = FALSE,
        description = "Configures how MS/MS spectra are selected for a feature (FALSE, NA, or numeric tolerance)"
    ),
    topMost = list(
        default = NULL,
        description = "Only extract MS peak lists from a maximum of topMost features (NULL for all)",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0)
    ),
    avgFeatParams = list(
        default = getDefAvgPListParams(),
        description = "Parameters for averaging MS peak lists of individual features",
        type = "list"
    ),
    avgFGroupParams = list(
        default = getDefAvgPListParams(),
        description = "Parameters for averaging MS peak lists of feature groups",
        type = "list"
    )
))

#' @export
MSPeakListsParam <- setClass("MSPeakListsParam", contains = "param")
setMethod("initialize", "MSPeakListsParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "MSPeakListsParam", baseName = "MSPeakListsParam",
                   description = "Parameters for generateMSPeakLists", version = "1.0",
                   definitions = getMSPeakListsParamDefs(), ...)
})

setValidity("MSPeakListsParam", function(object)
{
    parsFilled <- paramListFillDefaults(object@data, object@definitions)
    
    ac <- checkmate::makeAssertCollection()
    assertAvgPListParams(parsFilled$avgFeatParams, .var.name = "avgFeatParams", add = ac)
    assertAvgPListParams(parsFilled$avgFGroupParams, .var.name = "avgFGroupParams", add = ac)
    checkmate::assert(
        checkmate::checkFALSE(parsFilled$fixedIsolationWidth),
        checkmate::checkScalarNA(parsFilled$fixedIsolationWidth),
        checkmate::checkNumber(parsFilled$fixedIsolationWidth, lower = 0, finite = TRUE),
        .var.name = "fixedIsolationWidth", add = ac
    )
    OK <- tryCatch(checkmate::reportAssertions(ac), error = function(e) e)
    
    return(OK)
})


#' @export
getFilterMSPeakListsParamDefs <- paramConfigDefsFact(list(
    MSLevel = list(
        default = 2,
        description = "MS levels to filter (1 or 2)",
        type = "subset",
        typeCheckArgs = list(choices = 1:2)
    ),
    absMinIntensity = list(
        default = NULL,
        description = "Absolute minimum intensity threshold",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    relMinIntensity = list(
        default = 0.05,
        description = "Minimum relative intensity threshold",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    topMostPeaks = list(
        default = 25,
        description = "Only keep the top most intense peaks per spectrum (NULL to keep all)",
        type = "count",
        typeCheckArgs = list(positive = TRUE, null.ok = TRUE)
    ),
    maxMZOverPrec = list(
        default = 4,
        description = "Maximum m/z over precursor (relative) allowed in MS/MS",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    extraOpts = list(
        default = list(),
        description = "Extra filter options as a named list",
        type = "list",
        typeCheckArgs = list(null.ok = TRUE, any.missing = FALSE, names = "unique")
    )
))

#' @export
FilterMSPeakListsParam <- setClass("FilterMSPeakListsParam", contains = "param")
setMethod("initialize", "FilterMSPeakListsParam", function(.Object, ...)
{
    args <- list(...)
    defs <- getFilterMSPeakListsParamDefs()
    dotsConstr <- args[names(args) %in% names(defs)]
    dotsConstr$extraOpts <- args[setdiff(names(args), names(defs))]

    do.call(callNextMethod, c(list(.Object, name = "FilterMSPeakListsParam", baseName = "FilterMSPeakListsParam",
                                   description = "Parameters for MSPeakLists filtering", version = "1.0",
                                   definitions = defs), dotsConstr))
})
