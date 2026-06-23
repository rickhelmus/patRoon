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
