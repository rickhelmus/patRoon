# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
NULL

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
    callNextMethod(.Object, name = "MSPeakListsParam", baseName = "mslists", subDir = "annotation",
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
        default = 1:2,
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
        default = NULL,
        description = "Minimum relative intensity threshold",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    topMostPeaks = list(
        default = NULL,
        description = "Only keep the top most intense peaks per spectrum",
        type = "count",
        typeCheckArgs = list(positive = TRUE, null.ok = TRUE)
    ),
    minPeaks = list(
        default = NULL,
        description = "Minimum number of peaks required in an MS/MS peak list",
        type = "count",
        typeCheckArgs = list(positive = TRUE, null.ok = TRUE)
    ),
    maxMZOverPrec = list(
        default = NULL,
        description = "Maximum m/z over precursor allowed in MS/MS",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    absMinAbundanceFeat = list(
        default = NULL,
        description = "Minimum absolute abundance of a peak across feature spectra",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    relMinAbundanceFeat = list(
        default = NULL,
        description = "Minimum relative abundance of a peak across feature spectra",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    absMinAbundanceFGroup = list(
        default = NULL,
        description = "Minimum absolute abundance of a peak across feature group spectra",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    relMinAbundanceFGroup = list(
        default = NULL,
        description = "Minimum relative abundance of a peak across feature group spectra",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    maxRelCumIntensity = list(
        default = NULL,
        description = "Maximum relative cumulative intensity of peaks",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, upper = 1, finite = TRUE)
    ),
    isolatePrec = list(
        default = NULL,
        description = "Precursor isolation parameters, or TRUE for default settings"
    ),
    deIsotope = list(
        default = FALSE,
        description = "Remove isotopic peaks",
        type = "flag"
    ),
    removeMZs = list(
        default = NULL,
        description = "m/z values or background peaks to remove"
    ),
    withMSMS = list(
        default = FALSE,
        description = "Keep only results with MS/MS data",
        type = "flag"
    ),
    annotatedBy = list(
        # UNDONE: move
        default = NULL,
        description = "Formula or compound annotations used to filter MS/MS peaks"
    ),
    retainPrecursor = list(
        default = TRUE,
        description = "Retain precursor peaks in MS/MS peak lists",
        type = "flag"
    ),
    mzWindow = list(
        default = defaultLim("mz", "medium"),
        description = "m/z window used to find peaks for the removeMZs filter",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    reAverage = list(
        default = FALSE,
        description = "Re-average feature group peak lists after filtering",
        type = "flag"
    ),
    negate = list(
        default = FALSE,
        description = "Apply filters in the opposite manner",
        type = "flag"
    )
))

#' @export
FilterMSPeakListsParam <- setClass("FilterMSPeakListsParam", contains = "param")
setMethod("initialize", "FilterMSPeakListsParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "FilterMSPeakListsParam", baseName = "filter_mslists",
                   subDir = "annotation", description = "Parameters for MSPeakLists filtering",
                   version = "1.0", definitions = getFilterMSPeakListsParamDefs(), ...)
})

setValidity("FilterMSPeakListsParam", function(object)
{
    parsFilled <- paramListFillDefaults(object@data, object@definitions)
    ac <- checkmate::makeAssertCollection()
    if (!is.logical(parsFilled$isolatePrec))
        assertPListIsolatePrecParams(parsFilled$isolatePrec, .var.name = "isolatePrec", add = ac)

    if (is.data.frame(parsFilled$removeMZs)) # output of getBGMSMSPeaks()
    {
        assertHasNames(parsFilled$removeMZs, "mz", add = ac)
        parsFilled$removeMZs <- parsFilled$removeMZs[["mz"]]
    }
    checkmate::assertNumeric(parsFilled$removeMZs, lower = 0, finite = TRUE, any.missing = FALSE, null.ok = TRUE,
                             add = ac)
    
    OK <- tryCatch(checkmate::reportAssertions(ac), error = function(e) e)
    return(OK)
})
