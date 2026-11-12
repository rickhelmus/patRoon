# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
NULL

#' @export
getFilterFeatGroupsParamDefs <- paramConfigDefsFact(list(
    absMinIntensity = list(
        advanced = TRUE,
        default = NULL,
        description = "Minimum absolute intensity threshold",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    relMinIntensity = list(
        advanced = TRUE,
        default = NULL,
        description = "Minimum relative intensity threshold",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    preAbsMinIntensity = list(
        advanced = TRUE,
        default = NULL,
        description = "Pre-filter minimum absolute intensity threshold",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    preRelMinIntensity = list(
        advanced = TRUE,
        default = NULL,
        description = "Pre-filter minimum relative intensity threshold",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    absMinMaxIntensity = list(
        advanced = TRUE,
        default = NULL,
        description = "Minimum absolute intensity threshold for the most intense feature in a feature group",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    relMinMaxIntensity = list(
        advanced = TRUE,
        default = NULL,
        description = "Minimum relative intensity threshold for the most intense feature in a feature group",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    absMinAnalyses = list(
        advanced = TRUE,
        default = NULL,
        description = "Minimum number of analyses with feature data in a feature group",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    relMinAnalyses = list(
        advanced = TRUE,
        default = NULL,
        description = "Minimum relative number of analyses with feature data in a feature group",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    absMinReplicates = list(
        advanced = TRUE,
        default = NULL,
        description = "Minimum number of replicates with feature group data",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    relMinReplicates = list(
        advanced = TRUE,
        default = NULL,
        description = "Minimum relative number of replicates with feature group data",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    absMinFeatures = list(
        advanced = TRUE,
        default = NULL,
        description = "Minimum number of features per analysis",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    relMinFeatures = list(
        advanced = TRUE,
        default = NULL,
        description = "Minimum relative number of features per analysis",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    absMinReplicateAbundance = list(
        advanced = TRUE,
        default = NULL,
        description = "Minimum absolute abundance within a replicate",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    relMinReplicateAbundance = list(
        advanced = TRUE,
        default = NULL,
        description = "Minimum relative abundance within a replicate",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    absMinConc = list(
        advanced = TRUE,
        default = NULL,
        description = "Minimum absolute predicted concentration",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    relMinConc = list(
        advanced = TRUE,
        default = NULL,
        description = "Minimum relative predicted concentration",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    absMaxTox = list(
        advanced = TRUE,
        default = NULL,
        description = "Maximum absolute predicted toxicity",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    relMaxTox = list(
        advanced = TRUE,
        default = NULL,
        description = "Maximum relative predicted toxicity",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    absMinConcTox = list(
        advanced = TRUE,
        default = NULL,
        description = "Minimum absolute concentration-toxicity ratio",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    relMinConcTox = list(
        advanced = TRUE,
        default = NULL,
        description = "Minimum relative concentration-toxicity ratio",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    maxReplicateIntRSD = list(
        advanced = TRUE,
        default = NULL,
        description = "Maximum intensity relative standard deviation within a replicate",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    maxReplicateIntRSDPres = list(
        advanced = TRUE,
        default = NULL,
        description = "Maximum intensity relative standard deviation within a replicate, ignoring absent features",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    blankThreshold = list(
        advanced = TRUE,
        default = NULL,
        description = "Blank subtraction threshold relative to blank intensity",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    retentionRange = list(
        advanced = TRUE,
        default = NULL,
        description = "Retention time range (min, max)",
        type = "range",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    mzRange = list(
        default = NULL,
        description = "m/z range (min, max)",
        type = "range",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    mzDefectRange = list(
        advanced = TRUE,
        default = NULL,
        description = "m/z defect range (min, max)",
        type = "range",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    chromWidthRange = list(
        advanced = TRUE,
        default = NULL,
        description = "Chromatographic width range (min, max)",
        type = "range",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    # NOTE/UNDONE: we don't pass score names (NULL) here, since that needs fGroups data.
    featQualityRange = list(
        advanced = TRUE,
        default = NULL,
        description = "Feature quality or score ranges",
        type = "scoreRange",
        typeCheckArgs = list(scNames = NULL)
    ),
    groupQualityRange = list(
        advanced = TRUE,
        default = NULL,
        description = "Feature group quality or score ranges",
        type = "scoreRange",
        typeCheckArgs = list(scNames = NULL)
    ),
    replicates = list(
        advanced = TRUE,
        default = NULL,
        description = "Replicates to keep",
        type = "character",
        typeCheckArgs = list(null.ok = TRUE, min.len = 1, min.chars = 1, any.missing = FALSE)
    ),
    IMS = list(
        advanced = TRUE,
        default = NULL,
        description = "IMS filter selection",
        type = "IMS",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    withIMSPrecursor = list(
        advanced = TRUE,
        default = FALSE,
        description = "Keep only IMS feature groups with IMS precursors",
        type = "flag"
    ),
    IMSRangeParams = list(
        advanced = TRUE,
        default = NULL,
        description = "Parameters for filtering by ion mobility range",
        type = "IMSRangeParams",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    removeBlanks = list(
        advanced = TRUE,
        default = FALSE,
        description = "Remove blank analyses after filtering",
        type = "flag"
    ),
    removeISTDs = list(
        advanced = TRUE,
        default = FALSE,
        description = "Remove internal standard feature groups",
        type = "flag"
    ),
    checkFeaturesSession = list(
        advanced = TRUE,
        default = NULL,
        description = "Feature checking session file or TRUE to use the default session file"
    ),
    predAggrParams = list(
        advanced = TRUE,
        default = getDefPredAggrParams(),
        description = "Parameters for aggregating predicted concentrations and toxicities",
        type = "predAggrParams"
    ),
    removeNA = list(
        advanced = TRUE,
        default = FALSE,
        description = "Remove feature groups with missing aggregated predictions",
        type = "flag"
    ),
    sets = list(
        advanced = TRUE,
        default = NULL,
        description = "Sets to keep",
        type = "sets",
        typeCheckArgs = list(multiple = TRUE)
    ),
    absMinSets = list(
        advanced = TRUE,
        default = NULL,
        description = "Minimum absolute number of sets with feature group data",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    relMinSets = list(
        advanced = TRUE,
        default = NULL,
        description = "Minimum relative number of sets with feature group data",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    negate = list(
        advanced = TRUE,
        default = FALSE,
        description = "Apply filters in the opposite manner",
        type = "flag"
    ),
    applyIMS = list(
        advanced = TRUE,
        default = "both",
        description = "How filters are applied to IMS feature groups",
        type = "applyIMS"
    )
))

#' @export
FilterFeatGroupsParam <- setClass("FilterFeatGroupsParam", contains = "param")
setMethod("initialize", "FilterFeatGroupsParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "FilterFeatGroupsParam", baseName = "filter_featgroups", subDir = "featgroups",
                   description = "Parameters for featureGroups filtering", version = "1.0",
                   definitions = getFilterFeatGroupsParamDefs(), ...)
})

setValidity("FilterFeatGroupsParam", function(object)
{
    parsFilled <- paramListFillDefaults(object@data, object@definitions)
    ac <- checkmate::makeAssertCollection()
    if (!is.logical(parsFilled$checkFeaturesSession))
        assertCheckSession(parsFilled$checkFeaturesSession, mustExist = TRUE, null.ok = TRUE,
                           .var.name = "checkFeaturesSession", add = ac)
    OK <- tryCatch(checkmate::reportAssertions(ac), error = function(e) e)
    return(OK)
})
