# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
NULL

#' @export
getFilterFeatGroupsScreeningParamDefs <- paramConfigDefsFact(list(
    onlyHits = list(
        default = NULL,
        description = "Only keep feature groups with suspect hits",
        type = "flag",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    IMSMatchParams = list(
        advanced = TRUE,
        default = NULL,
        description = "IMS matching parameters",
        type = "IMSMatchParams",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    selectHitsBy = list(
        advanced = TRUE,
        default = NULL,
        description = "Select best suspect hits by intensity or identification level",
        type = "choice",
        typeCheckArgs = list(choices = c("intensity", "level"), null.ok = TRUE)
    ),
    selectBestFGroups = list(
        advanced = TRUE,
        default = FALSE,
        description = "Keep only the best feature group for each suspect",
        type = "flag"
    ),
    maxLevel = list(
        advanced = TRUE,
        default = NULL,
        description = "Maximum identification level",
        type = "count",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    maxFormRank = list(
        advanced = TRUE,
        default = NULL,
        description = "Maximum formula rank",
        type = "count",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    maxCompRank = list(
        advanced = TRUE,
        default = NULL,
        description = "Maximum compound rank",
        type = "count",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    minAnnSimForm = list(
        advanced = TRUE,
        default = NULL,
        description = "Minimum formula annotation similarity",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    minAnnSimComp = list(
        advanced = TRUE,
        default = NULL,
        description = "Minimum compound annotation similarity",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    minAnnSimBoth = list(
        advanced = TRUE,
        default = NULL,
        description = "Minimum combined annotation similarity",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    absMinFragMatches = list(
        advanced = TRUE,
        default = NULL,
        description = "Minimum absolute number of matched fragments",
        type = "count",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    relMinFragMatches = list(
        advanced = TRUE,
        default = NULL,
        description = "Minimum relative number of matched fragments",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    minRF = list(
        advanced = TRUE,
        default = NULL,
        description = "Minimum predicted response factor",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    maxLC50 = list(
        advanced = TRUE,
        default = NULL,
        description = "Maximum predicted LC50",
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
FilterFeatGroupsScreeningParam <- setClass("FilterFeatGroupsScreeningParam", contains = "param")
setMethod("initialize", "FilterFeatGroupsScreeningParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "FilterFeatGroupsScreeningParam", baseName = "filter_featgroups_screening",
                   subDir = "featgroups", description = "Parameters for featureGroupsScreening filtering",
                   version = "1.0", definitions = getFilterFeatGroupsScreeningParamDefs(), ...)
})
