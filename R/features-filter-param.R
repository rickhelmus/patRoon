# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
NULL

#' @export
getFilterFeaturesParamDefs <- paramConfigDefsFact(list(
    absMinIntensity = list(
        advanced = FALSE,
        default = NULL,
        description = "Minimum absolute intensity threshold",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    relMinIntensity = list(
        advanced = FALSE,
        default = NULL,
        description = "Minimum relative intensity threshold",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    retentionRange = list(
        advanced = FALSE,
        default = NULL,
        description = "Retention time range (min, max)",
        type = "range",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    mzRange = list(
        advanced = FALSE,
        default = NULL,
        description = "m/z range (min, max)",
        type = "range",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    mzDefectRange = list(
        advanced = FALSE,
        default = NULL,
        description = "m/z defect range (min, max)",
        type = "range",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    chromWidthRange = list(
        advanced = FALSE,
        default = NULL,
        description = "Chromatographic width range (min, max)",
        type = "range",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    IMSRangeParams = list(
        advanced = FALSE,
        default = NULL,
        description = "Parameters for filtering by ion mobility range",
        type = "IMSRangeParams",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    qualityRange = list(
        advanced = FALSE,
        default = NULL,
        description = "Feature quality or score ranges",
        type = "scoreRange",
        # Score names require a features object and are validated when filter() is called.
        typeCheckArgs = list(scNames = NULL)
    ),
    sets = list(
        advanced = TRUE,
        default = NULL,
        description = "Sets to keep",
        type = "sets",
        typeCheckArgs = list(multiple = TRUE)
    ),
    negate = list(
        advanced = TRUE,
        default = FALSE,
        description = "Apply filters in the opposite manner",
        type = "flag"
    )
))

#' @export
FilterFeaturesParam <- setClass("FilterFeaturesParam", contains = "param")
setMethod("initialize", "FilterFeaturesParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "FilterFeaturesParam", baseName = "filter_features", subDir = "features",
                   description = "Parameters for features filtering", version = "1.0",
                   definitions = getFilterFeaturesParamDefs(), ...)
})
