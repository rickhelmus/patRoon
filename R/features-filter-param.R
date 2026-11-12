# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
NULL

#' @export
getFilterFeaturesParamDefs <- paramConfigDefsFact(list(
    absMinIntensity = list(
        default = NULL,
        description = "Minimum absolute intensity threshold",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    relMinIntensity = list(
        default = NULL,
        description = "Minimum relative intensity threshold",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    retentionRange = list(
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
        default = NULL,
        description = "m/z defect range (min, max)",
        type = "range",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    chromWidthRange = list(
        default = NULL,
        description = "Chromatographic width range (min, max)",
        type = "range",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    IMSRangeParams = list(
        default = NULL,
        description = "Parameters for filtering by ion mobility range",
        type = "IMSRangeParams",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    qualityRange = list(
        default = NULL,
        description = "Feature quality or score ranges"
    ),
    negate = list(
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

setValidity("FilterFeaturesParam", function(object)
{
    parsFilled <- paramListFillDefaults(object@data, object@definitions)
    ac <- checkmate::makeAssertCollection()
    # Score names require a features object and are validated when filter() is called.
    assertScoreRange(parsFilled$qualityRange, NULL, .var.name = "qualityRange", add = ac)
    tryCatch(checkmate::reportAssertions(ac), error = function(e) e)
})
