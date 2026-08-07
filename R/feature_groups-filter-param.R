# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
NULL

#' @export
getFilterFeatGroupsParamDefs <- paramConfigDefsFact(list(
    preAbsMinIntensity = list(
        default = 100,
        description = "Pre-intensity filter: minimum absolute intensity threshold applied before other filters",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    absMinIntensity = list(
        default = 10000,
        description = "Minimum absolute intensity threshold",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    relMinReplicateAbundance = list(
        default = 1,
        description = "Minimum relative abundance within a replicate",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    maxReplicateIntRSD = list(
        default = 0.75,
        description = "Maximum relative standard deviation of intensities within a replicate",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    blankThreshold = list(
        default = 5,
        description = "Blank subtraction threshold (relative to blank intensity)",
        type = "number",
        typeCheckArgs = list(null.ok = TRUE, lower = 0, finite = TRUE)
    ),
    removeBlanks = list(
        default = TRUE,
        description = "Remove blank analyses after filtering",
        type = "flag"
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
    extraOpts = list(
        default = list(),
        description = "Extra filter options as a named list",
        type = "list",
        typeCheckArgs = list(null.ok = TRUE, any.missing = FALSE, names = "unique")
    )
))

#' @export
FilterFeatGroupsParam <- setClass("FilterFeatGroupsParam", contains = "param")
setMethod("initialize", "FilterFeatGroupsParam", function(.Object, ...)
{
    # split ... to args for constructor and extraOpts
    args <- list(...)
    defs <- getFilterFeatGroupsParamDefs()
    dotsConstr <- args[names(args) %in% names(defs)]
    dotsConstr$extraOpts <- args[setdiff(names(args), names(defs))]
    
    do.call(callNextMethod, c(list(.Object, name = "FilterFeatGroupsParam", baseName = "FilterFeatGroupsParam",
                                   description = "Parameters for featureGroups filtering", version = "1.0",
                                   definitions = defs), dotsConstr))
})
