# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
NULL

getFilterTPsParamDefs <- paramConfigDefsFact(list(
    properties = list(
        default = NULL,
        description = "Allowed property values",
        # NOTE: properties requires additional checks with data from the TPs object
        type = "list",
        typeCheckArgs = list(null.ok = TRUE, any.missing = FALSE)
    ),
    verbose = list(
        default = TRUE,
        description = "Verbose output",
        type = "flag"
    ),
    negate = list(
        default = FALSE,
        description = "Apply filters in the opposite manner",
        type = "flag"
    )
))

FilterTPsParam <- setClass("FilterTPsParam", contains = "param")
setMethod("initialize", "FilterTPsParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "FilterTPsParam", baseName = "filter_tps", subDir = "tps",
                   description = "Parameters for transformation product filtering", version = "1.0",
                   definitions = getFilterTPsParamDefs(), ...)
})


getFilterTPsStructureParamDefs <- paramConfigDefsFact(list(
    removeParentIsomers = list(
        default = FALSE,
        description = "Remove TPs isomeric with their parent",
        type = "flag"
    ),
    removeTPIsomers = list(
        default = FALSE,
        description = "Remove isomeric TPs",
        type = "flag"
    ),
    removeDuplicates = list(
        default = FALSE,
        description = "Remove duplicate structures",
        type = "flag"
    ),
    minSimilarity = list(
        default = NULL,
        description = "Minimum structural similarity",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE, null.ok = TRUE)
    ),
    verbose = list(
        default = TRUE,
        description = "Verbose output",
        type = "flag"
    ),
    negate = list(
        default = FALSE,
        description = "Apply filters in the opposite manner",
        type = "flag"
    )
))

FilterTPsStructureParam <- setClass("FilterTPsStructureParam", contains = "param")
setMethod("initialize", "FilterTPsStructureParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "FilterTPsStructureParam", baseName = "filter_tps_structure",
                   subDir = "tps",
                   description = "Parameters for TP filtering with structure data", version = "1.0",
                   definitions = getFilterTPsStructureParamDefs(), ...)
})


getFilterTPsAnnCompParamDefs <- paramConfigDefsFact(list(
    minFitFormula = list(
        default = 0,
        description = "Minimum formula fit",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    minFitCompound = list(
        default = 0,
        description = "Minimum compound fit",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    minSimSusp = list(
        default = 0,
        description = "Minimum suspect similarity",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    minFitCompOrSimSusp = list(
        default = c(0, 0),
        description = "Minimum compound fit OR suspect similarity",
        type = "numeric",
        typeCheckArgs = list(lower = 0, finite = TRUE, len = 2)
    ),
    minTPScore = list(
        default = 0,
        description = "Minimum TP score",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    topMost = list(
        default = NULL,
        description = "Maximum number of TPs to keep per group",
        type = "count",
        typeCheckArgs = list(positive = TRUE, null.ok = TRUE)
    ),
    verbose = list(
        default = TRUE,
        description = "Verbose output",
        type = "flag"
    ),
    negate = list(
        default = FALSE,
        description = "Apply filters in the opposite manner",
        type = "flag"
    )
))


FilterTPsAnnCompParam <- setClass("FilterTPsAnnCompParam", contains = "param")
setMethod("initialize", "FilterTPsAnnCompParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "FilterTPsAnnCompParam", baseName = "filter_tps_ann_comp",
                   subDir = "tps",
                   description = "Parameters for filtering TPs from compound annotations", version = "1.0",
                   definitions = getFilterTPsAnnCompParamDefs(), ...)
})

getFilterTPsAnnFormParamDefs <- paramConfigDefsFact(list(
    minFitFormula = list(
        default = 0,
        description = "Minimum formula fit",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    minTPScore = list(
        default = 0,
        description = "Minimum TP score",
        type = "number",
        typeCheckArgs = list(lower = 0, finite = TRUE)
    ),
    topMost = list(
        default = NULL,
        description = "Maximum number of TPs to keep per group",
        type = "count",
        typeCheckArgs = list(positive = TRUE, null.ok = TRUE)
    ),
    verbose = list(
        default = TRUE,
        description = "Verbose output",
        type = "flag"
    ),
    negate = list(
        default = FALSE,
        description = "Apply filters in the opposite manner",
        type = "flag"
    )
))

FilterTPsAnnFormParam <- setClass("FilterTPsAnnFormParam", contains = "param")
setMethod("initialize", "FilterTPsAnnFormParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "FilterTPsAnnFormParam", baseName = "filter_tps_ann_form",
                   subDir = "tps",
                   description = "Parameters for filtering TPs from formula annotations", version = "1.0",
                   definitions = getFilterTPsAnnFormParamDefs(), ...)
})
