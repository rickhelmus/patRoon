# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
NULL

#' @export
getFilterCompoundsParamDefs <- paramConfigDefsFact(list(
    minScore = list(
        default = NULL,
        description = "Minimum overall compound score",
        type = "number",
        typeCheckArgs = list(finite = TRUE, null.ok = TRUE)
    ),
    minFragScore = list(
        default = NULL,
        description = "Minimum fragmentation score",
        type = "number",
        typeCheckArgs = list(finite = TRUE, null.ok = TRUE)
    ),
    minFormulaScore = list(
        default = NULL,
        description = "Minimum formula score",
        type = "number",
        typeCheckArgs = list(finite = TRUE, null.ok = TRUE)
    ),
    IMSRangeParams = list(
        default = NULL,
        description = "Ion mobility range filter",
        type = "IMSRangeParams",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    IMSMatchParams = list(
        default = NULL,
        description = "Ion mobility match filter",
        type = "IMSMatchParams",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    negate = list(
        default = FALSE,
        description = "Apply filters in the opposite manner",
        type = "flag"
    ),
    sets = list(
        default = NULL,
        description = "Sets to keep",
        type = "sets",
        typeCheckArgs = list(multiple = TRUE)
    ),
    updateConsensus = list(
        default = FALSE,
        description = "Update set consensus after filtering",
        type = "flag"
    )
))

FilterCompoundsParam <- setClass("FilterCompoundsParam", contains = "param")
setMethod("initialize", "FilterCompoundsParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "FilterCompoundsParam", baseName = "filter_compounds",
                   subDir = "annotation", description = "Parameters for compound filtering", version = "1.0",
                   definitions = getFilterCompoundsParamDefs(), ...)
})
