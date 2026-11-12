# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
NULL

#' @export
getFilterFormulasParamDefs <- paramConfigDefsFact(list(
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

FilterFormulasParam <- setClass("FilterFormulasParam", contains = "param")
setMethod("initialize", "FilterFormulasParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "FilterFormulasParam", baseName = "filter_formulas",
                   subDir = "annotation", description = "Parameters for formula filtering", version = "1.0",
                   definitions = getFilterFormulasParamDefs(), ...)
})
