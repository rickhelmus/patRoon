# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
#' @include feature_groups-optimize.R
NULL

featureGroupsOptimizerOpenMS <- setRefClass("featureGroupsOptimizerOpenMS", contains = "featureGroupsOptimizer")

generateFGroupsOptPSetOpenMS <- function(...)
{
    return(list(maxAlignRT = c(15, 45),
                maxAlignMZ = c(0.002, 0.010),
                maxGroupRT = c(6, 18),
                maxGroupMZ = c(0.002, 0.008)))
}

getDefFGroupsOptParamRangesOpenMS <- function()
{
    return(list(maxAlignMZ = c(0.0001, Inf),
                maxGroupMZ = c(0.0001, Inf)))
}