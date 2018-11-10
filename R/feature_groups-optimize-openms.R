#' @include main.R
#' @include feature_groups-optimize.R
NULL

featureGroupsOptimizerOpenMS <- setRefClass("featureGroupsOptimizerOpenMS", contains = "featureGroupsOptimizer")

featureGroupsOptimizerOpenMS$methods(

    defaultParamRanges = function(params)
    {
        return(list(maxAlignMZ = c(0.0001, Inf),
                    maxGroupMZ = c(0.0001, Inf)))
    }
)

generateFGroupsOptPSetOpenMS <- function(...)
{
    return(list(maxAlignRT = c(15, 45),
                maxAlignMZ = c(0.002, 0.010),
                maxGroupRT = c(6, 18),
                maxGroupMZ = c(0.002, 0.008)))
}
