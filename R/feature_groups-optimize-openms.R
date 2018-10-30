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
