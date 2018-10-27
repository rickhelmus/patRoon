#' @include main.R
#' @include feature_groups-optimize.R
NULL

featureGroupsOptimizerOpenMS <- setRefClass("featureGroupsOptimizerOpenMS", contains = "featureGroupsOptimizer")

featureGroupsOptimizerOpenMS$methods(
    
    getMinOptSetting = function(settingName, params)
    {
        if (settingName %in% c("maxAlignMZ", "maxGroupMZ"))
            return(0.0001)

        return(1) # not reached
    }
    
)