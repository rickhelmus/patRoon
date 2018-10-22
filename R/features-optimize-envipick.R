#' @include main.R
#' @include features-optimize.R
NULL

featuresOptimizerEnviPick <- setRefClass("featuresOptimizerEnviPick", contains = "featuresOptimizer")

featuresOptimizerEnviPick$methods(
    
    fixOptParamBounds = function(param, bounds)
    {
        if (param %in% c("minpeak", "ended", "recurs"))
            return(round(bounds, 0))
        
        return(bounds)
    },
    
    fixOptParams = function(params)
    {
        return(fixOptParamRange(params, list(c("minint", "maxint"))))
    },
    
    getMinOptSetting = function(settingName, params)
    {
        return(1)
    }
)
