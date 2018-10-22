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



convertOptToCallParamsEnviPick <- function(params)
{
    return(params)
}

checkInitialOptParamsEnviPick <- function(params, isoIdent)
{
    return(params)
}

# based on part of optimizeXcmsSet() function from IPO
fixOptParamBoundsEnviPick <- function(param, bounds)
{
    if (param %in% c("minpeak", "ended", "recurs"))
        return(round(bounds, 0))
    
    return(bounds)
}

# based on part of optimizeXcmsSet() function from IPO
fixOptParamsEnviPick <- function(params)
{
    return(fixOptParamRange(params, list(c("minint", "maxint"))))
}

# based on part of optimizeXcmsSet() function from IPO
getMinOptSettingEnviPick <- function(settingName, params)
{
    # UNDONE: need this?
    
    # if (settingName == "min_peakwidth")
    #     return(3)
    # if (settingName == "mzdiff")
    #     return(if (params$no_optimization$method == "centWave") -100000000 else 0.001)
    # if (settingName == "step")
    #     return(0.0005)
    
    return(1)
}
