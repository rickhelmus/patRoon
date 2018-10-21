#' @include main.R
#' @include features-optimize.R
NULL

featuresOptimizerOpenMS <- setRefClass("featuresOptimizerOpenMS", contains = "featuresOptimizer")

featuresOptimizerOpenMS$methods(
    
    checkInitialParams = function(params)
    {
        if (isoIdent != "OpenMS")
        {
            warning("Isotopic detection during feature finding will be disabled (by setting localMZRange=0). Set isoIdent=\"OpenMS\" to avoid this (see ?optimizeFeatureFinding).")
            params[["localMZRange"]] <- 0
        }
        return(params)
    },
    
    fixOptParamBounds = function(param, bounds)
    {
        if (param == "traceTermOutliers")
            return(round(bounds, 0))
        
        return(bounds)
    },
    
    # based on part of optimizeXcmsSet() function from IPO
    fixOptParams = function(params)
    {
        return(fixOptParamRange(params, list(c("minFWHM", "maxFWHM"),
                                             c("minTraceLength", "maxTraceLength"))))
    }
    
)



convertOptToCallParamsOpenMS <- function(params)
{
    return(params)
}

checkInitialOptParamsOpenMS <- function(params, isoIdent)
{
    if (isoIdent != "OpenMS")
    {
        warning("Isotopic detection during feature finding will be disabled (by setting localMZRange=0). Set isoIdent=\"OpenMS\" to avoid this (see ?optimizeFeatureFinding).")
        params[["localMZRange"]] <- 0
    }
    return(params)
}

# based on part of optimizeXcmsSet() function from IPO
fixOptParamBoundsOpenMS <- function(param, bounds)
{
    if (param == "traceTermOutliers")
        return(round(bounds, 0))
    
    return(bounds)
}

# based on part of optimizeXcmsSet() function from IPO
fixOptParamsOpenMS <- function(params)
{
    return(fixOptParamRange(params, list(c("minFWHM", "maxFWHM"),
                                         c("minTraceLength", "maxTraceLength"))))
}

# based on part of optimizeXcmsSet() function from IPO
getMinOptSettingOpenMS <- function(settingName, params)
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
