#' @include main.R
#' @include features-optimize.R
NULL

featuresOptimizerEnviPick <- setRefClass("featuresOptimizerEnviPick", contains = "featuresOptimizer")

featuresOptimizerEnviPick$methods(

    checkInitialParams = function(params)
    {
        params[["parallel"]] <- !parallel
        return(params)
    },
    
    fixDesignParam = function(param, value) if (param %in% c("minpeak", "ended", "recurs")) round(value) else value,
    
    fixOptParamBounds = function(param, bounds)
    {
        if (param %in% c("minpeak", "ended", "recurs"))
            return(round(bounds, 0))

        return(bounds)
    },

    fixOptParams = function(params)
    {
        return(fixOptParamRange(params, list(c("minint", "maxint"))))
    }
)

generateFeatureOptPSetEnviPick <- function(...)
{
    return(list(dmzdens = c(3, 10),
                drtsmall = c(10, 30),
                minpeak = c(4, 10),
                drttotal = c(60, 150),
                minint = c(1E3, 1E4),
                ppm = TRUE))
}

getDefFeaturesOptParamRangesEnviPick <- function() list()
