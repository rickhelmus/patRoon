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

    fixOptParams = function(params)
    {
        return(fixOptParamRange(params, list(c("minFWHM", "maxFWHM"),
                                             c("minTraceLength", "maxTraceLength"))))
    }
)

generateFeatureOptPSetOpenMS <- function(...)
{
    return(list(chromFWHM = c(5, 10),
                mzPPM = c(3, 10),
                minFWHM = c(3, 6),
                maxFWHM = c(35, 65)))
}

getDefFeaturesOptParamRangesOpenMS <- function(params) list(localMZRange = c(0.00001, Inf))
