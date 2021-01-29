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

    fixDesignParam = function(param, value) if (param == "traceTermOutliers") round(value) else value,
    
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
    },
    
    calculateResponse = function(params, task, keepObject)
    {
        # disable (excessive) logging and parallelization if necessary
        opts <- list(patRoon.MP.logPath = FALSE)
        if (parallel)
            opts[c("patRoon.MP.method", "patRoon.MP.maxProcs")] <- list("classic", 1)
        withr::with_options(opts, callSuper(params, task, keepObject))
    }
)

generateFeatureOptPSetOpenMS <- function(...)
{
    return(list(chromFWHM = c(5, 10),
                mzPPM = c(3, 10),
                minFWHM = c(3, 6),
                maxFWHM = c(35, 65)))
}

getDefFeaturesOptParamRangesOpenMS <- function(params) list(localMZRange = c(0.00001, Inf),
                                                            traceTermOutliers = c(1, Inf))
