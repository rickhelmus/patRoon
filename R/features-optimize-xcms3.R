#' @include main.R
#' @include features-optimize-xcms.R
NULL


featuresOptimizerXCMS3 <- setRefClass("featuresOptimizerXCMS3", contains = "featuresOptimizerXCMS")

featuresOptimizerXCMS3$methods(
    
    convertOptToCallParams = function(params)
    {
        params <- callSuper(params)
        
        method <- params$method
        params[["method"]] <- NULL
        if (method == "centWave")
            return(list(param = do.call(xcms::CentWaveParam, params)))
        # if (method == "matchedFilter")
        return(list(param = do.call(xcms::MatchedFilterParam, params)))
    },
    
    calculateResponse = function(params, task, keepObject)
    {
        if (parallel)
            params$BPPARAM <- BiocParallel::SerialParam()
        callSuper(params, task, keepObject)
    }
)

generateFeatureOptPSetXCMS3 <- function(...)
{
    ret <- generateFeatureOptPSetXCMS(...)
    if (!is.null(ret[["step"]]))
        names(ret)[which(names(ret) == "step")] <- "binSize"
    return(ret)
}
getDefFeaturesOptParamRangesXCMS3 <- getDefFeaturesOptParamRangesXCMS
