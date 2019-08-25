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
    }
)

generateFeatureOptPSetXCMS3 <- generateFeatureOptPSetXCMS
getDefFeaturesOptParamRangesXCMS3 <- getDefFeaturesOptParamRangesXCMS
