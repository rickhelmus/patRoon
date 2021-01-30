#' @include main.R
#' @include features-optimize.R
NULL

featuresOptimizerKPIC2 <- setRefClass("featuresOptimizerKPIC2", contains = "featuresOptimizer")

featuresOptimizerKPIC2$methods(
    
    checkInitialParams = function(params)
    {
        params[["parallel"]] <- !parallel
        return(params)
    },
    
    fixDesignParam = function(param, value) if (param == "gap") round(value) else value,
    
    fixOptParamBounds = function(param, bounds)
    {
        if (param == "gap")
            return(round(bounds, 0))
        return(bounds)
    },
    
    fixOptParams = function(params)
    {
        return(fixOptParamRange(params, list(c("min_width", "max_width"))))
    },
    
    convertOptToCallParams = function(params)
    {
        if (!is.null(params[["min_width"]])) # also implies max_width
        {
            params$width <- c(params$min_width, params$peakwidth)
            params[c("min_width", "max_width")] <- NULL
        }
        return(params)
    }
)

generateFeatureOptPSetKPIC2 <- function(...)
{
    return(list(level = c(1E3, 1E4),
                mztol = c(0.001, 0.01),
                min_width = c(3, 8),
                max_width = c(35, 65),
                min_snr = c(2, 8),
                kmeans = TRUE))
}

getDefFeaturesOptParamRangesKPIC2 <- function() list(min_width = 3)
