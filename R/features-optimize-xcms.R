#' @include main.R
NULL

convertOptToCallParamsXCMS <- function(params)
{
    if (!is.null(params[["min_peakwidth"]])) # also implies max_peakwidth
    {
        params$peakwidth <- c(params$min_peakwidth, params$max_peakwidth)
        params[c("min_peakwidth", "max_peakwidth")] <- NULL
    }
    if (!is.null(params[["prefilter"]])) # also implies value_of_prefilter
    {
        params$prefilter <- c(params$prefilter, params$value_of_prefilter)
        params$value_of_prefilter <- NULL
    }
    return(params)
}

# based on part of optimizeXcmsSet() function from IPO
fixOptParamBoundsXCMS <- function(param, bounds)
{
    if (param == "steps" || param == "prefilter")
        return(round(bounds, 0))
    
    return(bounds)
}

# based on part of optimizeXcmsSet() function from IPO
checkOptParamsXCMS <- function(params)
{
    # UNDONE: move to algo specific functions
    
    if (params$no_optimization$method == "centWave")
    {
        #checking peakwidths plausiability
        if (!is.null(params$to_optimize$min_peakwidth) || 
            !is.null(params$to_optimize$max_peakwidth))
        {
            
            if (is.null(params$to_optimize$min_peakwidth))
                pw_min <- params$no_optimization$min_peakwidth
            else
                pw_min <- max(params$to_optimize$min_peakwidth)
            
            if (is.null(params$to_optimize$max_peakwidth))
                pw_max <- params$no_optimization$max_peakwidth
            else
                pw_max <- min(params$to_optimize$max_peakwidth)
            
            if (pw_min >= pw_max)
            {
                additional <- abs(pw_min-pw_max) + 1
                if (!is.null(params$to_optimize$max_peakwidth))
                    params$to_optimize$max_peakwidth <- params$to_optimize$max_peakwidth + additional
                else
                    params$no_optimization$max_peakwidth <- params$no_optimization$max_peakwidth + additional
            }
        }
    }
    
    return(params)
}

# based on part of optimizeXcmsSet() function from IPO
getMinOptSettingXCMS <- function(settingName, params)
{
    if (settingName == "min_peakwidth")
        return(3)
    if (settingName == "mzdiff")
        return(if (params$no_optimization$method == "centWave") 3 else 0.001)
    if (settingName == "step")
        return(0.0005)
    
    return(1)
}
