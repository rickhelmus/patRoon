#' @include main.R
NULL

convertOptToCallParamsOpenMS <- function(params)
{
    return(params)
}

# based on part of optimizeXcmsSet() function from IPO
fixOptParamBoundsOpenMS <- function(param, bounds)
{
    if (param == "trace_termination_outliers")
        return(round(bounds, 0))
    
    return(bounds)
}

# based on part of optimizeXcmsSet() function from IPO
checkOptParamsOpenMS <- function(params)
{
    # UNDONE: need this?
    
    # if (params$no_optimization$method == "centWave")
    # {
    #     #checking peakwidths plausiability
    #     if (!is.null(params$to_optimize$min_peakwidth) || 
    #         !is.null(params$to_optimize$max_peakwidth))
    #     {
    #         
    #         if (is.null(params$to_optimize$min_peakwidth))
    #             pw_min <- params$no_optimization$min_peakwidth
    #         else
    #             pw_min <- max(params$to_optimize$min_peakwidth)
    #         
    #         if (is.null(params$to_optimize$max_peakwidth))
    #             pw_max <- params$no_optimization$max_peakwidth
    #         else
    #             pw_max <- min(params$to_optimize$max_peakwidth)
    #         
    #         if (pw_min >= pw_max)
    #         {
    #             additional <- abs(pw_min-pw_max) + 1
    #             if (!is.null(params$to_optimize$max_peakwidth))
    #                 params$to_optimize$max_peakwidth <- params$to_optimize$max_peakwidth + additional
    #             else
    #                 params$no_optimization$max_peakwidth <- params$no_optimization$max_peakwidth + additional
    #         }
    #     }
    # }
    
    return(params)
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
