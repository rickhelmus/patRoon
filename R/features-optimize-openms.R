#' @include main.R
NULL

convertOptToCallParamsOpenMS <- function(params)
{
    return(params)
}

checkInitialOptParamsOpenMS <- function(params)
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
fixOptParamsOpenMS <- function(params)
{
    return(fixOptParamRange(params, list(c("minfwhm", "maxfwhm"),
                                         c("minlength", "maxlength"))))
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
