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

checkInitialOptParamsXCMS <- function(params, isoIdent)
{
    return(params)
}

# based on part of optimizeXcmsSet() function from IPO
fixOptParamBoundsXCMS <- function(param, bounds)
{
    if (param == "steps" || param == "prefilter")
        return(round(bounds, 0))
    
    return(bounds)
}

fixOptParamsXCMS <- function(params)
{
    if (params$no_optimization$method == "centWave")
        return(fixOptParamRange(params, list(c("min_peakwidth", "max_peakwidth"))))
    return(params)
}

# based on part of optimizeXcmsSet() function from IPO
getMinOptSettingXCMS <- function(settingName, params)
{
    if (settingName == "min_peakwidth")
        return(3)
    if (settingName == "mzdiff")
        return(if (params$no_optimization$method == "centWave") -100000000 else 0.001)
    if (settingName == "step")
        return(0.0005)
    
    return(1)
}
