#' @include main.R
#' @include features-optimize.R
NULL

featuresOptimizerXCMS <- setRefClass("featuresOptimizerXCMS", contains = "featuresOptimizer")

featuresOptimizerXCMS$methods(
    
    # Adapted from combineParams() function of IPO
    combineOptParams = function(params_1, params_2)
    {
        len <- max(unlist(sapply(params_1, length)))
        #num_params <- length(params_1)
        
        p_names <- c(names(params_1), names(params_2))
        matchedFilter <- ((!is.null(params_1[["method"]]) && params_1[["method"]] == "matchedFilter") ||
                          (!is.null(params_2[["method"]]) && params_2[["method"]] == "matchedFilter"))
        
        if (!matchedFilter)
            return(callSuper(params_1, params_2))
        
        for (i in seq_along(params_2))
        {
            new_index <- length(params_1) + 1
            fact <- params_2[[i]]
            params_1[[new_index]] <- fact

            if (p_names[new_index] == "sigma" && fact == 0)
            {
                # update values for sigma if zero
                if ("fwhm" %in% names(params_1))
                    params_1[[new_index]][1:len] <- params_1$fwhm/2.3548
                else
                    params_1[[new_index]][1:len] <- params_2$fwhm/2.3548
            }
            else if (p_names[new_index] == "mzdiff" && fact == 0)
            {
                # update values for mzdiff if zero
                if ("step" %in% names(params_1))
                {
                    if ("steps"  %in% names(params_1))
                        params_1[[new_index]][1:len] <- 0.8-params_1$step*params_1$steps
                    else
                        params_1[[new_index]][1:len] <- 0.8-params_1$step*params_2$steps
                }
                else
                {
                    if ("steps"  %in% names(params_1))
                        params_1[[new_index]][1:len] <- 0.8-params_2$step*params_1$steps
                    else
                        params_1[[new_index]][1:len] <- 0.8-params_2$step*params_2$steps
                }
            }
            else
            {  
                # standard: replicate value
                params_1[[new_index]][1:len] <- fact
            }
        } 
        names(params_1) <- p_names   
        return(params_1)
    },
    
    fixOptParamBounds = function(param, bounds)
    {
        if (param == "steps" || param == "prefilter")
            return(round(bounds, 0))
        
        return(bounds)
    },
    
    # based on part of optimizeXcmsSet() function from IPO
    fixOptParams = function(params)
    {
        return(fixOptParamRange(params, list(c("min_peakwidth", "max_peakwidth"))))
    },
    
    getMinOptSetting = function(settingName, params)
    {
        if (settingName == "min_peakwidth")
            return(3)
        if (settingName == "mzdiff")
            return(if (params$no_optimization$method == "centWave") -100000000 else 0.001)
        if (settingName == "step")
            return(0.0005)
        return(1)
    },
    
    convertOptToCallParams = function(params)
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
)
