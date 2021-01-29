#' @include main.R
#' @include features-optimize.R
NULL

featuresOptimizerXCMS <- setRefClass("featuresOptimizerXCMS", contains = "featuresOptimizer")

featuresOptimizerXCMS$methods(

    checkInitialParams = function(params)
    {
        if (is.null(params[["method"]]))
            params$method <- "centWave"
        return(params)
    },

    # Adapted from combineParams() function of IPO
    combineParams = function(params_1, params_2)
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

    fixDesignParam = function(param, value) if (param %in% c("steps", "prefilter")) round(value) else value,
    
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
    },
    
    calculateResponse = function(params, task, keepObject)
    {
        if (parallel)
            params$BPPARAM <- BiocParallel::SerialParam()
        callSuper(params, task, keepObject)
    }
)

generateFeatureOptPSetXCMS <- function(...)
{
    givenArgs <- list(...)
    ret <- list()

    # NOTE: method will be added to returned pset by generateFeatureOptPSet() if
    # already specified.
    if (!is.null(givenArgs[["method"]]))
        method <- givenArgs[["method"]]
    else
        method <- ret$method <- "centWave"
    
    if (method == "centWave")
    {
        # CHANGED: tightened ranges bit for modern equipment (ie smaller peakwidths and mz ranges)
        ret <- c(ret, list(min_peakwidth = c(4, 12), 
                           max_peakwidth = c(35, 65), 
                           ppm = c(5, 15),
                           mzdiff = c(-0.001, 0.01)))
    }
    else if (method == "matchedFilter")
    {
        ret <- c(ret, list(fwhm = c(25, 35), 
                           snthresh = c(3, 17), 
                           step = c(0.05, 0.15), 
                           steps = c(1, 3)))
    }
    else
        stop("Only centWave and matchedFilter methods supported")
    
    return(ret)
}

getDefFeaturesOptParamRangesXCMS <- function(method)
{
    return(list(min_peakwidth = c(3, Inf),
                mzdiff = c(if (method == "centWave") -100000000 else 0.001, Inf),
                step = c(0.0005, Inf)))
}
