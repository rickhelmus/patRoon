#' @include utils-IPO.R
#' @include main.R
NULL

callAlgoFunc <- function(func, algorithm, ...)
{
    func <- substitute(func)
    suffix <- switch(algorithm,
                     xcms = "XCMS",
                     openms = "OpenMS",
                     envipick = "enviPick",
                     stop("Invalid algorithm!"))
    do.call(paste0(func, suffix), list(...))
}

# Adapted from combineParams() function of IPO
combineOptParams = function(params_1, params_2, algorithm)
{
    len <- max(unlist(sapply(params_1, length)))
    #num_params <- length(params_1)
    
    p_names <- c(names(params_1), names(params_2))
    matchedFilter <- algorithm == "xcms" &&
        ((!is.null(params_1[["method"]]) && params_1[["method"]] == "matchedFilter") ||
         (!is.null(params_2[["method"]]) && params_2[["method"]] == "matchedFilter"))
    
    for(i in 1:length(params_2))
    {
        new_index <- length(params_1) + 1
        fact <- params_2[[i]]
        params_1[[new_index]] <- fact
        if (matchedFilter)
        {
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
        else
        {
            # standard: replicate value
            params_1[[new_index]][1:len] <- fact
        }
    } 
    names(params_1) <- p_names   
    return(params_1)
}


calculateFeatures <- function(anaInfo, algorithm, params, task)
{
    # UNDONE: do we want to keep caching this?

    # params is a data.frame when called from performOptimIteration() and a list
    # when called from performOptimIterationStat().
    if (is.data.frame(params))
        params <- as.list(params[task, ])
    else
        params <- as.list(params[task])
        
    # get rid of design specific params
    params[["run.order"]] <- NULL
    params[["std.order"]] <- NULL
    params[["Block"]] <- NULL
    
    params <- callAlgoFunc(convertOptToCallParams, algorithm, params)
    return(do.call(findFeatures, c(list(anaInfo, algorithm), params)))
}

# Heavily based on xcmsSetExperimentsCluster() from IPO
performOptimIteration <- function(anaInfo, algorithm, params, isoIdent)
{
    typParams <- utilsIPO$typeCastParams(params) 
    
    if (length(typParams$to_optimize) > 1)
    {
        design <- utilsIPO$getCcdParameter(typParams$to_optimize)  	
        designParams <- rsm::decode.data(design) 
    }
    else
    {
        design <- data.frame(run.order = 1:9, a = seq(-1,1,0.25))
        colnames(design)[2] <- names(typParams$to_optimize)
        designParams <- design
        designParams[,2] <- seq(min(typParams$to_optimize[[1]]),
                                max(typParams$to_optimize[[1]]), 
                                diff(typParams$to_optimize[[1]]) / 8)
    }
    
    designParams <- combineOptParams(designParams, typParams$no_optimization, algorithm)   
    tasks <- seq_len(nrow(design))
    
    response <- sapply(tasks, function(task)
    {
        # simplified optimizeSlaveCluster() from IPO
        feat <- calculateFeatures(anaInfo, algorithm, designParams, task)
        result <- utilsIPO$calcPPS(getXcmsSet(feat, TRUE), isoIdent)
        result[1] <- task
        return(result)
    })
    
    response <- t(response)
    colnames(response) <- c("exp", "num_peaks", "notLLOQP", "num_C13", "PPS")
    response <- response[order(response[,1]),]
    
    ret <- list()
    ret$params <- typParams
    ret$design <- design
    ret$response <- response
    
    return(ret)
}

# Heavily based on xcmsSetStatistic() from IPO
performOptimIterationStat <- function(anaInfo, algorithm, result, isoIdent)
{
    params <- result$params
    resp <- result$response[, "PPS"]
    
    result$model <- utilsIPO$createModel(result$design, params$to_optimize, resp)
    maxSettings <- utilsIPO$getMaximumLevels(result$model)

    # plotting rms --> UNDONE
    # tmp <- max_settings[1,-1] # first row without response
    # tmp[is.na(tmp)] <- 1 # if Na (i.e. -1, and 1), show +1
    # if (isTRUE(plot) & length(tmp) > 1) {
    #     if (!is.null(subdir)) {
    #         plotContours(
    #             model = result$model,
    #             maximum_slice = tmp,
    #             plot_name = file.path(subdir, paste("rsm_", iterator, sep = ""))
    #         )
    #     } else if (is.null(subdir))  {
    #         plotContours(result$model, tmp, plot_name = NULL)
    #     }
    # }
    
    result$max_settings <- maxSettings
    
    runParams <- as.list(utilsIPO$decodeAll(maxSettings[-1], params$to_optimize))      
    runParams <- combineOptParams(runParams, params$no_optimization, algorithm)
    
    if (!is.list(runParams))
        runParams <- as.list(runParams)
    
    result$features <- calculateFeatures(anaInfo, algorithm, runParams, 1)
    result$PPS <- utilsIPO$calcPPS(getXcmsSet(result$features, TRUE), isoIdent) # UNDONE: what about the "..." params?
    
    return(result)
}

fixOptParamBounds <- function(param, bounds, algorithm) callAlgoFunc(fixOptParamBounds, algorithm, param, bounds)
checkOptParams <- function(params, algorithm) callAlgoFunc(checkOptParams, algorithm, params)
getMinOptSetting <- function(settingName, algorithm, params) callAlgoFunc(getMinOptSetting, algorithm, settingName, params)

# heavily based on optimizeXcmsSet() from IPO
optimizeFeatureFinding <- function(anaInfo, algorithm, params, isoIdent = "IPO",
                                   maxIterations = 50)
{
    ac <- checkmate::makeAssertCollection()
    assertAnalysisInfo(anaInfo, add = ac)
    checkmate::assertChoice(algorithm, c("openms", "xcms", "enviPick"), add = ac)
    checkmate::assertList(params, add = ac)
    checkmate::assertChoice(isoIdent, c("IPO", "CAMERA"), add = ac)
    checkmate::assertInt(maxIterations, add = ac)
    checkmate::reportAssertions(ac)
    
    history <- list()
    bestRange <- 0.25
    
    for (iter in seq_len(maxIterations))
    {
        result <- performOptimIteration(anaInfo, algorithm, params, isoIdent)
        result <- performOptimIterationStat(anaInfo, algorithm, result, isoIdent)
        
        history[[iter]] <- result
        params <- result$params
        
        if (!utilsIPO$resultIncreased(history))
        {
            maxima <- 0
            maxIndex <- 1
            for (i in seq_len(length(history)))
            {
                if (history[[i]]$max_settings[1] > maxima)
                {
                    maxima <- history[[i]]$max_settings[1]
                    maxIndex <- i
                }
            }
            
            finalParams <- as.list(utilsIPO$decodeAll(history[[maxIndex]]$max_settings[-1],
                                                      history[[maxIndex]]$params$to_optimize))      
            finalParams <- combineOptParams(finalParams, params$no_optimization, algorithm)
            
            if (!is.list(finalParams))
                finalParams <- as.list(finalParams)
            
            bestSettings <- list()
            bestSettings$parameters <- finalParams
            
            bestSettings$features <- history[[maxIndex]]$features
            bestSettings$result <- history[[maxIndex]]$PPS
            history$bestSettings <- bestSettings
            
            return(history)   
        }
        
        for(i in seq_len(length(params$to_optimize)))
        {
            setting <- result$max_settings[i+1]
            bounds <- params$to_optimize[[i]] 
            settingName <- names(params$to_optimize)[i]
            
            minSetting <- getMinOptSetting(settingName, algorithm, params)

            # - if the parameter is NA, we increase the range by 20%, 
            # - if it was within the inner 25% of the previous range or
            #   at the minimum value we decrease the range by 20%
            if (is.na(setting))
                stepFactor <- 1.2
            else if (abs(setting) < bestRange ||
                     (setting == -1 && utilsIPO$decode(-1, params$to_optimize[[i]]) == minSetting))
                stepFactor <- 0.8
            else
                stepFactor <- 1
            
            step <- (diff(bounds) / 2) * stepFactor
            
            if (is.na(setting))
                setting <- 0
            
            newCenter <- utilsIPO$decode(setting, bounds)
            
            if ((newCenter-minSetting) > step)
                newBounds <- c(newCenter - step, newCenter + step) 
            else
                newBounds <- c(minSetting, 2*step+minSetting) 
            
            names(newBounds) <- NULL
            
            params$to_optimize[[i]] <- fixOptParamBounds(names(params$to_optimize)[i], newBounds, algorithm)
        }
        
        params <- checkOptParams(params, algorithm)
        params <- utilsIPO$attachList(params$to_optimize, params$no_optimization)
    }
    
    params <- utilsIPO$attachList(params$to_optimize, params$no_optimization)
    
    return(history)
}
