#' @include utils-IPO.R
#' @include main.R
NULL

DoEOptimizer <- setRefClass("DoEOptimizer", contains = "VIRTUAL",
                            fields = list(anaInfo = "data.frame", algorithm = "character"))

DoEOptimizer$methods(
    
    # dummy methods that may need to be overloaded
    checkInitialParams = function(params) params,
    getMinOptSetting = function(settingName, params) 1,
    fixParamBounds = function(params, bounds) bounds,
    fixOptParams = function(params) params,
    
    # "virtual" methods
    resultIncreased = function(history) stop("VIRTUAL"),
    calculateResponse = function(params, task, final = FALSE) stop("VIRTUAL"),
    
    # Adapted from combineParams() function of IPO
    combineParams = function(params_1, params_2)
    {
        len <- max(unlist(sapply(params_1, length)))

        p_names <- c(names(params_1), names(params_2))
        for (i in seq_along(params_2))
        {
            new_index <- length(params_1) + 1
            fact <- params_2[[i]]
            params_1[[new_index]] <- fact
            params_1[[new_index]][1:len] <- fact
        } 
        names(params_1) <- p_names   
        return(params_1)
    },
    
    # Heavily based on xcmsSetExperimentsCluster() from IPO
    performIteration = function(params)
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
        
        designParams <- combineParams(designParams, typParams$no_optimization)   
        tasks <- seq_len(nrow(design))
        
        response <- rbindlist(lapply(tasks, function(task)
        {
            # simplified optimizeSlaveCluster() from IPO
            result <- calculateResponse(designParams, task)
            result$experiment <- task
            return(result)
        }))
        
        ret <- list()
        ret$params <- typParams
        ret$design <- design
        ret$response <- response
        
        return(ret)
    },

    # Heavily based on xcmsSetStatistic() from IPO
    performIterationStat = function(result)
    {
        params <- result$params
        resp <- result$response$score
        
        result$model <- utilsIPO$createModel(result$design, params$to_optimize, resp)
        result$max_settings <- utilsIPO$getMaximumLevels(result$model)
        
        runParams <- as.list(utilsIPO$decodeAll(result$max_settings[-1], params$to_optimize))      
        runParams <- combineParams(runParams, params$no_optimization)
        
        if (!is.list(runParams))
            runParams <- as.list(runParams)
        
        result$finalResponse <- calculateResponse(runParams, task, final = TRUE)
        # UNDONE: what about the "..." params?
        
        return(result)
    },
    
    # heavily based on optimizeXcmsSet() from IPO
    optimize = function(params, maxIterations)
    {
        params <- startParams <- checkInitialParams(params)
        
        history <- list()
        bestRange <- 0.25
        
        for (iter in seq_len(maxIterations))
        {
            printf("\n\n===\n")
            printf("Starting new DoE with:\n")
            printf(paste0(rbind(paste0(names(params), ": "), 
                                paste0(params, "\n"))))
            printf("===\n\n")
            
            result <- performIteration(params)
            result <- performIterationStat(result)
            
            history[[iter]] <- result
            params <- result$params
            
            if (!resultIncreased(history))
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
                finalParams <- combineParams(finalParams, params$no_optimization)
                
                if (!is.list(finalParams))
                    finalParams <- as.list(finalParams)
                
                bestSettings <- list()
                bestSettings$parameters <- finalParams
                
                bestSettings$object <- history[[maxIndex]]$finalResponse$object
                bestSettings$result <- history[[maxIndex]]$finalResponse[names(history[[maxIndex]]$finalResponse) != "object"]
                bestSettings$score <- history[[maxIndex]]$finalResponse$score
                history$bestSettings <- bestSettings
                
                break  
            }
            
            for (i in seq_len(length(params$to_optimize)))
            {
                setting <- result$max_settings[i+1]
                bounds <- params$to_optimize[[i]] 
                settingName <- names(params$to_optimize)[i]
                
                minSetting <- getMinOptSetting(settingName, params)
                
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
                
                params$to_optimize[[i]] <- fixParamBounds(names(params$to_optimize)[i], newBounds)
            }
            
            params <- fixOptParams(params)
            params <- utilsIPO$attachList(params$to_optimize, params$no_optimization)
        }
        
        #params <- utilsIPO$attachList(params$to_optimize, params$no_optimization)
        
        return(list(startParams = startParams, finalResults = history$bestSettings, experiments = history[seq_len(iter)]))
    }
)

# based on part of optimizeXcmsSet() function from IPO
fixOptParamRange <- function(params, paramPairs)
{
    for (pp in paramPairs)
    {
        if (!is.null(params$to_optimize[[pp[1]]]) || !is.null(params$to_optimize[[pp[2]]]))
        {
            if (is.null(params$to_optimize[[pp[1]]]))
                pmin <- params$no_optimization[[pp[1]]]
            else
                pmin <- max(params$to_optimize[[pp[1]]])
            
            if (is.null(params$to_optimize[[pp[2]]]))
                pmax <- params$no_optimization[[pp[2]]]
            else
                pmax <- min(params$to_optimize[[pp[2]]])
            
            if (pmin >= pmax)
            {
                additional <- abs(pmin-pmax) + 1
                if (!is.null(params$to_optimize[[pp[2]]]))
                    params$to_optimize[[pp[2]]] <- params$to_optimize[[pp[2]]] + additional
                else
                    params$no_optimization[[pp[2]]] <- params$no_optimization[[pp[2]]] + additional
            }
        }
    }
    
    return(params)
}
