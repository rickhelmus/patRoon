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

# Adapted from IPO: add OpenMS isotope detection
calcPPS <- function(feat, isotopeIdentification = c("IPO", "CAMERA", "OpenMS"), ...)
{
    fTable <- featureTable(feat)
    
    isotopeIdentification <- match.arg(isotopeIdentification)
    
    ret <- vector(mode="numeric", 5) #array(0, dim=c(1,5)) 
    names(ret) <- c("ExpId", "#peaks", "#NonRP", "#RP", "PPS")
    
    if (length(feat) == 0)
        return(ret)

    ret[2] <- length(feat)
    
    doOpenMS <- isotopeIdentification == "OpenMS"
    if (!doOpenMS) # no need to find isotopes with OpenMS algo
    {
        xset <- getXcmsSet(feat, TRUE)
        peak_source <- utilsIPO$peaks_IPO(xset)[, c("mz", "rt", "sample", "into", "mzmin", 
                                                    "mzmax", "rtmin", "rtmax"), drop=FALSE]
        if(isotopeIdentification == "IPO")
            iso_mat <- utilsIPO$findIsotopes.IPO(xset, ...)  
        else
            iso_mat <- utilsIPO$findIsotopes.CAMERA(xset, ...)
    }
    
    isotope_abundance = 0.01108    
    
    #calculating low intensity peaks
    for (anai in seq_along(analyses(feat)))
    {
        if (doOpenMS)
        {
            # intensities/masses of features without assigned isotopes
            intensities <- fTable[[anai]][isocount == 1, intensity]
            masses <- fTable[[anai]][isocount == 1, mz]
        }
        else
        {
            non_isos_peaks <- peak_source
            
            if (nrow(iso_mat) > 0)
                non_isos_peaks <- peak_source[-unique(c(iso_mat)),,drop=FALSE] 
            
            speaks <- non_isos_peaks[non_isos_peaks[,"sample"]==anai,,drop=FALSE]
            intensities <- speaks[,"into"]
            na_int <- is.na(intensities)
            intensities <- intensities[!na_int]
            
            if (length(intensities) > 0)
            {
                masses <- speaks[!na_int, "mz"]
                #floor((masses-2*CH3)/CH2) + 2
            }            
        }
        
        if (length(intensities) > 0)
        {
            tmp <- intensities[order(intensities)]
            int_cutoff <- mean(tmp[1:max(round((length(tmp)/33),0),1)])
            
            maximum_carbon <- utilsIPO$calcMaximumCarbon(masses)
            carbon_probability <- maximum_carbon * isotope_abundance
            
            iso_int <- intensities * carbon_probability
            
            not_loq_peaks <- sum(iso_int > int_cutoff)
            ret[3] <- ret[3] + not_loq_peaks
        }
    }#end_for_sample    
    
    if (doOpenMS)
    {
        # isocount represent the number of isotopes collapsed in a feature. When
        # it's one, no other isotopes are present and hence its a NP. The
        # remaining are assumed to be RP.
        ret[4] <- sum(unlist(lapply(fTable, function(ft) ft[isocount > 1, isocount])))
    }
    else
        ret[4] <- length(unique(c(iso_mat)))
    
    if (ret[3] == 0)
        ret[5] <- (ret[4]+1)^2/(ret[3]+1)  
    else
        ret[5] <- ret[4]^2/ret[3]  
    
    return(ret)
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
    
    for (i in seq_along(params_2))
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

calculateFeatures <- function(anaInfo, algorithm, params, task, final = FALSE)
{
    # UNDONE: do we want to keep caching this?

    # params is a data.frame when not final and a list when it is final.
    if (!final)
        params <- as.list(params[task, ])

    # get rid of design specific params
    params[["run.order"]] <- NULL
    params[["std.order"]] <- NULL
    params[["Block"]] <- NULL

    if (!final)
        printf("---\nTask %d\n", task)
    else
        printf("---\nGetting features with final settings...\n")
    printf("%s: %s\n", names(params), params)
    printf("---\n")
    
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
        result <- calcPPS(feat, isoIdent)
        result[1] <- task
        return(result)
    })
    
    response <- t(response)
    colnames(response) <- c("exp", "num_peaks", "notLLOQP", "num_C13", "PPS")
    response <- response[order(response[, 1]), ]
    
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
    
    result$features <- calculateFeatures(anaInfo, algorithm, runParams, final = TRUE)
    result$PPS <- calcPPS(result$features, isoIdent) # UNDONE: what about the "..." params?
    
    return(result)
}

checkInitialOptParams <- function(params, algorithm) callAlgoFunc(checkInitialOptParams, algorithm, params)
fixOptParamBounds <- function(param, bounds, algorithm) callAlgoFunc(fixOptParamBounds, algorithm, param, bounds)
fixOptParams <- function(params, algorithm) callAlgoFunc(fixOptParams, algorithm, params)
getMinOptSetting <- function(settingName, algorithm, params) callAlgoFunc(getMinOptSetting, algorithm, settingName, params)

# heavily based on optimizeXcmsSet() from IPO
optimizeFeatureFinding <- function(anaInfo, algorithm, params, isoIdent = "IPO", maxIterations = 50)
{
    ac <- checkmate::makeAssertCollection()
    assertAnalysisInfo(anaInfo, add = ac)
    checkmate::assertChoice(algorithm, c("openms", "xcms", "envipick"), add = ac)
    checkmate::assertList(params, add = ac)
    checkmate::assertChoice(isoIdent, c("IPO", "CAMERA", "OpenMS"), add = ac)
    checkmate::assertCount(maxIterations, positive = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    params <- checkInitialOptParams(params, algorithm)
    
    history <- list()
    bestRange <- 0.25
    
    for (iter in seq_len(maxIterations))
    {
        printf("\n\n===\n")
        printf("Starting new DoE with:\n")
        printf(paste0(rbind(paste0(names(params), ": "), 
                            paste0(params, "\n"))))
        printf("===\n\n")
        
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
        
        for (i in seq_len(length(params$to_optimize)))
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
        
        params <- fixOptParams(params, algorithm)
        params <- utilsIPO$attachList(params$to_optimize, params$no_optimization)
    }
    
    #params <- utilsIPO$attachList(params$to_optimize, params$no_optimization)
    
    return(history)
}
