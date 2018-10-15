# This file contains IPO utility functions used for optimizing feature finding.
# The functions are copied unchanged from utils.R and
# optimizeXcmsSetParameters.R. All functions are placed inside a separate
# environment to keep the package namespace clean.


# utilsIPO <- list2env(hash = TRUE, envir = env, parent =  x = list(
utilsIPO <- setRefClass("utilsIPO", methods = list(

calcPPS = function(xset, isotopeIdentification=c("IPO", "CAMERA"), ...) {
    
    isotopeIdentification <- match.arg(isotopeIdentification)
    
    ret <- vector(mode="numeric", 5) #array(0, dim=c(1,5)) 
    names(ret) <- c("ExpId", "#peaks", "#NonRP", "#RP", "PPS")
    if(is.null(xset)) {
        return(ret)
    } 
    
    if(nrow(peaks_IPO(xset)) == 0) {
        return(ret)
    }
    
    peak_source <- peaks_IPO(xset)[,c("mz", "rt", "sample", "into", "mzmin", 
                                      "mzmax", "rtmin", "rtmax"),drop=FALSE]
    ret[2] <- nrow(peak_source)
    
    if(isotopeIdentification == "IPO")
        iso_mat <- findIsotopes.IPO(xset, ...)  
    else
        iso_mat <- findIsotopes.CAMERA(xset, ...)
    
    samples <- unique(peak_source[,"sample"])
    isotope_abundance = 0.01108    
    
    #calculating low intensity peaks
    for(sample in samples) {
        non_isos_peaks <- peak_source
        
        if(nrow(iso_mat) > 0) {
            non_isos_peaks <- peak_source[-unique(c(iso_mat)),,drop=FALSE] 
        } 
        
        speaks <- non_isos_peaks[non_isos_peaks[,"sample"]==sample,,drop=FALSE]
        intensities <- speaks[,"into"]
        na_int <- is.na(intensities)
        intensities <- intensities[!na_int]
        
        if(length(intensities)>0) {
            tmp <- intensities[order(intensities)]
            int_cutoff <- mean(tmp[1:max(round((length(tmp)/33),0),1)])
            
            masses <- speaks[!na_int, "mz"]
            #floor((masses-2*CH3)/CH2) + 2
            maximum_carbon <- calcMaximumCarbon(masses)
            carbon_probabilty <- maximum_carbon*isotope_abundance
            
            iso_int <- intensities * carbon_probabilty
            
            not_loq_peaks <- sum(iso_int>int_cutoff)
            ret[3] <- ret[3] + not_loq_peaks
        }
    }#end_for_sample    
    
    ret[4] <- length(unique(c(iso_mat)))
    if(ret[3] == 0) {
        ret[5] <- (ret[4]+1)^2/(ret[3]+1)  
    } else {    
        ret[5] <- ret[4]^2/ret[3]  
    }
    
    return(ret)
    
},




findIsotopes.IPO = function(xset, checkPeakShape=c("none", "borderIntensity", "sinusCurve", 
                                                   "normalDistr")) {
        
        checkPeakShape <- match.arg(checkPeakShape)
        
        iso_mat <- matrix(0, nrow=0, ncol=2)
        if(is.null(xset)) {
            return(iso_mat)
        }
        
        colnames(iso_mat) <- c("12C", "13C")
        peak_source <- peaks_IPO(xset)[,c("mz", "rt", "sample", "into", "maxo", "mzmin",
                                          "mzmax", "rtmin", "rtmax"), drop=FALSE]
        
        for(i in 1:ncol(peak_source)) {
            peak_source <- peak_source[!is.na(peak_source[,i]),,drop=FALSE]
        }
        
        peak_source <- cbind(1:nrow(peak_source), peak_source)
        colnames(peak_source)[1] <- "id"  
        
        #carbon = 12.0
        #hydrogen	= 1.0078250170
        #CH3 = carbon + 3 * hydrogen
        #CH2 = carbon + 2 * hydrogen
        isotope_mass = 1.0033548
        isotope_abundance = 0.01108
        
        samples <- max(peak_source[,"sample"])
        
        #start_sample
        for(sample in 1:samples) { 
            #only looking into peaks from current sample   
            speaks <- peak_source[peak_source[,"sample"]==sample,,drop=FALSE]
            split <- 250
            if(!(checkPeakShape=="none"))
                rawdata <- loadRaw(xcmsSource(filepaths(xset)[sample]))
            
            if(nrow(speaks)>1) {  		      
                #speaks <- speaks[,-c("sample")]
                speaks <- speaks[order(speaks[,"mz"]),]
                
                while(!is.null(nrow(speaks)) & length(speaks) > 3) {
                    part_peaks <- NULL
                    #splitting the data into smaller pieces to improve speed    
                    if(nrow(speaks) < split) {
                        part_peaks <- speaks
                    } else {          
                        upper_bound <- speaks[split,"mzmax"] + isotope_mass          
                        end_point <- sum(speaks[,"mz"] < upper_bound)
                        part_peaks <- speaks[1:end_point,,drop=FALSE]
                    }		
                    
                    rt <- part_peaks[,"rt"]
                    rt_window <- rt * 0.005
                    rt_lower <- part_peaks[,"rt"] - rt_window
                    rt_upper <- part_peaks[,"rt"] + rt_window
                    rt_matrix <-  
                        t(matrix(rep(rt, nrow(part_peaks)), ncol=nrow(part_peaks)))
                    rt_matrix_bool <- rt_matrix >= rt_lower & rt_matrix <= rt_upper
                    
                    mz <- part_peaks[,"mz"]
                    #isotope_masses - mz_window
                    mz_lower <- part_peaks[,"mzmin"] + isotope_mass
                    #isotope_masses + mz_window
                    mz_upper <- part_peaks[,"mzmax"] + isotope_mass
                    mz_matrix <-  
                        t(matrix(rep(mz, nrow(part_peaks)), ncol=nrow(part_peaks)))
                    mz_matrix_bool <- mz_matrix >= mz_lower & mz_matrix <= mz_upper
                    
                    rt_mz_matrix_bool <- rt_matrix_bool & mz_matrix_bool
                    
                    rt_mz_peak_ids <- which(rowSums(rt_mz_matrix_bool)>0)
                    calculations <- min(split, nrow(speaks))
                    rt_mz_peak_ids <- rt_mz_peak_ids[rt_mz_peak_ids < calculations]
                    
                    for(i in rt_mz_peak_ids) {
                        current <- part_peaks[i, ,drop=FALSE]
                        rt_mz_peaks <- part_peaks[rt_mz_matrix_bool[i,],,drop=FALSE]
                        rt_difference <- 
                            abs(current[,"rt"] - rt_mz_peaks[, "rt"]) / current[,"rt"]
                        rt_mz_peaks <- cbind(rt_mz_peaks, rt_difference)
                        #test intensity_window
                        #floor((current["mz"]-2*CH3)/CH2) + 2
                        maximum_carbon <- calcMaximumCarbon(current[,"mz"]) 
                        carbon_probabilty <- c(1,maximum_carbon)*isotope_abundance
                        iso_intensity <- current[,"into"] * carbon_probabilty
                        
                        int_bools <- 
                            rt_mz_peaks[,"into"] >= iso_intensity[1] & 
                            rt_mz_peaks[,"into"] <= iso_intensity[2]
                        
                        if(sum(int_bools) > 0) {
                            int_peaks <- rt_mz_peaks[int_bools,,drop=FALSE]
                            boundary_bool <- rep(TRUE, (nrow(int_peaks)+1))
                            if(!(checkPeakShape=="none")) {
                                if(checkPeakShape=="borderIntensity") {
                                    boundary_bool <- checkIntensitiesAtRtBoundaries(
                                        rawdata, 
                                        rbind(current,int_peaks[,-ncol(int_peaks), drop=FALSE]))
                                } else {
                                    if(checkPeakShape=="sinusCurve") {                
                                        boundary_bool <- checkSinusDistribution(
                                            rawdata, 
                                            rbind(current,int_peaks[,-ncol(int_peaks),drop=FALSE]))
                                    } else {                  
                                        boundary_bool <- checkNormalDistribution(
                                            rawdata, 
                                            rbind(current,int_peaks[,-ncol(int_peaks),drop=FALSE]))
                                    }
                                }
                            } #else {
                            #boundary_bool <- rep(TRUE, (nrow(int_peaks)+1))
                            #}              
                            if(boundary_bool[1] & sum(boundary_bool[-1])>0) {                 
                                iso_peaks <- int_peaks[boundary_bool[-1],,drop=FALSE]
                                iso_id <- 
                                    iso_peaks[which.min(iso_peaks[,"rt_difference"]), "id"]
                                #iso_list[[length(iso_list)+1]] <- c(current[,"id"], iso_id)            
                                iso_mat <- rbind(iso_mat, c(current[,"id"], iso_id))                
                            }
                        }
                    }
                    speaks <- speaks[-(1:calculations),]		    
                    
                }#end_while_sample_peaks 
            }
        }
        return(iso_mat)
    },

# checking intensities at rtmin and rtmax. peaks[,"maxo"] must be at least 
# double as high does not work for retention time corrected data 
checkIntensitiesAtRtBoundaries = function(rawdata, 
                                          peaks, 
                                          minBoundaryToMaxo=1/3, 
                                          ppmstep=15) {
        ret <- rep(TRUE, nrow(peaks))
        for(i in 1:nrow(peaks)) {
            peak <- peaks[i,]
            for(boundary in c("rtmin", "rtmax")) {
                rtIndex <- which(rawdata$rt==peak[boundary])
                if(length(rtIndex)>0) {
                    if(rtIndex==length(rawdata$scanindex)) {
                        rtIndices <- c(rawdata$scanindex[rtIndex], length(rawdata$mz))
                    } else {
                        rtIndices <- rawdata$scanindex[c(rtIndex, rtIndex+1)]
                    }
                    
                    #only relevant mz and intensity values regarding retention time
                    mz <- rawdata$mz[(rtIndices[1]+1):rtIndices[2]]	
                    intensities <- rawdata$intensity[(rtIndices[1]+1):rtIndices[2]]
                    
                    ppm <- peak[c("mzmin", "mzmax")]*ppmstep/1000000
                    mzIntensities <- 
                        c(0, intensities[mz>=peak["mzmin"]-ppm[1] & mz<=peak["mzmax"]+ppm[2]])
                    maxBoundaryIntensity <- max(mzIntensities)
                    ret[i] <- ret[i] & maxBoundaryIntensity<peak["maxo"]*minBoundaryToMaxo
                }
            }
        }
        
        return(ret)
        
    },

checkSinusDistribution = function(rawdata, peaks) {
    ret <- rep(TRUE, nrow(peaks))
    for(i in 1:nrow(peaks)) {
        ret[i] <- testSinusDistribution(rawdata, peaks[i,,drop=FALSE])
    }
    
    return(ret)
},

checkNormalDistribution = function(rawdata, peaks) {
    ret <- rep(TRUE, nrow(peaks))
    for(i in 1:nrow(peaks)) {
        ret[i] <- testNormalDistribution(rawdata, peaks[i,,drop=FALSE])
    }
    
    return(ret)
},

getIntensitiesFromRawdata = function(rawdata, peak) {
    rt <- rawdata$rt >= peak[,"rtmin"] & rawdata$rt <= peak[,"rtmax"]
    
    rtRange <- c(min(which(rt)), max(which(rt))+1)  
    scanIndices <- 
        rawdata$scanindex[rtRange[1]:min(rtRange[2], length(rawdata$scanindex))]
    #  scanIndices <- scanIndices[!is.na(scanIndices)]
    if(rtRange[2]>length(rawdata$scanindex)) {
        scanIndices <- c(scanIndices, length(rawdata$intensity))
    }
    
    if(length(scanIndices) < 3)
        return(FALSE)  
    
    y <- c()
    for(i in 1:(length(scanIndices)-1)) {
        scanRange <- c(scanIndices[i]+1, scanIndices[i+1])
        mz <- rawdata$mz[scanRange[1]:scanRange[2]]
        y <- 
            c(y, 
              max(0, (rawdata$intensity[scanRange[1]:scanRange[2]][
                  mz >= peak[,"mzmin"] & mz <= peak[,"mzmax"]])
              )
            )
    }
    
    y
},

testNormalDistribution = function(rawdata, peak) {
    
    y <- getIntensitiesFromRawdata(rawdata, peak)
    if(length(y) < 3) {
        return(FALSE)
    }
    
    if(max(y)==0) {
        return(FALSE)
    }
    
    normY <- (y-min(y))/(max(y)-min(y))
    
    mean=10; 
    sd=3;
    
    seqModel <- seq(-4,4,length=length(normY))*sd + mean
    yModel <- dnorm(seqModel,mean,sd)
    yModel = yModel* (1/max(yModel))
    correlation <- cor(yModel, normY)
    
    correlation > 0.7
    
    
},

testSinusDistribution = function(rawdata, peak) {
    
    y <- getIntensitiesFromRawdata(rawdata, peak)
    if(length(y) < 3) {
        return(FALSE)
    }
    if(max(y)==0) {
        return(FALSE)
    }
    
    normY <- (y-min(y))/(max(y)-min(y))
    sinCurve <- (sin(seq(-pi/2,pi+1.5,length=length(normY))) + 1) / 2
    correlation <- cor(sinCurve, normY)
    
    correlation > 0.7
    
},

findIsotopes.CAMERA = function(xset, ...) {
        
        iso_mat <- matrix(0, nrow=0, ncol=2)
        if(is.null(xset)) {
            return(iso_mat)
        }
        
        ids <- peaks_IPO(xset)[,"sample", drop=FALSE]
        ids <- cbind(1:length(ids), ids)
        
        xsets <- split(xset, unique(peaks_IPO(xset)[,"sample"]))
        samples <- unique(peaks_IPO(xset)[,"sample"])
        for(sample in samples) {
            an <- xsAnnotate(xset, sample=sample)
            isos <- findIsotopes(an, ...)@isoID[,c("mpeak", "isopeak"), drop=FALSE]
            #start_id <- ids[ids[,2]==sample,,drop=FALSE][1,1] - 1
            iso_mat <- rbind(iso_mat, matrix(ids[ids[,2]==sample,1][isos], ncol=2))
        }
        
        iso_mat
    },


calcMaximumCarbon = function(masses) {  
        
        carbon = 12.0
        hydrogen  = 1.0078250170
        CH3 = carbon + 3 * hydrogen
        CH2 = carbon + 2 * hydrogen  
        
        maximum_carbon <- floor((masses-2*CH3)/CH2) + 2
        
    },   


resultIncreased = function(history) {
        
        index = length(history)
        if(history[[index]]$PPS["PPS"] == 0 & index == 1)
            stop(paste("No isotopes have been detected,",
                       "peak picking not optimizable by IPO!"))
        
        if(index < 2)
            return(TRUE)
        
        if(history[[index-1]]$PPS["PPS"] >= history[[index]]$PPS["PPS"])
            return(FALSE)
        
        return(TRUE)
        
    },

### utils.cpp

attachList = function(params_1, params_2) {
    params <- params_1
    for(factor in params_2)
        params[[length(params)+1]] <- factor
    
    names(params) <- c(names(params_1), names(params_2))
    return(params)
},


checkParams = function(params, 
                        quantitative_parameters,
                        qualitative_parameters, 
                        unsupported_parameters) { 
    if(length(typeCastParams(params)$to_optimize)==0) {
        stop("No parameters for optimization specified; stopping!")  
    }
    
    for(i in 1:length(params)) {
        param <- params[[i]]
        name <- names(params)[i]
        if(name %in% unsupported_parameters) {
            stop(paste("The parameter", name, "is not supported! Please remove
                       from parameters; stopping!"))
        }
        if(name %in% qualitative_parameters) {
            if(length(param) == 0) {
                stop(paste("The parameter", name, "has no value set!
                           Please specify; stopping!"))
            }
            if(length(param) > 1) {
                stop(paste("Optimization of parameter", name, "not supported!
                           Please specify only one value; stopping!"))
            }
        }
        if(name %in% quantitative_parameters) {
            if(length(param) == 0) {
                stop(paste("The parameter", name, "has no value set!
                           Please specify between one and two; stopping!"))
            } 
            if(length(param) > 2) {
                stop(paste("Too many values for parameter", name, "!
                           Please specify only one or two; stopping!"))
            }
            if(!all(diff(param) > 0)) {
                stop(paste("Parameter", name, "has wrong order!",
                           "Please specify in increasing order; stopping!"))
            }
        }
    }
    missing_params <- 
        which(!(c(quantitative_parameters, qualitative_parameters) %in% 
                    names(params)))
    if(length(missing_params > 0)) {
        stop(paste("The parameter(s)", 
                   paste(c(quantitative_parameters,
                           qualitative_parameters)[missing_params], 
                         collapse=", "), 
                   "is/are missing! Please specify; stopping!"))
    }
    
},


combineParams = function(params_1, params_2) {
    len <- max(unlist(sapply(params_1, length)))
    #num_params <- length(params_1)
    
    p_names <- c(names(params_1), names(params_2))
    matchedFilter <- "fwhm" %in% p_names
    
    for(i in 1:length(params_2)) {
        new_index <- length(params_1) + 1
        fact <- params_2[[i]]
        params_1[[new_index]] <- fact
        if(matchedFilter) {
            if(p_names[new_index] == "sigma" && fact == 0) {
                # update values for sigma if zero
                if("fwhm" %in% names(params_1)) {
                    params_1[[new_index]][1:len] <- params_1$fwhm/2.3548
                } else {
                    params_1[[new_index]][1:len] <- params_2$fwhm/2.3548
                }
            } else if(p_names[new_index] == "mzdiff" && fact == 0) {
                # update values for mzdiff if zero
                if("step" %in% names(params_1)) {
                    if("steps"  %in% names(params_1)) {
                        params_1[[new_index]][1:len] <- 0.8-params_1$step*params_1$steps
                    } else {
                        params_1[[new_index]][1:len] <- 0.8-params_1$step*params_2$steps
                    }	
                } else {
                    if("steps"  %in% names(params_1)) {
                        params_1[[new_index]][1:len] <- 0.8-params_2$step*params_1$steps
                    } else {
                        params_1[[new_index]][1:len] <- 0.8-params_2$step*params_2$steps
                    }	
                }
            } else {  
                # standard: replicate value
                params_1[[new_index]][1:len] <- fact
            }
        } else {
            # standard: replicate value
            params_1[[new_index]][1:len] <- fact
        }
    } 
    names(params_1) <- p_names   
    return(params_1)
    
},

createModel = function(design, params, resp) {
    # add response to the design, which gives the data for the model
    design$resp <- resp
    if(length(params) > 1) {
        # create full second order (SO) model
        # use xi in model, instead of parameter names
        formula <- 
            as.formula(paste("resp ~ SO(", 
                             paste("x", 1:length(params), 
                                   sep="", collapse=","),
                             ")", sep=""))
        model <- rsm::rsm(formula, data=design) 
    } else {
        # create full second order model with one parameter
        # here: use parameter name in model
        param_name <- names(params)[1]
        formula <- as.formula(paste("resp ~ ", param_name, " + ", 
                                    param_name, " ^ 2", sep="")) 
        model <- lm(formula, data=design) 
        model$coding <- list(x1=as.formula(paste(param_name, "~ x1"))) 
        names(model$coding) <- param_name
        #attr(model, "class") <- c("rsm", "lm")
    }
    return(model)  
},


decode = function(value, bounds) {
    if(is.na(value))
        value <- 1
    x <- (value+1)/2 # from [-1,1] to [0, 1]
    x <- (x*(max(bounds)-min(bounds))) + min(bounds)
    
    return(x)
},


decodeAll = function(values, params) {
    
    ret <- rep(0, length(params))
    for(i in 1:length(params))
        ret[i] <- decode(values[i], params[[i]])
    
    names(ret) <- names(params)
    
    return(ret)
},

encode = function(value, bounds) {
    x <- (value - min(bounds)) / (max(bounds) - min(bounds))
    
    return(x*2-1)
},


getBbdParameter = function(params) {
    
    lower_bounds <- unlist(lapply(X=params, FUN=function(x) x[1]))
    higher_bounds <- unlist(lapply(X=params, FUN=function(x) x[2]))
    
    steps <- (higher_bounds - lower_bounds)/2
    
    x <- paste("x", 1:length(params), " ~ (", c(names(params)), " - ", 
               (lower_bounds + steps), ")/", steps, sep="")
    formulae <- list()
    for(i in 1:length(x))
        formulae[[i]] <- as.formula(x[i])  
    
    design <- bbd(length(params), n0 = 1, randomize = FALSE, coding = formulae)
    
    return(design)
    
},


getCcdParameter = function(params) {
    
    lower_bounds <- unlist(lapply(X=params, FUN=function(x) x[1]))
    higher_bounds <- unlist(lapply(X=params, FUN=function(x) x[2]))
    
    steps <- (higher_bounds - lower_bounds)/2
    
    # formula for each parameter, that transformes values from the range
    # to [0, 1]
    x <- paste("x", 1:length(params), " ~ (", c(names(params)), " - ", 
               (lower_bounds + steps), ")/", steps, sep="")
    
    # list with single formulas as entries
    formulae <- list()
    for(i in 1:length(x))
        formulae[[i]] <- as.formula(x[i])  
    
    design <- rsm::ccd(length(params), # number of variables
                       n0 = 1, # number of center points
                       alpha = "face", # position of the ‘star’ points
                       randomize = FALSE, 
                       inscribed = TRUE, # TRUE: axis points are at +/- 1 and the
                       # cube points are at interior positions
                       coding = formulae) # List of coding formulas for the design
    # variables
    return(design)
    
},



getMaximumLevels = function(model) {  
    # dimension of the modeled space
    dimensions <- length(model$coding)
    
    #slices <- getSlices(dimensions-2)
    #mat <- getResponses(slices, model)
    #testdata <- getTestData(dimensions)
    #if(dimensions==1)
    #  names(testdata) <- names(model$coding)
    #
    #return(getMaxSettings(testdata, model))
    
    # define grid, to test for maximum
    if(dimensions > 6) {
        testSpace <- seq(-1,1,0.2) # 11 points
    } else { 
        testSpace <- seq(-1,1,0.1) # 21 points
    }
    
    testSize <- 10^6 # maximum number of grid points for one test
    # amount for testing each point in testSpace
    testAmount <- length(testSpace)^dimensions 
    i <- 1
    max <- rep(-1, dimensions+1) # start maximum response + setting
    # if there are more test points (=testAmount), than testSize,
    # then the tests are split and each subset is tested seperately
    while(i < testAmount) {
        testdata <- expand.grid.subset(i:(i+testSize), testSpace, dimensions)
        if(dimensions==1)
            names(testdata) <- names(model$coding)
        max_tmp <- getMaxSettings(testdata, model)
        if(max_tmp[1]>max[1]) # if better solution (i.e. test response)
            max <- max_tmp
        i <- i + testSize + 1
    }
    
    return(max)
    
},

getMaxSettings = function(testdata, model) {
    
    response <- predict(model, testdata)
    max_response <- max(response)
    # select row(s) corresponding to max
    max_settings <- testdata[response==max_response,,drop=FALSE]
    ret <- max_response
    
    for(i in 1:ncol(testdata)) {
        levels <- max_settings[,i] # all settings of variable i
        if(all(c(-1,1) %in% levels)) # if both borders are in maximum settings
            ret <- cbind(ret, NA)
        else
            ret <- cbind(ret,levels[1]) # take first setting
    }
    
    colnames(ret) <- c("response", paste("x", 1:ncol(testdata), sep=""))
    return(ret)
},


expand.grid.subset  = function(subset, sequence, dimensions) { 
    # generate a list, with sequence for each dimension
    vars <- list()
    for(i in 1:dimensions) {
        vars[[i]] <- sequence
    }
    names(vars) <- paste("x", 1:dimensions, sep="")
    
    # number of points in sequence grid
    maximumSubset <- length(sequence)^dimensions 
    # from min(subset)) to min(maximumSubset, max(subset)) OR
    # from maximumSubset to maximumSubset
    subset <- min(maximumSubset,min(subset)):min(maximumSubset, max(subset))
    
    #nv <-  #length(vars) 
    # number of values per variable = length(sequence)
    lims <- sapply(vars,length) 
    stopifnot(length(lims) > 0, # i.e. dimensions > 0
              subset <= prod(lims), # i.e. subset <= maximumSubset
              length(names(vars)) == dimensions) # i.e. dimensions = dimensions
    # empty list of length names(vars)
    res <- structure(vector("list",dimensions), .Names = names(vars))
    
    if (dimensions > 1) {
        for(i in dimensions:2) { # count down to 2: set up grid top down
            # %% = mod, %/% = integer division
            f <- prod(lims[1:(i-1)]) # number of grid points up to variable nr. (i-1)
            # repeat each element on grid 1:f
            res[[i]] <- vars[[i]][(subset - 1)%/%f + 1] 
            subset <- (subset - 1)%%f + 1 
        } 
    }
    res[[1]] <- vars[[1]][subset] 
    as.data.frame(res) 
},


peaks_IPO = function(xset) {
    # peaks function, to work with older xcms-version (<2.99.7)
    # sample column is missing, if first the first sample processed has no peaks
    # see https://github.com/sneumann/xcms/issues/220
    peaks_act <- xcms::peaks(xset)
    if (!("sample" %in% colnames(peaks_act))) {
        colnames(peaks_act)[colnames(peaks_act) == ""] <- "sample"
    }
    peaks_act
},

plotContours = function(model, maximum_slice, plot_name = NULL) {
    # generate for all variable combinations formulas
    # (which will give the plots)
    plots <- c()
    for(i in 1:(length(maximum_slice)-1)) {
        for(j in (i+1):length(maximum_slice)) {
            plots <- c(plots, as.formula(paste("~ x", i, "* x", j, sep="")))
        } 
    }
    
    # determine number of plot rows and column on single device
    plot_rows <- round(sqrt(length(plots)))
    plot_cols <- 
        if(plot_rows==1){
            length(plots)
        } else {
            ceiling(sqrt(length(plots)))
        }
    
    # save as jpeg, if plot_name is given
    if(!is.null(plot_name)) {
        plot_name = paste(plot_name, ".jpg", sep="")
        jpeg(plot_name, width=4*plot_cols, height=2*plot_rows+2, 
             units="in", res=c(200,200))
    } # otherwise plot on device
    
    op <- par(mfrow = c(plot_rows, plot_cols), oma = c(0,0,0,0))  
    # contour.lm is called
    ret_tr <- try({#may produce errors, if figure margins are too small
        contour(model, plots, image = TRUE, at = maximum_slice)
    })
    if (class(ret_tr) == "try-error") {
        message("Error in plotting (see above), but IPO continues to work normally!")
    }
    
    if (!is.null(plot_name)) {
        dev.off()
    }
    par(op)
},


toMatrix = function(data) {
    
    if(!is.matrix(data)) {
        tmp <- names(data)
        data <- matrix(data, nrow=1)
        colnames(data) <- tmp  
    } 
    
    return(data)  
    
},


typeCastParams = function(params) {
    ret_1 <- list()
    ret_2 <- list()  
    ret <- list()
    for(i in  1:length(params)) {
        factor <- params[[i]]
        if(length(factor) == 2) {
            ret_1[[(length(ret_1)+1)]] <- factor
            names(ret_1)[length(ret_1)] <- names(params)[i]
        } else {	
            ret_2[[(length(ret_2)+1)]] <- factor
            names(ret_2)[length(ret_2)] <- names(params)[i]
        }
    }	
    ret$to_optimize <- ret_1
    ret$no_optimization <- ret_2
    
    return(ret)
}


))()

