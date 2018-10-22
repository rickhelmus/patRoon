#' @include utils-IPO.R
#' @include doe-optimizer.R
#' @include features.R
#' @include main.R
NULL

featuresOptimizer <- setRefClass("featuresOptimizer", contains = c("DoEOptimizer" ,"VIRTUAL"),
                                 fields = list(isoIdent = "character"))

featuresOptimizer$methods(
 
    # dummy methods to be potentially overrided
    convertOptToCallParams = function(params) params,
    
    # Adapted from IPO: add OpenMS isotope detection
    calcPPS = function(feat, ...) # UNDONE: handle ...
    {
        fTable <- featureTable(feat)
        
        ret <- list(featureCount = length(feat), nonRP = 0, RP = 0, PPS = 0)
        
        if (length(feat) == 0)
            return(ret)
        
        doOpenMS <- isoIdent == "OpenMS"
        if (!doOpenMS) # no need to find isotopes with OpenMS algo
        {
            xset <- getXcmsSet(feat, TRUE)
            peak_source <- utilsIPO$peaks_IPO(xset)[, c("mz", "rt", "sample", "into", "mzmin", 
                                                        "mzmax", "rtmin", "rtmax"), drop = FALSE]
            if(isoIdent == "IPO")
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
                    non_isos_peaks <- peak_source[-unique(c(iso_mat)), , drop = FALSE] 
                
                speaks <- non_isos_peaks[non_isos_peaks[,"sample"]==anai, , drop = FALSE]
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
                ret$nonRP <- ret$nonRP + not_loq_peaks
            }
        }#end_for_sample    
        
        if (doOpenMS)
        {
            # isocount represent the number of isotopes collapsed in a feature. When
            # it's one, no other isotopes are present and hence its a NP. The
            # remaining (except the last, hence - 1) are assumed to be RP.
            ret$RP <- sum(unlist(lapply(fTable, function(ft) ft[isocount > 1, isocount - 1])))
        }
        else
            ret$RP <- length(unique(c(iso_mat)))
        
        if (ret[3] == 0)
            ret$PPS <- (ret$RP+1)^2/(ret$nonRP+1)  
        else
            ret$PPS <- ret$RP^2/ret$nonRP
        
        ret$score <- ret$PPS
        
        return(ret)
    },
    
    calculateResponse = function(params, task, final = FALSE)
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
        
        params <- convertOptToCallParams(params)
        feat <- do.call(findFeatures, c(list(anaInfo, algorithm), params))
        ret <- calcPPS(feat)
        
        if (final) # store optimized features object
            ret$object <- feat
        
        return(ret)
    },
    
    resultIncreased = function(history)
    {
        index <- length(history)
        if(history[[index]]$finalResponse$PPS == 0 & index == 1)
            stop("No isotopes have been detected!")
        if(index < 2)
            return(TRUE)
        if(history[[index-1]]$finalResponse$PPS >= history[[index]]$finalResponse$PPS)
            return(FALSE)
        return(TRUE)
    }
)



#' @export
featuresOptimization <- setClass("featuresOptimization",
                                 slots = c(algorithm = "character",
                                           startParams = "list", finalResults = "list",
                                           experiments = "list"))

#' @describeIn featuresOptimization Returns the algorithm that was used for finding features.
setMethod("algorithm", "featuresOptimization", function(obj) obj@algorithm)

#' @describeIn featuresOptimization Obtain total number of experimental design iteratations performed.
#' @export
setMethod("length", "featuresOptimization", function(x) length(x@experiments))

#' @describeIn featuresOptimization Shows summary information for this object.
#' @export
setMethod("show", "featuresOptimization", function(object)
{
    printf("A features optimization object ('%s')\n", class(object))
    printf("Algorithm: %s\n", algorithm(object))
    printf("Experimental designs performed: %d\n", length(object))
    printf("Starting params:\n"); printf("- %s: %s\n", names(object@startParams), object@startParams)
    printf("Optimized params:\n"); printf("- %s: %s\n", names(object@finalResults$parameters), object@finalResults$parameters)
    
    br <- object@finalResults$result
    br <- br[!names(br) %in% "ExpId"]
    printf("Best results: "); cat(paste(names(br), br, sep = ": ", collapse = "; ")); cat("\n")
    
    printf("\nOptimized object:\n---\n"); show(object@finalResults$object); cat("---\n")
    
    showObjectSize(object)
})

#' @describeIn featuresOptimization Returns the \code{\link{features}} object obtained with the optimized parameters.
#' @export
setMethod("getFeatures", "featuresOptimization", function(obj) obj@finalResults$features)

#' @export
setMethod("plot", "featuresOptimization", function(x, index, paramsToPlot = NULL, maxCols = NULL, type = "contour",
                                                   image = TRUE, contours = "colors", ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertInt(index, lower = 1, upper = length(x))
    checkmate::assert(checkmate::checkList(paramsToPlot, types = "character", any.missing = FALSE, min.len = 1, null.ok = TRUE),
                      checkmate::checkCharacter(paramsToPlot, min.chars = 1, len = 2),
                      checkmate::checkNull(paramsToPlot),
                      .var.name = "paramsToPlot")
    checkmate::assertCount(maxCols, positive = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertChoice(type, c("contour", "image", "persp"), add = ac)
    checkmate::assertFlag(image, add = ac)
    checkmate::assert(checkmate::checkFlag(contours),
                      checkmate::checkCharacter(contours),
                      checkmate::checkList(contours),
                      .var.name = "contours")
    checkmate::reportAssertions(ac)

    ex <- x@experiments[[index]]
    
    if (is.null(paramsToPlot))
    {
        paramsToPlot <- list()
        optNames <- names(ex$params$to_optimize)
        for (i in seq(1, length(optNames)-1))
        {
            for (j in seq(i+1, length(optNames)))
                paramsToPlot <- c(paramsToPlot, list(c(optNames[i], optNames[j])))
        }
    }
    else if (is.character(paramsToPlot))
        paramsToPlot <- list(paramsToPlot)
    
    codedNames <- names(ex$design)
    decodedNames <- rsm::truenames(ex$design)
    
    forms <- lapply(paramsToPlot, function(pn)
    {
        # change to coded names
        pn <- sapply(pn, function(n) codedNames[decodedNames == n])
        return(as.formula(paste(pn[2], "~", pn[1])))
    })
    
    maxSlice <- ex$max_settings[1, -1]
    maxSlice[is.na(maxSlice)] <- 1
    
    formsLen <- length(forms)
    if (formsLen > 1) # multiple plots?
    {
        if (is.null(maxCols))
            maxCols <- ceiling(sqrt(formsLen))
        
        if (formsLen <= maxCols)
        {
            cols <- formsLen
            rows <- 1
        }
        else
        {
            cols <- maxCols
            rows <- ceiling(formsLen / cols)
        }
        
        withr::local_par(list(mfrow = c(rows, cols)))
    }
    
    switch(type,
           contour = contour(ex$model, forms, image = image, at = maxSlice, ...),
           image = image(ex$model, forms, at = maxSlice, ...),
           persp = persp(ex$model, forms, contours = contours, at = maxSlice, ...))
})

