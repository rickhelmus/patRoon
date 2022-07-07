#' @include main.R
NULL

findPeaks <- function(EICs, algorithm, ...)
{
    # UNDONE: export? If yes, add checkmate's
    
    f <- switch(algorithm,
                openms = findPeaksOpenMS,
                xcms3 = findPeaksXCMS3,
                envipick = findPeaksEnviPick)
    f(EICs, ...)
}

findPeaksOpenMS <- function(EICs, minRTDistance = 10, minNumPeaks = 5, minSNRatio = 2, resampleTraces = FALSE,
                            extraOpts = NULL, smoothWidth = NULL, intSearchRTWindow = 3)
{
    # EICs should be a named list of data.tables
    
    chromFile <- tempfile(fileext = ".mzML")
    featsFile <- tempfile(fileext = ".featureXML")
    
    printf("Exporting EICs... ")
    writeChromsToMzML(EICs, chromFile)
    printf("Done!\n")
    
    if (!is.null(smoothWidth))
    {
        printf("Smoothing traces with OpenMS...\n-----------\n")
        executeCommand(getCommandWithOptPath("NoiseFilterGaussian", "OpenMS"),
                       OpenMSArgListToOpts(c("-in" = chromFile,
                                             "-out" = chromFile,
                                             "-algorithm:gaussian_width" = smoothWidth)))
        printf("\n-----------\n")
    }
    
    printf("Finding peaks with OpenMS...\n-----------\n")
    settings <- c("-in" = chromFile,
                  "-out" = featsFile,
                  "-algorithm:min_rt_distance" = minRTDistance,
                  "-algorithm:min_num_peaks_per_feature" = minNumPeaks,
                  "-algorithm:min_signal_to_noise_ratio" = minSNRatio)
    if (resampleTraces)
        settings <- c(settings, "-algorithm:resample_traces")
    if (!is.null(extraOpts))
        settings <- modifyList(settings, extraOpts)
    executeCommand(getCommandWithOptPath("FeatureFinderMRM", "OpenMS"), OpenMSArgListToOpts(settings))
    printf("\n-----------\n")
    
    printf("Importing peaks... ")
    peaks <- setDT(parseFeatureXMLFile(featsFile))
    printf("Done!\n")
    
    # NOTE: the mz parameter is (ab)used to set the EIC index
    ret <- split(peaks, by = "mz", keep.by = FALSE)
    names(ret) <- names(EICs)
    
    ret <- pruneList(ret, checkZeroRows = TRUE)
    
    printf("Filling in intensities... ")
    # subset columns, convert to data.tables & fill in intensities (not reported by OpenMS)
    ret <- Map(ret, EICs[names(ret)], f = function(p, eic)
    {
        p <- setDT(p[, c("ret", "retmin", "retmax", "area")])
        p[, intensity := sapply(ret, function(r)
        {
            e <- eic[time %between% c(r - intSearchRTWindow, r + intSearchRTWindow)]
            return(max(e$intensity))
        })]
        return(p)
    })
    printf("Done!\n")
    
    return(ret)
}

findPeaksEnviPick <- function(EICs, ...)
{
    checkPackage("enviPick", "blosloos/enviPick")
    
    # make a dummy input list for enviPick::mzpick()
    dummyPL <- list(State = data.frame("Raw?" = TRUE, "Partioned?" = TRUE, "Clustered?" = TRUE, "Filtered?" = FALSE,
                                       "Picked?" = FALSE),
                    Parameters = data.frame(parameters = character(34), value = character(34)),
                    Results = numeric(5),
                    Scans = list(), # filled in later
                    Partition_index = matrix(),
                    EIC_index = matrix(), # filled in later
                    Peak_index = 0,
                    Peaklist = 0)
    
    # UNDONE: this could probably be faster if called with all EICs at once?
    ret <- sapply(EICs, function(eic)
    {
        dummyPL$Scans[[1]] <- eic$time
        sc <- copy(eic)
        setnames(sc, "time", "RT")
        sc[, "m/z" := 100] # dummy, no need for this
        sc[, measureID := seq_len(.N)]
        sc[, c("partID", "clustID") := 1] # fake that this trace was partitioned & clustered
        sc[, peakID := 0]
        setcolorder(sc, c("m/z", "intensity", "RT"))
        # setorderv(sc, "m/z") # needs to be ordered by mz (but we are using a single dummy value, so no need)
        dummyPL$Scans[[2]] <- as.matrix(sc)
        dummyPL$EIC_index <- as.matrix(data.table(start_ID = 1, end_ID = nrow(EIC), number_peaks = nrow(EIC)))
        
        p <- as.data.table(enviPick::mzpick(dummyPL, ...)$Peaklist[, c("RT", "minRT", "maxRT", "sum_int", "max_int")])
        setnames(p, c("ret", "retmin", "retmax", "area", "intensity"))
        return(p)
    }, simplify = FALSE)
    names(ret) <- names(EICs)
    ret <- pruneList(ret, checkZeroRows = TRUE)
    
    return(ret)
}
