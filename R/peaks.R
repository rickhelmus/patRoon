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
