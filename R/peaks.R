#' @include main.R
NULL

findPeaks <- function(EICs, params, verbose = TRUE)
{
    # UNDONE: export? If yes, add checkmate's
    
    f <- switch(params$algorithm,
                openms = findPeaksOpenMS,
                xcms3 = findPeaksXCMS3,
                envipick = findPeaksEnviPick,
                dietrich = findPeaksDietrich)
    f(EICs, params[setdiff(names(params), "algorithm")], verbose = verbose)
}

findPeaksOpenMSOld <- function(EICs, minRTDistance = 10, minNumPeaks = 5, minSNRatio = 2, resampleTraces = FALSE,
                            extraOpts = NULL, smoothWidth = NULL, intSearchRTWindow = 3, scaleTimeFactor = NULL,
                            verbose = TRUE)
{
    # EICs should be a named list of data.tables
    
    # HACK HACK HACK: OpenMS errors if the time range is very small. For instance, this is a problem if IMS data is used
    # with findMobilities() --> just multiply everything 'time' related with 100 for now.
    if (!is.null(scaleTimeFactor))
    {
        EICs <- lapply(EICs, function(eic) copy(eic)[, time := time * scaleTimeFactor])
        minRTDistance <- minRTDistance * scaleTimeFactor
        intSearchRTWindow <- intSearchRTWindow * scaleTimeFactor
    }
    
    chromFile <- tempfile(fileext = ".mzML")
    featsFile <- tempfile(fileext = ".featureXML")
    
    maybePrintf <- if (verbose) printf else function(...) NULL
    
    maybePrintf("Exporting EICs... ")
    writeChromsToMzML(EICs, chromFile)
    maybePrintf("Done!\n")
    
    if (!is.null(smoothWidth))
    {
        maybePrintf("Smoothing traces with OpenMS...\n-----------\n")
        executeCommand(getCommandWithOptPath("NoiseFilterGaussian", "OpenMS"),
                       OpenMSArgListToOpts(c("-in" = chromFile,
                                             "-out" = chromFile,
                                             "-algorithm:gaussian_width" = smoothWidth)),
                       stdout = if (verbose) "" else FALSE)
        maybePrintf("\n-----------\n")
    }
    
    maybePrintf("Finding peaks with OpenMS...\n-----------\n")
    settings <- c("-in" = chromFile,
                  "-out" = featsFile,
                  "-algorithm:min_rt_distance" = minRTDistance,
                  "-algorithm:min_num_peaks_per_feature" = minNumPeaks,
                  "-algorithm:min_signal_to_noise_ratio" = minSNRatio)
    if (resampleTraces)
        settings <- c(settings, "-algorithm:resample_traces")
    if (!is.null(extraOpts))
        settings <- modifyList(settings, extraOpts)
    executeCommand(getCommandWithOptPath("FeatureFinderMRM", "OpenMS"), OpenMSArgListToOpts(settings),
                   stdout = if (verbose) "" else FALSE)
    maybePrintf("\n-----------\n")
    
    maybePrintf("Importing peaks... ")
    peaks <- setDT(parseFeatureXMLFile(featsFile))
    maybePrintf("Done!\n")
    
    # NOTE: the mz parameter is (ab)used to set the EIC index
    peaks[, name := names(EICs)[mz]][, mz := NULL]

    # BUG: sometimes RT reported is wrong (https://github.com/OpenMS/OpenMS/issues/6239). For now simply omit the
    # results...
    peaks <- peaks[ret %between% list(retmin, retmax)]
    
    peaksList <- split(peaks, by = "name", keep.by = FALSE)
    peaksList <- pruneList(peaksList, checkZeroRows = TRUE)
    
    maybePrintf("Post-processing... ")
    
    # subset columns, convert to data.tables & fill in intensities (not reported by OpenMS)
    peaksList <- Map(peaksList, EICs[names(peaksList)], f = function(p, eic)
    {
        p <- setDT(p[, c("ret", "retmin", "retmax", "area")])
        p[, intensity := sapply(ret, function(r)
        {
            e <- eic[time %between% c(r - intSearchRTWindow, r + intSearchRTWindow)]
            return(max(e$intensity))
        })]
        
        if (!is.null(scaleTimeFactor))
        {
            p[, c("ret", "retmin", "retmax") := .(ret / scaleTimeFactor,
                                                  retmin / scaleTimeFactor,
                                                  retmax / scaleTimeFactor)]
        }
        
        return(p)
    })
    maybePrintf("Done!\n")
    
    return(peaksList)
}

findPeaksOpenMS <- function(EICs, params, scaleTimeFactor = NULL, verbose = TRUE)
{
    # UNDONE: more parameters, check what are sensible defaults
    
    # EICs should be a named list
    
    # HACK HACK HACK: OpenMS errors if the time range is very small. For instance, this is a problem if IMS data is used
    # with findMobilities() --> just increase the scale by scaleTimeFactor for now.
    # UNDONE: do we still need this? If so, move to params
    if (!is.null(scaleTimeFactor))
    {
        EICs <- lapply(EICs, function(eic) { eic$time <- eic$time * scaleTimeFactor; return(eic) })
    }
    
    TraMLFile <- tempfile(fileext = ".TraML")
    chromFile <- tempfile(fileext = ".mzML")
    featsFile <- tempfile(fileext = ".featureXML")
    
    maybePrintf <- if (verbose) printf else function(...) NULL
    boolToChr <- function(b) if (b) "true" else "false"
    
    maybePrintf("Exporting EICs... ")
    writeTraML(names(EICs), TraMLFile)
    writeChromsToMzML(EICs, names(EICs), chromFile)
    maybePrintf("Done!\n")
    
    maybePrintf("Finding peaks with OpenMS...\n-----------\n")
    settings <- c("-in" = chromFile,
                  "-tr" = TraMLFile,
                  "-out" = featsFile,
                  "-algorithm:min_peak_width" = params$minPeakWidth,
                  "-algorithm:background_subtraction" = params$backgroundSubtraction,
                  #"-algorithm:compute_peak_shape_metrics" = boolToChr(TRUE), UNDONE: this gives an error with OpenMS 2.7 --> maybe fixed with newer versions?
                  "-algorithm:PeakPickerMRM:sgolay_frame_length" = params$SGolayFrameLength,
                  "-algorithm:PeakPickerMRM:sgolay_polynomial_order" = params$SGolayPolyOrder,
                  "-algorithm:PeakPickerMRM:use_gauss" = boolToChr(params$useGauss),
                  "-algorithm:PeakPickerMRM:signal_to_noise" = params$SN,
                  "-algorithm:PeakPickerMRM:sn_win_len" = params$SNWinLen,
                  "-algorithm:PeakPickerMRM:sn_bin_count" = params$SNBinCount,
                  "-algorithm:PeakPickerMRM:method" = params$method,
                  "-algorithm:PeakIntegrator:integration_type" = params$integrationType,
                  "-algorithm:PeakIntegrator:baseline_type" = params$baselineType,
                  "-algorithm:PeakIntegrator:fit_EMG" = boolToChr(params$fitEMG))

    if (!is.null(params[["extraOpts"]]) && length(params$extraOpst) > 0)
        settings <- modifyList(settings, params$extraOpts)
    executeCommand(getExtDepPath("openms", "MRMTransitionGroupPicker", "OpenMS"), OpenMSArgListToOpts(settings),
                   stdout = if (verbose) "" else FALSE)
    maybePrintf("\n-----------\n")
    
    maybePrintf("Importing peaks... ")
    peaks <- setDT(parseFeatureMRMXMLFile(featsFile))
    maybePrintf("Done!\n")

    peaks[, ID := NULL]
    setnames(peaks, "chromID", "name")
    setcolorder(peaks, c("name", "ret", "retmin", "retmax", "area", "intensity"))

    if (!is.null(scaleTimeFactor))
    {
        peaks[, c("ret", "retmin", "retmax") := .(ret / scaleTimeFactor,
                                                  retmin / scaleTimeFactor,
                                                  retmax / scaleTimeFactor)]
    }
    
    peaksList <- split(peaks, by = "name", keep.by = FALSE)
    peaksList <- pruneList(peaksList, checkZeroRows = TRUE)
    
    return(peaksList)
}

findPeaksXCMS3 <- function(EICs, params, verbose = TRUE)
{
    ret <- sapply(EICs, function(eic)
    {
        p <- as.data.table(do.call(xcms::peaksWithCentWave, c(list(eic$intensity, eic$time), params)))
        cols <- c("ret", "retmin", "retmax", "area", "intensity")
        setnames(p, c("rt", "rtmin", "rtmax", "into", "maxo"), cols)
        setcolorder(p, cols)
        return(p)
    }, simplify = FALSE)
    ret <- pruneList(ret, checkZeroRows = TRUE)
}

findPeaksEnviPick <- function(EICs, params, verbose = TRUE)
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
        sc <- as.data.table(eic)
        setnames(sc, "time", "RT")
        sc[, "m/z" := 100] # dummy, no need for this
        sc[, measureID := seq_len(.N)]
        sc[, c("partID", "clustID") := 1] # fake that this trace was partitioned & clustered
        sc[, peakID := 0]
        setcolorder(sc, c("m/z", "intensity", "RT"))
        # setorderv(sc, "m/z") # needs to be ordered by mz (but we are using a single dummy value, so no need)
        dummyPL$Scans[[2]] <- as.matrix(sc)
        dummyPL$EIC_index <- as.matrix(data.table(start_ID = 1, end_ID = nrow(eic), number_peaks = nrow(eic)))
        
        p <- do.call(enviPick::mzpick, c(list(dummyPL), params))$Peaklist
        if (isTRUE(all.equal(p, 0))) # no results
            return(data.table(ret = numeric(), retmin = numeric(), retmax = numeric(), area = numeric(),
                              intensity = numeric()))
        p <- as.data.table(p)[, c("RT", "minRT", "maxRT", "sum_int", "max_int"), with = FALSE]
        setnames(p, c("ret", "retmin", "retmax", "area", "intensity"))
        return(p)
    }, simplify = FALSE)
    ret <- pruneList(ret, checkZeroRows = TRUE)
    
    return(ret)
}

findPeaksDietrich <- function(EICs, params, verbose = TRUE)
{
    setOMPThreads()
    peaks <- doFindPeaksDietrich(EICs, minIntensity = params$minIntensity, SN = params$SN,
                                 peakWidthMin = params$peakWidth[1], peakWidthMax = params$peakWidth[2],
                                 RTMin = params$RTRange[1],
                                 RTMax = if (!is.finite(params$RTRange[2])) 0 else params$RTRange[2],
                                 maxPeaksPerSignal = params$maxPeaksPerSignal, verbose = verbose)
    names(peaks) <- names(EICs)
    peaks <- lapply(peaks, setDT)
    peaks <- pruneList(peaks, checkZeroRows = TRUE)
    return(peaks)
}
