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
    ret <- f(EICs, params[setdiff(names(params), c("algorithm", "forcePeakRange", "relMinIntensity"))], verbose = verbose)
    
    if (any(params$forcePeakRange != 0) || params$relMinIntensity > 0)
    {
        ret <- lapply(ret, function(pl)
        {
            if (nrow(pl) == 0)
                return(pl)
            
            pl <- copy(pl)
            
            if (params$relMinIntensity > 0)
            {
                minInt <- params$relMinIntensity * max(pl$intensity)
                pl <- pl[numGTE(intensity, minInt) == TRUE]
            }
            
            if (params$forcePeakRange[1] > 0)
                pl[, retmin := pmax(ret - params$forcePeakRange[1], retmin)]
            if (params$forcePeakRange[2] > 0)
                pl[, retmax := pmin(ret + params$forcePeakRange[2], retmax)]
        })
    }
    
    return(ret)
}

findPeaksOpenMS <- function(EICs, params, verbose = TRUE)
{
    # EICs should be a named list
    
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
                  "-algorithm:PeakPickerMRM:sgolay_frame_length" = params$SGolayFrameLength,
                  "-algorithm:PeakPickerMRM:sgolay_polynomial_order" = params$SGolayPolyOrder,
                  "-algorithm:PeakPickerMRM:use_gauss" = boolToChr(params$useGauss),
                  "-algorithm:PeakPickerMRM:gauss_width" = params$gaussWidth,
                  "-algorithm:PeakPickerMRM:signal_to_noise" = params$SN,
                  "-algorithm:PeakPickerMRM:sn_win_len" = params$SNWinLen,
                  "-algorithm:PeakPickerMRM:sn_bin_count" = params$SNBinCount,
                  "-algorithm:PeakPickerMRM:method" = params$method,
                  "-algorithm:PeakIntegrator:integration_type" = params$integrationType,
                  "-algorithm:PeakIntegrator:baseline_type" = params$baselineType,
                  "-algorithm:PeakIntegrator:fit_EMG" = boolToChr(params$fitEMG))

    if (!is.null(params[["extraOpts"]]) && length(params$extraOpst) > 0)
        settings <- modifyList(settings, params$extraOpts)
    osettings <- OpenMSArgListToOpts(settings)
    osettings <- c(osettings, "-algorithm:compute_peak_shape_metrics") # not a value, just a flag
    executeCommand(getExtDepPath("openms", "MRMTransitionGroupPicker", "OpenMS"), osettings,
                   stdout = if (verbose) "" else FALSE)
    maybePrintf("\n-----------\n")
    
    maybePrintf("Importing peaks... ")
    peaks <- setDT(parseFeatureMRMXMLFile(featsFile))
    maybePrintf("Done!\n")

    peaks[, ID := NULL]
    setnames(peaks, "chromID", "name")
    setcolorder(peaks, c("name", "ret", "retmin", "retmax", "area", "intensity"))

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
