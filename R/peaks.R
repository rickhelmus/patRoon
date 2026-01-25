# SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
NULL

findPeaks <- function(EICs, fillEICs, params, logPath)
{
    # UNDONE: export? If yes, add checkmate's
 
    f <- switch(params$algorithm,
                openms = findPeaksOpenMS,
                xcms3 = findPeaksXCMS3,
                envipick = findPeaksEnviPick,
                piek = findPeaksPiek)
    
    if (params$algorithm == "openms") # UNDONE: change in case we get additional algos with logging
        mkdirp(dirname(logPath))
    
    ret <- f(EICs, fillEICs, params[setdiff(names(params), c("algorithm", "forcePeakWidth", "relMinIntensity",
                                                             "calcCentroid"))],
             logPath = logPath)
    
    if (any(params$forcePeakWidth != 0) || params$relMinIntensity > 0 || params$calcCentroid != "algorithm")
    {
        ret <- Map(ret, EICs[names(ret)], f = function(pl, eic)
        {
            if (nrow(pl) == 0)
                return(pl)
            
            pl <- copy(pl)
            
            if (params$relMinIntensity > 0)
            {
                minInt <- params$relMinIntensity * max(pl$intensity)
                pl <- pl[numGTE(intensity, minInt) == TRUE]
            }
            
            if (params$forcePeakWidth[1] > 0)
            {
                pl[(retmax - retmin) < params$forcePeakWidth[1],
                   c("retmin", "retmax") := .(pmin(ret - params$forcePeakWidth[1] / 2, retmin),
                                              pmax(ret + params$forcePeakWidth[1] / 2, retmax))]
            }
            if (params$forcePeakWidth[2] > 0)
            {
                pl[(retmax - retmin) > params$forcePeakWidth[2],
                   c("retmin", "retmax") := .(pmax(ret - params$forcePeakWidth[2] / 2, retmin),
                                              pmin(ret + params$forcePeakWidth[2] / 2, retmax))]
            }
            
            if (params$calcCentroid != "algorithm")
            {
                for (row in seq_len(nrow(pl)))
                {
                    eicS <- eic[numGTETol(eic[, "time"], pl$retmin[row]) & numLTETol(eic[, "time"], pl$retmax[row]), , drop = FALSE]
                    if (nrow(eicS) > 1)
                    {
                        if (params$calcCentroid == "max")
                            set(pl, row, "ret", eicS[which.max(eicS[, "intensity"]), "time"])
                        else if (params$calcCentroid == "weighted.mean")
                            set(pl, row, "ret", weighted.mean(eicS[, "time"], eicS[, "intensity"]))
                        else # if (params$calcCentroid == "centerOfMass")
                            set(pl, row, "ret", calcCenterOfMass(eicS[, "time"], eicS[, "intensity"]))
                    }
                    
                }
            }
            
            return(pl)
        })
    }
    
    return(ret)
}

findPeaksOpenMS <- function(EICs, fillEICs, params, logPath)
{
    # EICs should be a named list
    
    TraMLFile <- tempfile(fileext = ".TraML")
    chromFile <- tempfile(fileext = ".mzML")
    featsFile <- tempfile(fileext = ".featureXML")

    unlink(logPath, force = TRUE)
    logPrintf <- function(...) cat(sprintf(...), file = logPath, append = TRUE, sep = "\n")
        
    boolToChr <- function(b) if (b) "true" else "false"
    
    logPrintf("Exporting EICs... ")
    writeTraML(names(EICs), TraMLFile)
    writeChromsToMzML(EICs, fillEICs, names(EICs), chromFile)
    logPrintf("Done!\n")
    
    logPrintf("Finding peaks with OpenMS...\n-----------\n")
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
                  "-algorithm:PeakIntegrator:fit_EMG" = boolToChr(params$fitEMG),
                  "-debug" = 10)

    if (!is.null(params[["extraOpts"]]) && length(params$extraOpst) > 0)
        settings <- modifyList(settings, params$extraOpts)
    osettings <- OpenMSArgListToOpts(settings)
    osettings <- c(osettings, "-algorithm:compute_peak_shape_metrics") # not a value, just a flag
    logPrintf("%s", executeCommand(getExtDepPath("openms", "MRMTransitionGroupPicker", "OpenMS"), osettings,
                                   stdout = TRUE, stderr = TRUE))
    logPrintf("\n-----------\n")
    
    logPrintf("Importing peaks... ")
    peaks <- setDT(parseFeatureMRMXMLFile(featsFile))
    logPrintf("Done!\n")

    peaks[, ID := NULL]
    setnames(peaks, "chromID", "name")
    setcolorder(peaks, c("name", "ret", "retmin", "retmax", "area", "intensity"))

    peaksList <- split(peaks, by = "name", keep.by = FALSE)
    peaksList <- pruneList(peaksList, checkZeroRows = TRUE)
    
    return(peaksList)
}

findPeaksXCMS3 <- function(EICs, fillEICs, params, logPath)
{
    allTime <- attr(EICs, "allXValues")
    ret <- sapply(EICs, function(eic)
    {
        inten <- if (fillEICs) doFillEIXIntensities(allTime, eic[, "time"], eic[, "intensity"]) else eic[, "intensity"]
        tim <- if (fillEICs) allTime else eic[, "time"]
        p <- as.data.table(suppressWarnings(do.call(xcms::peaksWithCentWave, c(list(inten, tim), params))))
        cols <- c("ret", "retmin", "retmax", "area", "intensity")
        setnames(p, c("rt", "rtmin", "rtmax", "into", "maxo"), cols)
        setcolorder(p, cols)
        return(p)
    }, simplify = FALSE)
    ret <- pruneList(ret, checkZeroRows = TRUE)
}

findPeaksEnviPick <- function(EICs, fillEICs, params, logPath)
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
    
    allTime <- if (fillEICs) attr(EICs, "allXValues")
    
    # UNDONE: this could probably be faster if called with all EICs at once?
    ret <- sapply(EICs, function(eic)
    {
        if (fillEICs)
        {
            dummyPL$Scans[[1]] <- allTime
            sc <- data.table(RT = allTime, intensity = doFillEIXIntensities(allTime, eic[, "time"], eic[, "intensity"]))
        }
        else
        {
            dummyPL$Scans[[1]] <- eic[, "time"]
            sc <- data.table(RT = eic[, "time"], intensity = eic[, "intensity"])
        }
        sc[, "m/z" := 100] # dummy, no need for this
        sc[, measureID := seq_len(.N)]
        sc[, c("partID", "clustID") := 1] # fake that this trace was partitioned & clustered
        sc[, peakID := 0]
        setcolorder(sc, c("m/z", "intensity", "RT"))
        # setorderv(sc, "m/z") # needs to be ordered by mz (but we are using a single dummy value, so no need)
        dummyPL$Scans[[2]] <- as.matrix(sc)
        dummyPL$EIC_index <- as.matrix(data.table(start_ID = 1, end_ID = nrow(sc), number_peaks = nrow(sc)))
        
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

findPeaksPiek <- function(EICs, fillEICs, params, logPath)
{
    setOMPThreads()
    # UNDONE: log output?
    peaks <- doFindPeaksPiek(EICs, fillEICs = fillEICs, minIntensity = params$minIntensity, SN = params$SN,
                             peakWidthMin = params$peakWidth[1], peakWidthMax = params$peakWidth[2],
                             RTMin = params$RTRange[1],
                             RTMax = if (!is.finite(params$RTRange[2])) 0 else params$RTRange[2],
                             maxPeaksPerSignal = params$maxPeaksPerSignal, verbose = FALSE)
    names(peaks) <- names(EICs)
    peaks <- lapply(peaks, setDT)
    peaks <- pruneList(peaks, checkZeroRows = TRUE)
    return(peaks)
}
