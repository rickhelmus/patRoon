#' @include main.R
NULL

loadSpectra <- function(path, rtRange = NULL, verbose = TRUE, cacheDB = NULL)
{
    hash <- makeHash(makeFileHash(path), rtRange)
    ret <- loadCacheData("specData", hash, cacheDB)
    if (is.null(ret))
    {
        if (verbose)
            printf("Loading raw spectra for '%s'...\n", path)
        msf <- mzR::openMSfile(path)
        hd <- as.data.table(mzR::header(msf))

        if (is.null(rtRange))
            ps <- mzR::peaks(msf) # load all
        else
            ps <- mzR::peaks(msf, hd[numGTE(retentionTime, rtRange[1]) & numLTE(retentionTime, rtRange[2]), seqNum])

        ret <- list(header = hd, spectra = ps)
        mzR::close(msf)
        saveCacheData("specData", ret, hash, cacheDB)
    }

    return(ret)
}

getSpectraHeader <- function(spectra, rtRange, MSLevel, precursor, precursorMzWindow)
{
    hd <- spectra$header[numGTE(retentionTime, rtRange[1]) & numLTE(retentionTime, rtRange[2]) & msLevel == MSLevel]

    if (!is.null(precursor) && !is.null(precursorMzWindow))
        hd <- hd[numLTE(abs(precursorMZ - precursor), precursorMzWindow)]

    return(hd)
}

getEIC <- function(spectra, rtRange, mzRange, MSLevel = 1, precursor = NULL, precursorMzWindow = NULL)
{
    hd <- getSpectraHeader(spectra, rtRange, MSLevel, precursor, precursorMzWindow)
    return(data.table(time = hd$retentionTime,
                      intensity = sapply(spectra$spectra[hd$seqNum],
                                         function(s) sum(s[numGTE(mz, mzRange[1]) & numLTE(mz, mzRange[2]), intensity]))))
}

getEICsForFGroups <- function(fGroups, rtWindow, mzWindow, topMost, onlyPresent)
{
    if (length(fGroups) == 0)
        return(list())

    gNames <- names(fGroups)
    gTable <- groups(fGroups)
    gInfo <- groupInfo(fGroups)
    ftind <- groupFeatIndex(fGroups)
    fTable <- featureTable(fGroups)
    anaInfo <- analysisInfo(fGroups)

    # load EICs per analysis: we don't want to load multiple potentially large
    # analysis files simultaneously. Before that, it's more efficient to first
    # figure out for which feature groups EICs have to be generated per
    # analysis. Furthermore, collect all retention ranges for EICs as these also
    # have to be checked on a per group basis.


    doEICs <- list()
    rtRanges <- list()

    getFTCol <- function(fg, anai, col)
    {
        if (ftind[[fg]][anai] == 0)
            NA
        else
            fTable[[anaInfo$analysis[anai]]][ftind[[fg]][anai]][[col]]
    }

    for (fg in gNames)
    {
        if (!is.null(topMost))
        {
            oint <- order(-gTable[[fg]])
            topAnalysesInd <- oint[seq_len(topMost)]
        }
        else
            topAnalysesInd <- seq_len(nrow(anaInfo))
        
        if (onlyPresent)
            topAnalysesInd <- topAnalysesInd[gTable[[fg]][topAnalysesInd] != 0]

        rtrMins <- sapply(topAnalysesInd, function(anai) getFTCol(fg, anai, "retmin"))
        rtrMins <- rtrMins[!is.na(rtrMins)]
        rtrMaxs <- sapply(topAnalysesInd, function(anai) getFTCol(fg, anai, "retmax"))
        rtrMaxs <- rtrMaxs[!is.na(rtrMaxs)]
        
        if (length(rtrMins) > 0 && length(rtrMaxs) > 0)
            rtr <- c(min(rtrMins), max(rtrMaxs))
        else
            rtr <- gInfo[fg, "rts"]
        
        rtr <- rtr + c(-rtWindow, rtWindow)

        doEICs[[fg]] <- anaInfo$analysis[topAnalysesInd]
        rtRanges[[fg]] <- rtr
    }

    if (!is.null(topMost))
        topMost <- min(topMost, nrow(anaInfo))

    cacheDB <- openCacheDBScope()

    EICs <- list()
    for (anai in seq_len(nrow(anaInfo)))
    {
        ana <- anaInfo$analysis[anai]
        anaHash <- NULL
        dfile <- getMzMLOrMzXMLAnalysisPath(ana, anaInfo$path[anai])

        doGroups <- gNames[sapply(gNames, function(fg) ana %in% doEICs[[fg]])]
        rtr <- rtRanges[doGroups]
        mzr <- sapply(doGroups, function(fg) gInfo[fg, "mzs"] + c(-mzWindow, mzWindow), simplify = FALSE)
        names(mzr) <- doGroups

        EICHashes <- sapply(doGroups, function(fg)
        {
            if (is.null(anaHash))
                anaHash <<- makeFileHash(dfile)
            makeHash(anaHash, rtr[[fg]], mzr[[fg]])
        })

        gEICs <- pruneList(sapply(doGroups, function(fg) loadCacheData("mzREIC", EICHashes[[fg]], cacheDB),
                                  simplify = FALSE))

        EICGNames <- doGroups[!doGroups %in% names(gEICs)]
        if (length(EICGNames) > 0)
        {
            spectra <- loadSpectra(dfile, verbose = FALSE, cacheDB = cacheDB)
            rtr <- rtr[EICGNames]
            mzr <- mzr[EICGNames]
            gEICs <- c(gEICs, setNames(loadEICs(spectra, rtr, mzr), EICGNames))

            for (fg in EICGNames)
                saveCacheData("mzREIC", gEICs[[fg]], EICHashes[[fg]], cacheDB)
        }

        EICs[[ana]] <- gEICs
    }

    return(EICs)
}

averageSpectraMZR <- function(spectra, hd, clusterMzWindow, topMost, minIntensityPre,
                              minIntensityPost, avgFun, method, precursor,
                              pruneMissingPrecursor, retainPrecursor)
{
    if (nrow(hd) == 0) # no spectra, return empty spectrum
        return(emptyMSPeakList())

    sp <- spectra$spectra[hd$seqNum]
    # convert to peaklist format
    sp <- lapply(sp, function(spec) setnames(as.data.table(spec), c("mz", "intensity")))
    sp <- lapply(sp, assignPrecursorToMSPeakList, precursor)

    return(averageSpectra(sp, clusterMzWindow, topMost, minIntensityPre, minIntensityPost,
                          avgFun, method, pruneMissingPrecursor, retainPrecursor))
}
