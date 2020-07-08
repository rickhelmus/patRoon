#' @include main.R
NULL

loadSpectra <- function(path, rtRange = NULL, verbose = TRUE, cacheDB = NULL)
{
    hash <- makeHash(makeFileHash(path), rtRange)
    ret <- loadCacheData("specData", hash, cacheDB)
    if (!is.null(ret) && length(ret$spectra) > 1 && is.data.table(ret$spectra[[1]]))
        ret <- NULL # old (pre v1.1) format, ignore cache to avoid crashes with Rcpp interface
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

setMethod("getEICsForFGroups", "featureGroups", function(fGroups, rtWindow, mzExpWindow, topMost, topMostByRGroup,
                                                         onlyPresent)
{
    if (length(fGroups) == 0)
        return(list())

    gNames <- names(fGroups)
    gTable <- groupTable(fGroups)
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

    if (!is.null(topMost))
        topMost <- min(topMost, nrow(anaInfo))
    
    EICInfo <- rbindlist(sapply(gNames, function(fg)
    {
        if (!is.null(topMost))
        {
            if (topMostByRGroup)
            {
                tbl <- data.table(int = gTable[[fg]], group = anaInfo$group, anaInd = seq_len(nrow(anaInfo)))
                tbl[, rank := frank(-int, ties.method = "first"), by = "group"]
                topAnalysesInd <- tbl[rank <= topMost]$anaInd
            }
            else
            {
                oint <- order(gTable[[fg]], decreasing = TRUE)
                topAnalysesInd <- oint[seq_len(topMost)]
            }
        }
        else
            topAnalysesInd <- seq_len(nrow(anaInfo))
        
        if (onlyPresent)
            topAnalysesInd <- topAnalysesInd[gTable[[fg]][topAnalysesInd] != 0]

        rtMins <- sapply(topAnalysesInd, getFTCol, fg = fg, col = "retmin")
        rtMaxs <- sapply(topAnalysesInd, getFTCol, fg = fg, col = "retmax")
        rtMins[is.na(rtMins)] <- min(rtMins, na.rm = TRUE); rtMaxs[is.na(rtMaxs)] <- max(rtMaxs, na.rm = TRUE)
        rtMins <- rtMins - rtWindow; rtMaxs <- rtMaxs + rtWindow
        
        mzMins <- sapply(topAnalysesInd, getFTCol, fg = fg, col = "mzmin")
        mzMaxs <- sapply(topAnalysesInd, getFTCol, fg = fg, col = "mzmax")
        mzMins[is.na(mzMins)] <- min(mzMins, na.rm = TRUE) - mzExpWindow
        mzMaxs[is.na(mzMaxs)] <- max(mzMaxs, na.rm = TRUE) + mzExpWindow
        
        return(data.table(analysis = anaInfo$analysis[topAnalysesInd],
                          retmin = rtMins, retmax = rtMaxs,
                          mzmin = mzMins, mzmax = mzMaxs))
    }, simplify = FALSE), idcol = "group")
    
    cacheDB <- openCacheDBScope()

    EICs <- Map(anaInfo$analysis, anaInfo$path, f = function(ana, path)
    {
        EICInfoAna <- EICInfo[analysis == ana]
        if (nrow(EICInfoAna) == 0)
            return(NULL)
        
        dfile <- getMzMLOrMzXMLAnalysisPath(ana, path)
        anaHash <- makeFileHash(dfile)
        
        EICInfoAna[, hash := makeHash(anaHash, .SD), by = seq_len(nrow(EICInfoAna)),
                   .SDcols = c("retmin", "retmax", "mzmin", "mzmax")]

        gEICs <- pruneList(setNames(lapply(EICInfoAna$hash, loadCacheData, category = "mzREIC", dbArg = cacheDB),
                                    EICInfoAna$group))

        EICInfoAnaTODO <- EICInfoAna[!group %in% names(gEICs)]
        if (nrow(EICInfoAnaTODO) > 0)
        {
            spectra <- loadSpectra(dfile, verbose = FALSE, cacheDB = cacheDB)
            eics <- setNames(loadEICs(spectra, EICInfoAnaTODO$retmin, EICInfoAnaTODO$retmax,
                                      EICInfoAnaTODO$mzmin, EICInfoAnaTODO$mzmax), EICInfoAnaTODO$group)

            for (i in seq_along(eics))
                saveCacheData("mzREIC", eics[[i]], EICInfoAnaTODO$hash[i], cacheDB)
            
            gEICs <- c(gEICs, eics)
        }

        return(gEICs)
    })
    
    return(pruneList(EICs))
})

setMethod("getEICsForFGroups", "featureGroupsSet", function(fGroups, rtWindow, mzExpWindow, topMost, topMostByRGroup,
                                                            onlyPresent)
{
    ionizedFGroupsList <- sapply(sets(fGroups), ionize, obj = fGroups, simplify = FALSE)
    EICList <- sapply(ionizedFGroupsList, getEICsForFGroups, rtWindow = rtWindow, mzExpWindow = mzExpWindow,
                      topMost = topMost, topMostByRGroup = topMostByRGroup, onlyPresent = onlyPresent, simplify = FALSE)
    EICs <- unlist(EICList, recursive = FALSE, use.names = FALSE) # use.names gives combined set/ana name, we just want ana
    names(EICs) <- unlist(lapply(EICList, names))
    EICs <- EICs[intersect(analyses(fGroups), names(EICs))] # sync order
    
    return(EICs)
})

getEICsForFeatures <- function(features)
{
    if (length(features) == 0)
        return(list())
    
    fTable <- featureTable(features)
    anaInfo <- analysisInfo(features)
    
    cacheDB <- openCacheDBScope()
    
    EICs <- Map(anaInfo$analysis, anaInfo$path, fTable, f = function(ana, path, ft)
    {
        dfile <- getMzMLOrMzXMLAnalysisPath(ana, path)
        anaHash <- makeFileHash(dfile)
        
        hashes <- ft[, makeHash(anaHash, .SD), by = seq_len(nrow(ft)),
                     .SDcols = c("retmin", "retmax", "mzmin", "mzmax")][[2]]
        
        fEICs <- lapply(hashes, loadCacheData, category = "mzREIC", dbArg = cacheDB)
        isNotCached <- sapply(fEICs, is.null)
        
        if (any(isNotCached) > 0)
        {
            spectra <- loadSpectra(dfile, verbose = FALSE, cacheDB = cacheDB)
            ftTODO <- ft[isNotCached]
            fEICs[isNotCached] <- loadEICs(spectra, ftTODO$retmin, ftTODO$retmax,
                                           ftTODO$mzmin, ftTODO$mzmax)
            
            for (i in which(isNotCached))
                saveCacheData("mzREIC", fEICs[[i]], hashes[i], cacheDB)
        }
        
        return(fEICs)
    })
    
    return(pruneList(EICs))
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
