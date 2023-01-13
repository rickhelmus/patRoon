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

doGetEICs <- function(file, ranges, cacheDB = NULL)
{
    anaHash <- makeFileHash(file)
    
    # NOTE: subset columns here, so any additional columns from e.g. feature tables are not considered
    hashes <- ranges[, makeHash(anaHash, .SD), by = seq_len(nrow(ranges)),
                     .SDcols = c("retmin", "retmax", "mzmin", "mzmax")][[2]]
    
    EICs <- lapply(hashes, loadCacheData, category = "mzREIC", dbArg = cacheDB)
    isNotCached <- sapply(EICs, is.null)
    
    if (any(isNotCached))
    {
        spectra <- loadSpectra(file, verbose = FALSE, cacheDB = cacheDB)
        rangesToDo <- ranges[isNotCached]
        EICs[isNotCached] <- loadEICs(spectra, rangesToDo$retmin, rangesToDo$retmax,
                                      rangesToDo$mzmin, rangesToDo$mzmax)
        
        for (i in which(isNotCached))
            saveCacheData("mzREIC", EICs[[i]], hashes[i], cacheDB)
    }
    
    return(EICs)
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

    verifyDataCentroided(anaInfo)
    
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
    
    EICInfoTab <- sapply(gNames, function(fg)
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
    }, simplify = FALSE)
    
    cacheDB <- openCacheDBScope()
    
    EICInfo <- split(rbindlist(EICInfoTab, idcol = "group"), by = "analysis")
    EICInfo <- EICInfo[intersect(anaInfo$analysis, names(EICInfo))] # sync order
    anaInfoEICs <- anaInfo[anaInfo$analysis %in% names(EICInfo), ]
    anaPaths <- mapply(anaInfoEICs$analysis, anaInfoEICs$path, FUN = getMzMLOrMzXMLAnalysisPath,
                       MoreArgs = list(mustExist = TRUE))
    
    EICs <- Map(anaPaths, EICInfo, f = doGetEICs, MoreArgs = list(cacheDB = cacheDB))
    EICs <- Map(EICs, lapply(EICInfo, "[[", "group"), f = setNames)
    
    return(pruneList(EICs))
})

setMethod("getEICsForFGroups", "featureGroupsSet", function(fGroups, rtWindow, mzExpWindow, topMost, topMostByRGroup,
                                                            onlyPresent)
{
    unsetFGroupsList <- sapply(sets(fGroups), unset, obj = fGroups, simplify = FALSE)
    EICList <- sapply(unsetFGroupsList, getEICsForFGroups, rtWindow = rtWindow, mzExpWindow = mzExpWindow,
                      topMost = topMost, topMostByRGroup = topMostByRGroup, onlyPresent = onlyPresent, simplify = FALSE)
    EICs <- unlist(EICList, recursive = FALSE, use.names = FALSE) # use.names gives combined set/ana name, we just want ana
    names(EICs) <- unlist(lapply(EICList, names))
    EICs <- EICs[intersect(analyses(fGroups), names(EICs))] # sync order
    
    if (!is.null(topMost))
    {
        # topMost is applied per set, make sure that the final result also doesn't contain >topMost results
        
        gTable <- groupTable(fGroups)
        gNames <- names(fGroups)
        anaInfo <- analysisInfo(fGroups)
        anasInEICs <- names(EICs)
        anaIndsInEICs <- match(anasInEICs, anaInfo$analysis)
        topMost <- min(topMost, length(anasInEICs))
        
        for (fg in gNames)
        {
            if (topMostByRGroup)
            {
                tbl <- data.table(int = gTable[[fg]], group = anaInfo$group, ana = anaInfo$analysis)
                tbl <- tbl[ana %in% anasInEICs]
                tbl[, rank := frank(-int, ties.method = "first"), by = "group"]
                topAnalyses <- tbl[rank <= topMost]$ana
            }
            else
            {
                oint <- order(gTable[[fg]][anaIndsInEICs], decreasing = TRUE)
                topAnalyses <- anasInEICs[oint[seq_len(topMost)]]
            }
            
            # clearout any analysis results not being in topMost
            otherAnas <- setdiff(anasInEICs, topAnalyses)
            EICs[otherAnas] <- lapply(EICs[otherAnas], function(e) e[setdiff(names(e), fg)])
        }
        EICs <- pruneList(EICs, checkEmptyElements = TRUE) # in case all EICs were removed from analyses
    }

    return(EICs)
})

setMethod("getEICsForFeatures", "features", function(features)
{
    if (length(features) == 0)
        return(list())
    
    fTable <- featureTable(features)
    anaInfo <- analysisInfo(features)
    
    verifyDataCentroided(anaInfo)
    
    cacheDB <- openCacheDBScope()
    anaPaths <- mapply(anaInfo$analysis, anaInfo$path, FUN = getMzMLOrMzXMLAnalysisPath,
                       MoreArgs = list(mustExist = TRUE))
    EICs <- Map(anaPaths, fTable, f = doGetEICs, MoreArgs = list(cacheDB = cacheDB))
    
    return(pruneList(EICs))
})

setMethod("getEICsForFeatures", "featuresSet", function(features)
{
    unsetFeatList <- sapply(sets(features), unset, obj = features, simplify = FALSE)
    EICList <- sapply(unsetFeatList, getEICsForFeatures, simplify = FALSE)
    EICs <- unlist(EICList, recursive = FALSE, use.names = FALSE) # use.names gives combined set/ana name, we just want ana
    names(EICs) <- unlist(lapply(EICList, names))
    EICs <- EICs[intersect(analyses(features), names(EICs))] # sync order
    return(EICs)
})

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

verifyDataCentroided <- function(anaInfo)
{
    cacheDB <- openCacheDB()
    
    printf("Verifying if your data is centroided...\n")
    
    isCentroided <- mapply(anaInfo$analysis, anaInfo$path, FUN = function(ana, path)
    {
        fpath <- getMzMLOrMzXMLAnalysisPath(ana, path, mustExist = TRUE)
        
        # hash <- makeFileHash(fpath) complete file hash is most accurate, but slow with many analyses. For now assume
        # user keeps centroided data centroided, and only check again if it wasn't
        hash <- makeHash(fpath)
        cd <- loadCacheData("dataCentroided", hash, cacheDB)
        if (!is.null(cd) && cd)
            return(TRUE)

        msf <- mzR::openMSfile(fpath)
        # NOTE: don't check more than first 100 spectra: most often the first will tell us enough, and loading
        # everything takes some time, _especially_ if it profile data.
        hd <- mzR::header(msf, seq_len(min(100, length(msf))))
        mzR::close(msf)
        
        # UNDONE: for now we just don't put out any warnings if there is no centroided flag available
        isCentr <- is.null(hd[["centroided"]]) || all(hd$centroided)
        saveCacheData("dataCentroided", isCentr, hash, cacheDB)
        
        return(isCentr)
    })
    
    if (!all(isCentroided))
        warning("It seems that the following files may not be centroided: ",
                paste0(anaInfo$analysis[!isCentroided], collapse = ", "),
                ". Please ensure that your MS data is centroided, for instance by using convertMSFiles()",
                call. = FALSE)
    
    invisible(NULL)
}
