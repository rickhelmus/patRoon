# SPDX-FileCopyrightText: 2016 - 2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
NULL

getHeaders <- function(path, rtRange, MSLevel, verbose = FALSE)
{
    if (verbose)
        printf("Loading spectra headers for '%s'...\n", path)
    
    msf <- mzR::openMSfile(path)
    hd <- as.data.table(mzR::header(msf))
    mzR::close(msf)
    
    if (nrow(hd) > 0) {
        if (!is.null(MSLevel))
            hd <- hd[msLevel %in% MSLevel, ]
        
        if (!is.null(rtRange))
            hd <- hd[numGTE(retentionTime, rtRange[1]) & numLTE(retentionTime, rtRange[2]), ]
    }
    
    return(hd)
}

loadSpectra <- function(path, rtRange = NULL, verbose = TRUE, cacheDB = NULL)
{
    # NOTE: limit length as this function may be called frequently
    hash <- makeHash(makeFileHash(path, length = 8192), rtRange)
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

averageSpectraMZR <- function(spectra, hd, clusterMzWindow, topMost, minIntensityPre, minIntensityPost, minAbundance,
                              avgFun, method, precursor, pruneMissingPrecursor, retainPrecursor)
{
    if (nrow(hd) == 0) # no spectra, return empty spectrum
        return(emptyMSPeakList("feat_abundance", NULL))

    sp <- spectra$spectra[hd$seqNum]
    # convert to peaklist format
    sp <- lapply(sp, function(spec) setnames(as.data.table(spec), c("mz", "intensity")))
    sp <- lapply(sp, assignPrecursorToMSPeakList, precursor)

    return(averageSpectra(sp, clusterMzWindow, topMost, minIntensityPre, minIntensityPost, minAbundance,
                          "feat_abundance", avgFun, NULL, method, TRUE, pruneMissingPrecursor, retainPrecursor))
}

verifyDataCentroided <- function(anaInfo)
{
    if (!isTRUE(getOption("patRoon.checkCentroided", default = TRUE)))
        return(invisible(NULL))
    
    cacheDB <- openCacheDB()
    
    printf("Verifying if your data is centroided... ")
    
    filePaths <- getCentroidedMSFilesFromAnaInfo(anaInfo)
    
    isCentroided <- sapply(filePaths, function(fpath)
    {
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
    
    printf("Done!\n")
    invisible(NULL)
}
