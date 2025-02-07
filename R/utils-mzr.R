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
