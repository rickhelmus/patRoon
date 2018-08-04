#' @include main.R

loadSpectra <- function(path, rtRange = NULL, verbose = TRUE)
{
    hash <- makeHash(makeFileHash(path), rtRange)
    ret <- loadCacheData("specData", hash)
    if (is.null(ret))
    {
        if (verbose)
            printf("Loading raw spectra for '%s'...\n", path)
        msf <- openMSfile(path)
        hd <- as.data.table(header(msf))

        if (is.null(rtRange))
            ps <- peaks(msf) # load all
        else
            ps <- peaks(msf, hd[numGTE(retentionTime, rtRange[1]) & numLTE(retentionTime, rtRange[2]), seqNum])

        spectra <- lapply(ps, function(spec) setnames(as.data.table(spec), c("mz", "intensity")))
        ret <- list(header = hd, spectra = spectra)
        close(msf)
        saveCacheData("specData", ret, hash)
    }

    return(ret)
}

loadAllSpectra <- function(analyses, paths)
{
    ret <- lapply(seq_along(analyses), function(i) loadSpectra(getMzMLOrMzXMLAnalysisPath(analyses[i], paths[i])))
    names(ret) <- analyses
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

# align & average spectra by clustering or between peak distances
# code inspired from msProcess R package: https://github.com/zeehio/msProcess
averageSpectra <- function(spectra, rtRange, clusterMzWindow = 0.003, maxPeaks = 50, minIntensity, MSLevel = 1, precursor = NULL,
                           precursorMzWindow = NULL, avgMassFun = mean, method = "hclust")
{
    hd <- getSpectraHeader(spectra, rtRange, MSLevel, precursor, precursorMzWindow)

    if (nrow(hd) == 0) # no spectra, return empty spectrum
        return(data.table(mz = integer(0), intensity = integer(0)))

    sp <- spectra$spectra[hd$seqNum]
    sp <- lapply(sp, function(s) s[intensity >= minIntensity])

    # limit amount of peaks per spectra, otherwise clustering will clog all memory/cpu
    sp <- lapply(sp, function(s)
    {
        if (nrow(s) > maxPeaks)
        {
            ord <- order(-s$intensity)
            s <- s[ord[seq_len(maxPeaks)]]
        }
        return(s)
    })

    spcomb <- rbindlist(sp)
    spcomb[, ret := rep(hd$retentionTime, sapply(sp, nrow))]
    setorder(spcomb, mz)

    if (nrow(spcomb) < 2)
        return(spcomb)

    # UNDONE: why the hell does this happen??
    if (MSLevel == 2 && !is.null(precursorMzWindow) && max(spcomb$mz) > (precursor + precursorMzWindow))
        warning(sprintf("Found fragment masses above precursor isolation window for m/z %.4f!", precursor))

    if (method == "hclust")
    {
        mzd <- dist(spcomb$mz)
        hc <- hclust(mzd)
        spcomb[, cluster := cutree(hc, h = clusterMzWindow)]

    }
    else if (method == "distance")
    {
        mzdiff <- abs(diff(spcomb$mz))
        spcomb[, cluster := 1 + c(0, cumsum(mzdiff > clusterMzWindow))]
    }

    if (any(spcomb[, .(dup = anyDuplicated(ret)), key = cluster][["dup"]] > 0))
        warning("Clustered multiple masses from same spectrum, consider lowering clusterMzWindow!\n")

    spcount <- nrow(hd)
    return(spcomb[, .(mz = avgMassFun(mz), intensity = sum(intensity) / spcount), by = cluster][, cluster := NULL])
}
