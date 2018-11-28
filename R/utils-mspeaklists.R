# align & average spectra by clustering or between peak distances
# code inspired from msProcess R package: https://github.com/zeehio/msProcess
averageSpectra <- function(spectra, clusterMzWindow, maxPeaks, minIntensity, avgMassFun, method)
{
    if (length(spectra) == 0) # no spectra, return empty spectrum
        return(data.table(mz = integer(0), intensity = integer(0)))
    
    spectra <- lapply(spectra, function(s)
    {
        s <- s[intensity >= minIntensity]
        if (nrow(s) > maxPeaks)
        {
            ord <- order(-s$intensity)
            s <- s[ord[seq_len(maxPeaks)]]
        }
        return(s)
    })
    
    spcomb <- rbindlist(spectra, idcol = "spid")
    setorder(spcomb, mz)
    
    if (nrow(spcomb) < 2 || length(spectra) < 2)
        return(spcomb[, c("mz", "intensity")])
    
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
    
    if (any(spcomb[, .(dup = anyDuplicated(spid)), key = cluster][["dup"]] > 0))
        warning("Clustered multiple masses from same spectrum, consider lowering clusterMzWindow!\n")
    
    spcount <- length(spectra)
    return(spcomb[, .(mz = avgMassFun(mz), intensity = sum(intensity) / spcount), by = cluster][, cluster := NULL])
}

averageMSPeakLists <- function(peakLists, clusterMzWindow, maxPeaks, minIntensity, avgMassFun, method)
{
    # UNDONE: use cache sets?
    
    cat("Generating averaged peak lists for all feature groups...\n")
    
    hash <- makeHash(peakLists, clusterMzWindow, maxPeaks, minIntensity, avgMassFun, method)
    avgPLists <- loadCacheData("MSPeakListsAvg", hash)
    
    if (is.null(avgPLists))
    {
        # figure out feature groups from (non-averaged) peak lists
        gNames <- unique(unlist(sapply(peakLists, names, simplify = FALSE), use.names = FALSE))
        gCount <- length(gNames)
        
        prog <- txtProgressBar(0, gCount, style = 3)
        
        avgPLists <- lapply(seq_len(gCount), function(grpi)
        {
            plistsMS <- lapply(peakLists, function(pl) pl[[gNames[grpi]]][["MS"]])
            plistsMS <- plistsMS[!sapply(plistsMS, is.null)]
            plistsMS <- averageSpectra(plistsMS, clusterMzWindow, maxPeaks, minIntensity, avgMassFun, method)
            
            plistsMSMS <- lapply(peakLists, function(pl) pl[[gNames[grpi]]][["MSMS"]])
            plistsMSMS <- plistsMSMS[!sapply(plistsMSMS, is.null)]
            plistsMSMS <- averageSpectra(plistsMSMS, clusterMzWindow, maxPeaks, minIntensity, avgMassFun, method)
            
            setTxtProgressBar(prog, grpi)
            return(list(MS = if (nrow(plistsMS) > 0) plistsMS else NULL,
                        MSMS = if (nrow(plistsMSMS) > 0) plistsMSMS else NULL))
        })
        names(avgPLists) <- gNames
        
        setTxtProgressBar(prog, gCount)
        close(prog)
        
        saveCacheData("MSPeakListsAvg", avgPLists, hash)
    }
    else
        cat("Done! (cached)")
    
    return(avgPLists)
}

deIsotopeMSPeakList <- function(MSPeakList)
{
    if (nrow(MSPeakList) == 0)
        return(MSPeakList)
    
    # make sure most intense ions top the table
    MSPeakList <- MSPeakList[order(mz, -intensity)]
    
    unique_iso <- sapply(seq_along(MSPeakList$cmp), function(i)
    {
        # first and unassigned compounds are always unique
        if (i == 1 || MSPeakList$cmp[i] == "")
            return(TRUE)
        
        # peak may belong to multiple isotope compounds (separated by whitespace)
        cmp <- strsplit(MSPeakList$cmp[i], "\\s+")
        
        # isotope compound present before this entry?
        othercmp <- MSPeakList[seq_len(i - 1)][cmp != ""]$cmp
        for (ocmp in othercmp)
        {
            if (any(cmp %in% strsplit(ocmp, "\\s+")))
                return(FALSE)
        }
        
        return(TRUE)
    }, USE.NAMES = FALSE)
    
    return(MSPeakList[unique_iso])
}
