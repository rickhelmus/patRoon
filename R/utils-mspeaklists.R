
#' @details The \code{getDefAvgPListParams} is used to create a parameter list
#'   for peak list averaging (discussed below).
#' @rdname MSPeakLists-generation
#' @export
getDefAvgPListParams <- function(...)
{
    def <- list(clusterMzWindow = 0.005,
                topMost = 50,
                minIntensityPre = 500,
                minIntensityPost = 500,
                avgFun = mean,
                method = "hclust")
    return(modifyList(def, list(...)))
}

# align & average spectra by clustering or between peak distances
# code inspired from msProcess R package: https://github.com/zeehio/msProcess
averageSpectra <- function(spectra, clusterMzWindow, topMost, minIntensityPre, minIntensityPost, avgFun, method)
{
    if (length(spectra) == 0) # no spectra, return empty spectrum
        return(data.table(mz = integer(0), intensity = integer(0)))
    
    spectra <- lapply(spectra, function(s)
    {
        s <- s[intensity >= minIntensityPre]
        if (nrow(s) > topMost)
        {
            ord <- order(-s$intensity)
            s <- s[ord[seq_len(topMost)]]
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
        hc <- fastcluster::hclust(mzd)
        spcomb[, cluster := cutree(hc, h = clusterMzWindow)]
    }
    else if (method == "distance")
    {
        mzdiff <- abs(diff(spcomb$mz))
        spcomb[, cluster := 1 + c(0, cumsum(mzdiff > clusterMzWindow))]
    }
    
    if (any(spcomb[, .(dup = anyDuplicated(spid)), key = cluster][["dup"]] > 0))
        warning("During spectral averaging multiple masses from the same spectrum were clustered, consider tweaking clusterMzWindow!\n")
    
    spcount <- length(spectra)
    ret <- spcomb[, .(mz = avgFun(mz), intensity = sum(intensity) / spcount), by = cluster][, cluster := NULL]
    
    # post intensity filter
    ret <- ret[intensity >= minIntensityPost]
    
    return(ret)
}

averageMSPeakLists <- function(peakLists, clusterMzWindow, topMost, minIntensityPre, minIntensityPost, avgFun, method)
{
    # UNDONE: use cache sets?
    
    cat("Generating averaged peak lists for all feature groups...\n")
    
    hash <- makeHash(peakLists, clusterMzWindow, topMost, minIntensityPre, avgFun, method)
    avgPLists <- loadCacheData("MSPeakListsAvg", hash)
    
    # figure out feature groups from (non-averaged) peak lists
    gNames <- unique(unlist(sapply(peakLists, names, simplify = FALSE), use.names = FALSE))
    gCount <- length(gNames)
    
    if (gCount == 0)
        avgPLists <- list()
    else if (is.null(avgPLists))
    {
        prog <- txtProgressBar(0, gCount, style = 3)
        
        avgPLists <- lapply(seq_len(gCount), function(grpi)
        {
            plistsMS <- lapply(peakLists, function(pl) pl[[gNames[grpi]]][["MS"]])
            plistsMS <- plistsMS[!sapply(plistsMS, is.null)]
            plistsMS <- averageSpectra(plistsMS, clusterMzWindow, topMost, minIntensityPre, minIntensityPost, avgFun, method)
            
            plistsMSMS <- lapply(peakLists, function(pl) pl[[gNames[grpi]]][["MSMS"]])
            plistsMSMS <- plistsMSMS[!sapply(plistsMSMS, is.null)]
            plistsMSMS <- averageSpectra(plistsMSMS, clusterMzWindow, topMost, minIntensityPre, minIntensityPost, avgFun, method)
            
            setTxtProgressBar(prog, grpi)
            return(pruneList(list(MS = if (nrow(plistsMS) > 0) plistsMS else NULL,
                                  MSMS = if (nrow(plistsMSMS) > 0) plistsMSMS else NULL)))
        })
        names(avgPLists) <- gNames
        
        setTxtProgressBar(prog, gCount)
        close(prog)
        
        saveCacheData("MSPeakListsAvg", avgPLists, hash)
    }
    else
        cat("Done!\n")
    
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

doMSPeakListFilter <- function(pList, absMSIntThr, absMSMSIntThr, relMSIntThr,
                               relMSMSIntThr, topMSPeaks, topMSMSPeaks,
                               deIsotopeMS, deIsotopeMSMS)
{
    if (!is.null(pList[["MS"]]))
    {
        if (!is.null(absMSIntThr))
            pList[["MS"]] <- pList[["MS"]][intensity >= absMSIntThr]
        
        if (!is.null(relMSIntThr) && nrow(pList[["MS"]]) > 0)
        {
            thr <- max(pList[["MS"]]$intensity) * relMSIntThr
            pList[["MS"]] <- pList[["MS"]][intensity >= thr]
        }
        
        if (!is.null(topMSPeaks) && nrow(pList[["MS"]]) > topMSPeaks)
        {
            ord <- order(-pList[["MS"]]$intensity)
            pList[["MS"]] <- pList[["MS"]][ord[seq_len(topMSPeaks)]]
        }
        
        if (deIsotopeMS)
            pList[["MS"]] <- deIsotopeMSPeakList(pList[["MS"]])
    }
    
    if (!is.null(pList[["MSMS"]]))
    {
        if (!is.null(absMSMSIntThr) && !is.null(pList[["MSMS"]]))
            pList[["MSMS"]] <- pList[["MSMS"]][intensity >= absMSMSIntThr]
        
        if (!is.null(relMSMSIntThr) && nrow(pList[["MS"]]) > 0)
        {
            thr <- max(pList[["MSMS"]]$intensity) * relMSMSIntThr
            pList[["MSMS"]] <- pList[["MSMS"]][intensity >= thr]
        }
        
        if (!is.null(topMSMSPeaks) && nrow(pList[["MSMS"]]) > topMSMSPeaks)
        {
            ord <- order(-pList[["MSMS"]]$intensity)
            pList[["MSMS"]] <- pList[["MSMS"]][ord[seq_len(topMSMSPeaks)]]
        }
        
        if (deIsotopeMSMS)
            pList[["MSMS"]] <- deIsotopeMSPeakList(pList[["MSMS"]])
    }
    
    # prune empty
    if (!is.null(pList[["MS"]]) && nrow(pList[["MS"]]) == 0)
        pList[["MS"]] <- NULL
    if (!is.null(pList[["MSMS"]]) && nrow(pList[["MSMS"]]) == 0)
        pList[["MSMS"]] <- NULL
    
    return(pList)
}
