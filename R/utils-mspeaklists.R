emptyMSPeakList <- function() data.table(mz = numeric(), intensity = numeric(), precursor = logical())


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
                method = "hclust",
                pruneMissingPrecursorMS = TRUE,
                retainPrecursorMSMS = TRUE)
    return(modifyList(def, list(...)))
}

# For docs
getDefAvgPListParamsRD <- function()
{
    def <- getDefAvgPListParams()
    def <- sapply(def, function(v) if (is.character(v)) paste0("\"", v, "\"") else v)
    def$avgFun <- "mean" # UNDONE?
    return(paste0("\\code{", names(def), "=", def, "}", collapse = "; "))
}

#' @details The \code{getDefIsolatePrecParams} is used to create a parameter
#'   list for isolating the precursor and its isotopes (see \verb{Isolating precursor data}).
#' @rdname MSPeakLists-class
#' @export
getDefIsolatePrecParams <- function(...)
{
    def <- list(maxIsotopes = 5,
                mzDefectRange = c(-0.01, 0.01),
                intRange = c(0.001, 2),
                z = 1,
                maxGap = 2)
    return(modifyList(def, list(...)))
}

# For docs
getDefIsolatePrecParamsRD <- function()
{
    def <- getDefIsolatePrecParams()
    def <- sapply(def, function(v) if (length(v) == 2) sprintf("c(%s, %s)", v[1], v[2]) else as.character(v))
    return(paste0("\\code{", names(def), "=", def, "}", collapse = "; "))
}

# align & average spectra by clustering or between peak distances
# code inspired from msProcess R package: https://github.com/zeehio/msProcess
averageSpectra <- function(spectra, clusterMzWindow, topMost, minIntensityPre, minIntensityPost,
                           avgFun, method, pruneMissingPrecursor, retainPrecursor)
{
    if (length(spectra) == 0) # no spectra, return empty spectrum
        return(emptyMSPeakList())

    spectra <- lapply(spectra, function(s)
    {
        s <- s[intensity >= minIntensityPre]

        if (nrow(s) > topMost)
        {
            ord <- order(s$intensity, decreasing = TRUE)
            keep <- ord[seq_len(topMost)]
            if (retainPrecursor)
                keep <- unique(c(keep, s[precursor == TRUE, which = TRUE]))
            s <- s[keep]
        }
        return(s)
    })

    spcomb <- rbindlist(spectra, idcol = "spid")
    setorder(spcomb, mz)

    if (nrow(spcomb) < 2 || length(spectra) < 2)
        return(spcomb[, -"spid"])

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

    # update precursor flags
    precMZs <- spcomb[precursor == TRUE, mz]
    if (length(precMZs) > 0)
        precMZ <- mean(precMZs)
    else
        precMZ <- -1
    ret <- assignPrecursorToMSPeakList(ret, precMZ)

    # post intensity filter
    ret <- ret[intensity >= minIntensityPost]

    if (pruneMissingPrecursor && !any(ret$precursor))
        return(emptyMSPeakList())

    return(ret)
}

# get corresponding mz of feature from MS peaklist
getMZIndexFromMSPeakList <- function(featMZ, plist)
{
    ret <- which.min(abs(plist$mz - featMZ))
    if (abs(plist$mz[ret] - featMZ) <= 0.01) # UNDONE: make range configurable?
        return(ret)
    return(NA)
}

assignPrecursorToMSPeakList <- function(MSPeakList, precursorMZ)
{
    if (nrow(MSPeakList) == 0)
        MSPeakList[, precursor := logical()]
    else
    {
        MSPeakList[, precursor := FALSE]
        ind <- getMZIndexFromMSPeakList(precursorMZ, MSPeakList)
        if (!is.na(ind))
            MSPeakList[ind, precursor := TRUE]
    }
    return(MSPeakList)
}

deIsotopeMSPeakList <- function(MSPeakList, negate)
{
    if (nrow(MSPeakList) == 0)
        return(MSPeakList)

    if (is.null(MSPeakList[["cmp"]]))
        stop("No isotope information available. Note that this is currently only implemented for DataAnalysis peak lists (if configured properly, see ?generateMSPeakLists.")

    # make sure most intense ions top the table
    MSPeakList <- MSPeakList[order(mz, -intensity)]

    unique_iso <- sapply(seq_along(MSPeakList$cmp), function(i)
    {
        # first and unassigned compounds are always unique
        if (i == 1 || !nzchar(MSPeakList$cmp[i]))
            return(TRUE)

        # peak may belong to multiple isotope compounds (separated by whitespace)
        cmp <- strsplit(MSPeakList$cmp[i], "\\s+")

        # isotope compound present before this entry?
        othercmp <- MSPeakList[seq_len(i - 1)][nzchar(cmp)]$cmp
        for (ocmp in othercmp)
        {
            if (any(cmp %in% strsplit(ocmp, "\\s+")))
                return(FALSE)
        }

        return(TRUE)
    }, USE.NAMES = FALSE)

    if (negate)
        unique_iso <- !unique_iso

    return(MSPeakList[unique_iso])
}

doMSPeakListFilter <- function(pList, absIntThr, relIntThr, topMost, minPeaks, deIsotope, retainPrecursor, negate)
{
    ret <- pList

    intPred <- if (negate) function(i, t) i < t else function(i, t) i >= t

    if (!is.null(absIntThr))
        ret <- ret[intPred(intensity, absIntThr)]

    if (!is.null(relIntThr) && nrow(ret) > 0)
    {
        thr <- max(ret$intensity) * relIntThr
        ret <- ret[intPred(intensity, thr)]
    }

    if (!is.null(topMost) && nrow(ret) > topMost)
    {
        ord <- order(ret$intensity, decreasing = !negate)
        ret <- ret[seq_len(.N) %in% ord[seq_len(topMost)]] # NOTE: keep order
    }

    if (deIsotope)
        ret <- deIsotopeMSPeakList(ret, negate)

    # re-add precursor if necessary
    if (retainPrecursor && !any(ret$precursor))
    {
        prec <- pList[precursor == TRUE]
        if (nrow(prec) > 0)
        {
            ret <- rbind(ret, prec)
            setorderv(ret, "mz")
        }
    }

    if (!is.null(minPeaks))
    {
        notEnough <- (nrow(ret)-1) < minPeaks # -1: don't count precursor peaks
        if (negate != notEnough)
            ret <- ret[0, ]
    }
    
    return(ret)
}

isolatePrecInMSPeakList <- function(plist, isolatePrec, negate)
{
    prec <- plist[precursor == TRUE]
    if (nrow(prec) == 1) # 0 if no precursor
    {
        s <- seq_len(isolatePrec$maxIsotopes) / isolatePrec$z
        mzranges <- matrix(c(s + isolatePrec$mzDefectRange[1] * s,
                             s + isolatePrec$mzDefectRange[2] * s), ncol = 2) + prec$mz
        keep <- plist[, precursor | (
            # only keep peaks with reasonable intensities
            intensity %between% (isolatePrec$intRange * prec$intensity) &

            # only keep with reasonably close m/z to precursor, taking
            # in to account larger windows for larger distances
            inrange(mz, mzranges[, 1], mzranges[, 2]))]
        plist <- plist[if (negate) !keep else keep]

        # remove gaps
        gaps <- plist[round(mz - shift(mz)) > (isolatePrec$maxGap / isolatePrec$z), which = TRUE]
        if (length(gaps) > 0)
        {
            sq <- seq_len(gaps[1] - 1)
            plist <- plist[if (negate) -sq else sq]
        }
    }

    return(plist)
}

getSpec <- function(MSPeakLists, groupName, MSLevel, analysis)
{
    MSInd <- if (MSLevel == 1) "MS" else "MSMS"
    if (!is.null(analysis))
        spec <- MSPeakLists[[analysis, groupName]][[MSInd]]
    else
        spec <- MSPeakLists[[groupName]][[MSInd]]
    return(spec)
}

getMSPeakListPlotTitle <- function(MSLevel, analysis, groupName)
{
    MSInd <- if (MSLevel == 1) "MS" else "MSMS"
    if (!is.null(analysis))
        return(sprintf("%s (%s) %s", groupName, analysis, MSInd))
    return(paste(groupName, MSInd))
}

# NOTE: This function is mainly for debugging the C++ version
binPeakLists <- function(pl1, pl2, shift, absMzDev)
{
    prep <- function(pl)
    {
        pl <- copy(pl)
        pl[, c("index", "low", "high") := .(seq_len(nrow(pl)), mz - absMzDev, mz + absMzDev)]
        setkeyv(pl, c("low", "high"))
        return(pl)
    }
    pl1 <- prep(pl1); pl2 <- prep(pl2)
    
    if (shift != "none")
    {
        if (!any(pl1$precursor) || !any(pl2$precursor))
            stop("Cannot shift spectra: one or both lack precursor ion!")
        precDiff <- pl2[precursor == TRUE]$mz - pl1[precursor == TRUE]$mz
        if (shift == "precursor")
            pl2[, mz := mz - precDiff]
        else # both
        {
            # first bin as normal (recursive call)
            binNone <- binPeakLists(pl1, pl2, "none", absMzDev)
            pl1Unique <- setnames(binNone[intensity_2 == 0, -"intensity_2"], "intensity_1", "intensity")
            pl2Unique <- setnames(binNone[intensity_1 == 0, -"intensity_1"], "intensity_2", "intensity")
            
            # bin missing with shift
            pl2Unique[, mz := mz - precDiff]
            binShift <- binPeakLists(pl1Unique, pl2Unique, "none", absMzDev)
            
            # merge both: add missing from binNone
            ret <- rbind(binNone[intensity_1 != 0 & intensity_2 != 0], binShift)
            setorderv(ret, "mz")
            return(ret)
        }
    }
    
    ov <- foverlaps(pl1, pl2)
    # NOTE: i.xx cols are from left
    
    ret <- ov[, c("i.mz", "intensity", "i.intensity"), with = FALSE]
    setnames(ret, c("i.mz", "intensity", "i.intensity"), c("mz", "intensity_2", "intensity_1"))

    # all left entries should be in there, right can miss and should be added manually
    stopifnot(all(pl1$index %in% ov$i.index))
    ret <- rbind(ret, setnames(pl2[!index %in% ov$index, c("mz", "intensity")], "intensity", "intensity_2"), fill = TRUE)
    
    setnafill(ret, fill = 0)
    setorderv(ret, "mz")
    setcolorder(ret, c("mz", "intensity_1", "intensity_2"))
    
    return(ret)
}

# NOTE: This function is mainly for debugging the C++ version
specSimilarityR <- function(pl1, pl2, method, shift = "none", removePrecursor = FALSE, mzWeight = 0, intWeight = 1,
                           absMzDev = 0.005, relMinIntensity = 0.1)
{
    # UNDONE: refs
    
    # code contributed by Bas van de Velde
    # cosine similarity from OrgMassSpecR
    
    binnedPL <- as.data.table(binSpecCPP(pl1, pl2, shift, absMzDev))
    # binnedPL <- binPeakLists(pl1, pl2, shift, absMzDev)
    
    # remove precursor
    if (removePrecursor)
        binnedPL <- binnedPL[precursor == FALSE]
    
    # normalize
    binnedPL[, c("intensity_1", "intensity_2") :=
                 .(normalize(intensity_1, FALSE), normalize(intensity_2, FALSE))]
    
    # intensity filter
    binnedPL[intensity_1 < relMinIntensity, intensity_1 := 0]
    binnedPL[intensity_2 < relMinIntensity, intensity_2 := 0]
    binnedPL <- binnedPL[intensity_1 != 0 | intensity_2 != 0]
    
    u <- binnedPL$mz^mzWeight * binnedPL$intensity_1^intWeight
    v <- binnedPL$mz^mzWeight * binnedPL$intensity_2^intWeight

    return(switch(method,
                  cosine = as.vector((u %*% v) / (sqrt(sum(u^2)) * sqrt(sum(v^2)))),
                  pearson = cor(u, v, method = "pearson"),
                  spearman = cor(u, v, method = "spearman"),
                  jaccard = binnedPL[intensity_1 != 0 & intensity_2 != 0, .N] / nrow(binnedPL)))
}

prepSpecSimilarity <- function(pl1, pl2, removePrecursor, relMinIntensity, minPeaks)
{
    if (removePrecursor)
    {
        pl1 <- pl1[precursor == FALSE]
        pl2 <- pl2[precursor == FALSE]
    }
    
    if (relMinIntensity > 0)
    {
        filterPL <- function(pl)
        {
            minInt <- relMinIntensity * max(pl$intensity)
            return(pl[numGTE(intensity, minInt)])
        }
        pl1 <- filterPL(pl1); pl2 <- filterPL(pl2)
    }
    
    if (nrow(pl1) == 0 || nrow(pl2) == 0)
        return(list())
    
    if (nrow(pl1) < minPeaks || nrow(pl2) < minPeaks)
        return(list())
    
    return(list(pl1 = pl1, pl2 = pl2))
}

specSimilarity <- function(pl1, pl2, method, shift = "none", precDiff = 0, removePrecursor = FALSE,
                           mzWeight = 0, intWeight = 1, absMzDev = 0.005, relMinIntensity = 0.1,
                           minPeaks = 0)
{
    # UNDONE: refs
    # UNDONE: export? add asserts if yes, otherwise remove defaults
    
    # code contributed by Bas van de Velde
    # cosine similarity from OrgMassSpecR

    prep <- prepSpecSimilarity(pl1, pl2, removePrecursor, relMinIntensity, minPeaks)
    if (length(prep) == 0)
        return(0)
    
    return(calcSpecSimilarity(prep$pl1, prep$pl2, method, shift, precDiff, mzWeight, intWeight, absMzDev))
}

specSimilaritySets <- function(pl1, pl2, method, shift = "none", precDiff = 0, removePrecursor = FALSE,
                               mzWeight = 0, intWeight = 1, absMzDev = 0.005, relMinIntensity = 0.1,
                               minPeaks = 0, combMethod)
{
    # UNDONE: refs
    # UNDONE: export? add asserts if yes
    
    # code contributed by Bas van de Velde
    # cosine similarity from OrgMassSpecR

    if (combMethod == "average")
    {
        allSets <- union(pl1$set, pl2$set)
        sims <- sapply(allSets, function(s)
        {
            spl1 <- pl1[set == s]; spl2 <- pl2[set == s]
            if (nrow(spl1) == 0 || nrow(spl2) == 0)
                return(NA)
            
            prep <- prepSpecSimilarity(spl1, spl2, removePrecursor, relMinIntensity, minPeaks)
            if (length(prep) == 0)
                return(0)
            
            return(calcSpecSimilarity(prep$pl1, prep$pl2, method, shift, precDiff, mzWeight, intWeight, absMzDev))
        })
        if (all(is.na(sims)))
            return(NA_real_)
        return(mean(sims, na.rm = TRUE))
    }
    
    # combMethod == "merge"
    setsBoth <- intersect(pl1$set, pl2$set) # only consider data from sets present in both peaklists
    spl1 <- pl1[set %in% setsBoth]; spl2 <- pl2[set %in% setsBoth]
    if (nrow(spl1) == 0 || nrow(spl2) == 0)
        return(NA)
    
    prep <- prepSpecSimilarity(spl1, spl2, removePrecursor, relMinIntensity, minPeaks)
    if (length(prep) == 0)
        return(0)
    
    return(calcSpecSimilarity(prep$pl1, prep$pl2, method, shift, precDiff, mzWeight, intWeight, absMzDev))
}
