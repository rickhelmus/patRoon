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

setMethod("getEICFGroupInfo", "featureGroups", function(fGroups, analysis, groupName, EICParams)
{
    takeAnalysis <- analysis # copy name to workaround DT access below
    
    anaInfo <- analysisInfo(fGroups)[analysis %chin% takeAnalysis]
    featTab <- as.data.table(getFeatures(fGroups))
    
    topMost <- if (!is.null(EICParams$topMost))
        min(EICParams$topMost, nrow(anaInfo))
    else
        NULL
    
    # subset relevant things in advance    
    featTab <- featTab[group %chin% groupName, c("group", "analysis", "intensity", "retmin", "retmax", "mzmin", "mzmax"),
                       with = FALSE]
    
    # NOTE: we subset and split here in advance, as doing it in the loop below gets quite slow with many fGroups
    featTabAnaSub <- featTab[analysis %chin% takeAnalysis]
    featTabSplitGrp <- split(featTab, by = "group", keep.by = FALSE)
    featTabAnaSubSplitGrp <- split(featTabAnaSub, by = "group", keep.by = FALSE)
    
    return(sapply(groupName, function(fg)
    {
        ret <- featTabAnaSubSplitGrp[[fg]]
        if (is.null(ret))
            ret <- featTabSplitGrp[[fg]][0] # not present for this analysis, take full table to get all columns
        ret <- copy(ret)
        
        # add missing analyses if needed
        if (!EICParams$onlyPresent)
        {
            if (any(!analysis %chin% ret$analysis))
            {
                ftAllAna <- featTabSplitGrp[[fg]]
                ret <- rbind(ret, data.table(analysis = setdiff(analysis, ret$analysis), intensity = 0,
                                             retmin = min(ftAllAna$retmin), retmax = max(ftAllAna$retmax),
                                             mzmin = min(ftAllAna$mzmin) - EICParams$mzExpWindow,
                                             mzmax = max(ftAllAna$mzmax) + EICParams$mzExpWindow))
            }
        }
        
        if (!is.null(topMost))
        {
            if (EICParams$topMostByRGroup)
            {
                ret[, rGroup := anaInfo$group[match(analysis, anaInfo$analysis)]]
                ret[, rank := frank(-intensity, ties.method = "first"), by = "rGroup"]
                ret <- ret[rank <= topMost]
            }
            else
            {
                setorderv(ret, "intensity", order = -1L)
                ret <- ret[seq_len(topMost)]
            }
        }
        
        ret[, c("retmin", "retmax") := .(retmin - EICParams$rtWindow, retmax + EICParams$rtWindow)]
        return(ret)
    }, simplify = FALSE))
})

setMethod("getEICFGroupInfo", "featureGroupsSet", function(fGroups, analysis, groupName, EICParams)
{
    ret <- callNextMethod()
    
    anaInfo <- analysisInfo(fGroups)
    featTab <- as.data.table(getFeatures(fGroups))
    
    # HACK: since feature tables store the character form, it's easier to keep it all the same
    EICParams$setsAdductPos <- as.character(EICParams$setsAdductPos)
    EICParams$setsAdductNeg <- as.character(EICParams$setsAdductNeg)
    
    # 'ionize' m/zs
    return(Map(names(ret), ret, f = function(grp, ranges)
    {
        featTabGrp <- featTab[group == grp]
        ranges[, adduct := featTabGrp[match(ranges$analysis, analysis)]$adduct]
        
        if (!EICParams$onlyPresent && any(is.na(ranges$adduct))) # adduct will be NA for 'missing' features
        {
            # First try to get adduct from other features in the same set: assume that adduct per set for a single
            # feature group is always the same
            adductSets <- unique(featTabGrp[, c("adduct", "set"), with = FALSE])
            ranges[is.na(adduct), adduct := {
                s <- anaInfo$set[match(analysis, anaInfo$analysis)]
                adductSets[set == s]$adduct
            }, by = "analysis"]
            
            # Then fallback to default adducts for a polarity. For this we get the polarity from another feature in the
            # same analysis (if present)
            ranges[is.na(adduct), adduct := sapply(analysis, function(ana)
            {
                t <- featTab[analysis == ana]
                if (nrow(t) == 0)
                    NA # all features were removed
                else if (as.adduct(t$adduct[1])@charge > 0)
                    EICParams$setsAdductPos
                else
                    EICParams$setsAdductNeg
            })]
            
            if (any(is.na(ranges$adduct)))
            {
                # If all failed then simply omit
                warning(sprintf("Cannot get adduct information for group '%s' for features in analyses %s", grp,
                                paste0(ranges[is.na(adduct)]$analysis, collapse = ", ")), call. = FALSE)
                ranges <- ranges[!is.na(adduct)]
            }
        }

        # NOTE: this is mostly copied from unset features method        
        allAdducts <- sapply(unique(ranges$adduct), as.adduct)
        mzmins <- calculateMasses(ranges$mzmin, allAdducts[ranges$adduct], type = "mz")
        nmd <- mzmins - ranges$mzmin
        set(ranges, j = c("mzmin", "mzmax"), value = list(mzmins, ranges$mzmax + nmd))
        
        return(ranges)
    }))
})

setMethod("getEICsForFGroups", "featureGroups", function(fGroups, analysis, groupName, EICParams)
{
    if (length(fGroups) == 0 || length(analysis) == 0 || length(groupName) == 0)
        return(list())
    
    takeAnalysis <- analysis # for DT subset below
    anaInfo <- analysisInfo(fGroups)[analysis %chin% takeAnalysis]

    EICInfoTab <- getEICFGroupInfo(fGroups, analysis, groupName, EICParams)
    EICInfo <- split(rbindlist(EICInfoTab, idcol = "group"), by = "analysis")
    EICInfo <- EICInfo[intersect(anaInfo$analysis, names(EICInfo))] # sync order
    anaInfoEICs <- anaInfo[analysis %in% names(EICInfo)]

    EICs <- doGetEICs(anaInfoEICs, EICInfo)
    EICs <- Map(EICs, lapply(EICInfo, "[[", "group"), f = setNames)
    
    return(pruneList(EICs))
})

setMethod("getEICsForFeatures", "features", function(features)
{
    if (length(features) == 0)
        return(list())
    
    fTable <- featureTable(features)
    anaInfo <- analysisInfo(features)
    
    return(pruneList(doGetEICs(anaInfo, fTable)))
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
