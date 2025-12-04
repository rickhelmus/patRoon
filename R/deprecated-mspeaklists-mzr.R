# SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
#' @include mspeaklists.R
NULL

# align & average spectra by clustering or between peak distances
# code inspired from msProcess R package: https://github.com/zeehio/msProcess
averageSpectra <- function(spectra, clusterMzWindow, topMost, minIntensityPre, minIntensityPost, minAbundanceRel,
                           minAbundanceAbs, abundanceColumn, avgFun, avgCols, method, assignPrecursor,
                           pruneMissingPrecursor, retainPrecursor)
{
    if (length(spectra) == 0) # no spectra, return empty spectrum
        return(emptyMSPeakList(abundanceColumn, avgCols))
    
    abundanceColumnRel <- paste0(abundanceColumn, "_rel"); abundanceColumnAbs <- paste0(abundanceColumn, "_abs")
    
    spectra <- lapply(spectra, function(s)
    {
        s <- s[intensity >= minIntensityPre]
        if (nrow(s) > topMost)
        {
            ord <- order(s$intensity, decreasing = TRUE)
            keep <- ord[seq_len(topMost)]
            if (retainPrecursor)
                keep <- union(keep, which(s$precursor))
            s <- s[keep]
        }
        return(s)
    })
    
    spcount <- length(spectra) # un-filtered count
    spcomb <- rbindlist(spectra, idcol = "spid")
    setorderv(spcomb, "mz")
    
    if (nrow(spcomb) < 2 || length(spectra) < 2)
    {
        if (!is.null(abundanceColumn))
        {
            spcomb[, (abundanceColumnRel) := 1 / spcount]
            spcomb[, (abundanceColumnAbs) := 1]
        }
        return(spcomb[, names(emptyMSPeakList(abundanceColumn, avgCols)), with = FALSE])
    }
    
    if (method == "hclust")
    {
        mzd <- dist(spcomb$mz)
        hc <- fastcluster::hclust(mzd)
        spcomb[, cluster := cutree(hc, h = clusterMzWindow)]
    }
    else if (method == "distance_point")
    {
        mzdiff <- abs(diff(spcomb$mz))
        spcomb[, cluster := 1 + c(0, cumsum(mzdiff > clusterMzWindow))]
    }
    
    if (any(spcomb[, .(dup = anyDuplicated(spid)), key = cluster][["dup"]] > 0))
        warning("During spectral averaging multiple masses from the same spectrum were clustered, consider tweaking clusterMzWindow!\n")
    
    doAvgPL <- function(sps, by, intN)
    {
        if (!is.null(avgCols))
            sps[, (avgCols) := lapply(.SD, mean), .SDcols = avgCols, by = by]
        sps[, c("mz", "intensity") := .(avgFun(mz), sum(intensity) / intN), by = by]
        return(unique(sps, by = by))
    }
    
    # average in two steps:
    # 1. same mass peaks from same analysis (regular intensity average)
    # 2. same mass peaks (intensity average over all considered spectra)
    
    ret <- doAvgPL(spcomb, c("cluster", "spid"), 1) # set intN==1 to sum merged peaks
    
    # do abundance calculation _after_ removing duplicated mass peaks from same spectra
    if (!is.null(abundanceColumn))
    {
        ret[, (abundanceColumnRel) := .N / spcount, by = "cluster"]
        ret[, (abundanceColumnAbs) := .N, by = "cluster"]
        ret <- ret[get(abundanceColumnRel) >= minAbundanceRel & get(abundanceColumnAbs) >= minAbundanceAbs]
    }
    
    ret <- doAvgPL(ret, "cluster", spcount)
    ret[, c("cluster", "spid") := NULL]
    setcolorder(ret, intersect(c("mz", "intensity", abundanceColumnRel, abundanceColumnAbs, avgCols), names(ret)))
    
    if (nrow(ret) == 0)
        return(ret)
    
    if (assignPrecursor)
    {
        precMZs <- spcomb[precursor == TRUE, mz]
        if (length(precMZs) > 0)
            ret <- assignPrecursorToMSPeakList(ret, mean(precMZs))
        else
            ret[, precursor := FALSE]
    }
    else
        ret[, precursor := FALSE]
    
    # post intensity filter
    ret <- ret[intensity >= minIntensityPost]
    
    if (pruneMissingPrecursor && !any(ret$precursor))
        return(emptyMSPeakList(abundanceColumn, avgCols))
    
    return(ret)
}

#' @include main.R
NULL

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

averageSpectraMZR <- function(spectra, hd, clusterMzWindow, topMost, minIntensityPre, minIntensityPost, minAbundanceRel,
                              minAbundanceAbs, avgFun, method, precursor, pruneMissingPrecursor, retainPrecursor)
{
    if (nrow(hd) == 0) # no spectra, return empty spectrum
        return(emptyMSPeakList("abundance", NULL))
    
    sp <- spectra$spectra[hd$seqNum]
    # convert to peaklist format
    sp <- lapply(sp, function(spec) setnames(as.data.table(spec), c("mz", "intensity")))
    sp <- lapply(sp, assignPrecursorToMSPeakList, precursor)
    
    return(averageSpectra(sp, clusterMzWindow, topMost, minIntensityPre, minIntensityPost, minAbundanceRel,
                          minAbundanceAbs, "abundance", avgFun, NULL, method, TRUE, pruneMissingPrecursor,
                          retainPrecursor))
}

# use mzR to generate MS peaklists.
# limitations compared to DA: no bg subtraction, no isotope information

#' Generate peak lists with mzR (deprecated)
#'
#' Uses the \pkg{mzR} package to read the MS data needed for MS peak lists. This function is now deprecated, please use
#' \code{\link{generateMSPeakLists}} instead.
#'
#' @details The MS data files should be either in \file{.mzXML} or \file{.mzML} format.
#'
#' @template centroid_note_mandatory
#'
#' @param precursorMzWindow The \emph{m/z} window (in Da) to find MS/MS spectra of a precursor. This is typically used
#'   for Data-Dependent like MS/MS data and should correspond to the isolation \emph{m/z} window (\emph{i.e.} +/- the
#'   precursor \emph{m/z}) that was used to collect the data. For Data-Independent MS/MS experiments, where precursor
#'   ions are not isolated prior to fragmentation (\emph{e.g.} bbCID, MSe, all-ion, ...) the value should be
#'   \code{NULL}.
#' @param topMost Only extract MS peak lists from a maximum of \code{topMost} analyses with highest intensity. If
#'   \code{NULL} all analyses will be used.
#' @param avgFeatParams Parameters used for averaging MS peak lists of individual features. Analogous to
#'   \code{avgFGroupParams}.
#'
#' @template mspl_algo-args
#' @inheritParams generateMSPeakLists
#'
#' @inherit generateMSPeakLists return
#'
#' @references \addCitations{mzR}
#'
#' @templateVar what generateMSPeakListsMzR
#' @template main-rd-method
#' @keywords internal
#' @export
setMethod("generateMSPeakListsMzR", "featureGroups", function(fGroups, maxMSRtWindow = 5,
                                                              precursorMzWindow = 4, topMost = NULL,
                                                              avgFeatParams = getDefAvgPListParams(clusterMzWindow = 0.005),
                                                              avgFGroupParams = getDefAvgPListParams(clusterMzWindow = 0.005,
                                                                                                     withPrecursorMS = FALSE))
{
    .Deprecated(old = "generateMSPeakListsMzR", new = "generateMSPeakLists")
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertNumber(maxMSRtWindow, lower = 1, finite = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertNumber(precursorMzWindow, lower = 0, finite = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertCount(topMost, positive = TRUE, null.ok = TRUE, add = ac)
    assertAvgPListParams(avgFeatParams, add = ac)
    assertAvgPListParams(avgFGroupParams, add = ac)
    checkmate::reportAssertions(ac)

    ftindex <- groupFeatIndex(fGroups)
    fTable <- featureTable(fGroups)
    anaInfo <- analysisInfo(fGroups)
    filePaths <- getCentroidedMSFilesFromAnaInfo(anaInfo)
    gTable <- groupTable(fGroups)
    gCount <- length(fGroups)
    gNames <- names(fGroups)
    anaCount <- nrow(anaInfo)

    if (gCount == 0)
        return(MSPeakLists(algorithm = "mzr"))

    cacheDB <- openCacheDBScope()
    setHash <- makeHash(fGroups, maxMSRtWindow, precursorMzWindow, topMost, avgFeatParams)
    cachedSet <- loadCacheSet("MSPeakListsMzR", setHash, cacheDB)
    resultHashes <- vector("character", anaCount * gCount)
    resultHashCount <- 0

    warning("The withPrecursorMS, minIntensityIMS, minAbundanceIMSRel, and minAbundanceIMSAbs parameters are not ",
            "supported with mzR, they will be ignored. ",
            "Furthermore, the avgFun parameter was removed, and now default to 'mean'.", call. = FALSE)
    avgFeatParams <- avgFeatParams[setdiff(names(avgFeatParams),
                                           c("withPrecursorMS", "minIntensityIMS", "minAbundanceIMSRel",
                                             "minAbundanceIMSAbs"))]
    avgFeatParams$avgFun <- mean
    
    avgFeatParamsMS <- avgFeatParamsMSMS <-
        avgFeatParams[setdiff(names(avgFeatParams), c("pruneMissingPrecursorMS", "retainPrecursorMSMS"))]
    avgFeatParamsMS$retainPrecursor <- TRUE;
    avgFeatParamsMS$pruneMissingPrecursor <- avgFeatParams$pruneMissingPrecursorMS
    avgFeatParamsMSMS$pruneMissingPrecursor <- FALSE
    avgFeatParamsMSMS$retainPrecursor <- avgFeatParams$retainPrecursorMSMS

    # structure: [[analysis]][[fGroup]][[MSType]][[MSPeak]]
    plists <- list()
    metadata <- list()

    # if topMost is specified, make list (topAna) of topMost intense analyses
    if (!is.null(topMost) && topMost > 0)
        topAna <- sapply(seq_along(gTable), function(grpi) anaInfo$analysis[order(gTable[[grpi]], decreasing = TRUE)[seq_len(topMost)]])

    for (anai in seq_len(anaCount))
    {
        ana <- anaInfo$analysis[anai]
        spectra <- NULL

        baseHash <- makeHash(ana, maxMSRtWindow, precursorMzWindow, topMost, avgFeatParams)

        printf("Loading all MS peak lists for %d feature groups in analysis '%s'...\n", gCount, ana)
        prog <- openProgBar(0, gCount)

        for (grpi in seq_along(ftindex))
        {
            if (!is.null(topMost) && !ana %in% topAna[[grpi]])
                next # not intense enough

            grp <- gNames[grpi]

            fti <- ftindex[[grpi]][anai]
            if (fti == 0)
                next
            ft <- fTable[[ana]][fti]

            hash <- makeHash(baseHash, ft)
            resultHashCount <- resultHashCount + 1
            resultHashes[resultHashCount] <- hash

            results <- NULL
            if (!is.null(cachedSet))
                results <- cachedSet[[hash]]
            if (is.null(results))
                results <- loadCacheData("MSPeakListsMzR", hash, cacheDB)

            if (is.null(results))
            {
                results <- list(plists = list(), metatadata = list())

                rtRange <- c(ft$retmin, ft$retmax)
                if (!is.null(maxMSRtWindow) && diff(rtRange) > maxMSRtWindow*2)
                    rtRange <- c(max(rtRange[1], ft$ret - maxMSRtWindow), min(rtRange[2], ft$ret + maxMSRtWindow))

                if (is.null(spectra))
                    spectra <- loadSpectra(filePaths[anai], verbose = FALSE)

                # NOTE: precursor is set here only for precursor assignment,
                # keeping precursorMzWindow unset for the header makes sure that
                # no spectra selection is made.
                results$metadata$MS <- getSpectraHeader(spectra, rtRange, 1, NULL, NULL)
                results$plists$MS <- do.call(averageSpectraMZR,
                                             c(list(spectra = spectra, hd = results$metadata$MS, precursor = ft$mz), avgFeatParamsMS))

                hdMSMS <- getSpectraHeader(spectra, rtRange, 2, ft$mz, precursorMzWindow)
                MSMS <- do.call(averageSpectraMZR, c(list(spectra = spectra, hd = hdMSMS, precursor = ft$mz),
                                                     avgFeatParamsMSMS))
                if (nrow(MSMS) > 0)
                {
                    results$plists$MSMS <- MSMS
                    results$metadata$MSMS <- hdMSMS
                }

                saveCacheData("MSPeakListsMzR", results, hash, cacheDB)
            }

            plists[[ana]][[grp]] <- results$plists
            metadata[[ana]][[grp]] <- results$metadata

            setTxtProgressBar(prog, grpi)
        }

        setTxtProgressBar(prog, gCount)
        close(prog)
    }

    if (is.null(cachedSet))
        saveCacheSet("MSPeakListsMzR", resultHashes[seq_len(resultHashCount)], setHash, cacheDB)

    return(MSPeakLists(peakLists = plists, metadata = metadata, avgPeakListArgs = avgFGroupParams,
                       origFGNames = gNames, algorithm = "mzr"))
})

#' @rdname generateMSPeakListsMzR
#' @export
setMethod("generateMSPeakListsMzR", "featureGroupsSet", function(fGroups, ...)
{
    generateMSPeakListsSet(fGroups, generateMSPeakListsMzR, ...)
})
