#' @include main.R
#' @include mspeaklists.R
NULL

# UNDONE: merge with regular function
getDefAvgPListParamsTIMS <- function(...)
{
    def <- list(clusterMzWindow = 0.01,
                topMost = 50,
                minIntensityPre = 50,
                minIntensityPost = 500,
                minIntensityFinal = 500,
                minAbundance = 2,
                method = "diff",
                pruneMissingPrecursorMS = TRUE,
                retainPrecursorMSMS = TRUE)
    return(modifyList(def, list(...)))
}

setMethod("generateMSPeakListsTIMS", "featureGroups", function(fGroups, maxMSRtWindow = 5, topMost = NULL,
                                                               avgFeatParams = getDefAvgPListParamsTIMS(),
                                                               avgFGroupParams = getDefAvgPListParams())
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertNumber(maxMSRtWindow, lower = 1, finite = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertCount(topMost, positive = TRUE, null.ok = TRUE, add = ac)
    # UNDONE: update assertAvgPListParams(avgFeatParams, add = ac)
    assertAvgPListParams(avgFGroupParams, add = ac)
    checkmate::reportAssertions(ac)

    ftindex <- groupFeatIndex(fGroups)
    fTable <- featureTable(fGroups)
    anaInfo <- analysisInfo(fGroups)
    gTable <- groupTable(fGroups)
    gCount <- length(fGroups)
    gNames <- names(fGroups)
    anaCount <- nrow(anaInfo)

    if (gCount == 0)
        return(MSPeakLists(algorithm = "tims"))

    cacheDB <- openCacheDBScope()

    avgFeatParamsMS <- avgFeatParamsMSMS <-
        avgFeatParams[setdiff(names(avgFeatParams), c("pruneMissingPrecursorMS", "retainPrecursorMSMS"))]
    avgFeatParamsMS$retainPrecursor <- TRUE
    avgFeatParamsMS$pruneMissingPrecursor <- avgFeatParams$pruneMissingPrecursorMS
    avgFeatParamsMSMS$pruneMissingPrecursor <- FALSE
    avgFeatParamsMSMS$retainPrecursor <- avgFeatParams$retainPrecursorMSMS

    # UNDONE: make filter?
    if (!is.null(topMost) && topMost < anaCount)
    {
        fGroups <- delete(fGroups, j = function(gInds, ...)
        {
            sinds <- order(gInds, decreasing = TRUE)
            return(sinds[-seq_len(topMost)])
        })
    }
    
    printf("Loading all MS peak lists for %d feature groups in %d analyses...\n", gCount, anaCount)
    prog <- openProgBar(0, anaCount)
    
    # structure: [[analysis]][[fGroup]][[MSType]][[MSPeak]]
    
    plists <- lapply(seq_len(anaCount), function(anai)
    {
        ana <- anaInfo$analysis[anai]
        fp <- getBrukerAnalysisPath(ana, anaInfo$path[anai])
        fTableAna <- copy(fTable[[ana]])
        
        TIMSDB <- openTIMSMetaDBScope(f = fp)
        frames <- getTIMSMetaTable(TIMSDB, "Frames", c("Id", "Time", "MsMsType"))
        framesMS <- frames[MsMsType == 0]
        pasefInfo <- getTIMSMetaTable(TIMSDB, "PasefFrameMsMsInfo", c("Frame", "ScanNumBegin", "ScanNumEnd",
                                                                      "IsolationMz", "isolationWidth"))

        baseHash <- makeHash(ana, frames, pasefInfo, maxMSRtWindow, topMost, avgFeatParams)
        setHash <- makeHash(fTableAna, baseHash)
        cachedSet <- loadCacheSet("MSPeakListsTIMS", setHash, cacheDB)
        
        # UNDONE: file hash? (also for other algos) Or is frames/pasefInfo enough in this case?
        baseHash <- makeHash(ana, frames, pasefInfo, maxMSRtWindow, topMost, avgFeatParams)
        
        fTableAna[, hash := makeHash(baseHash, .SD), by = seq_len(nrow(fTableAna)), .SDcols = c("mz", "ret", "retmin",
                                                                                                "retmax")]
        cachedData <- setNames(lapply(fTableAna$hash, function(h)
        {
            if (!is.null(cachedSet) && !is.null(cachedSet[[h]]))
                return(cachedSet[[h]])
            return(loadCacheData("MSPeakListsTIMS", h, cacheDB))
        }), fTableAna$group)
        cachedData <- pruneList(cachedData)
        
        fTableAna <- fTableAna[!group %chin% names(cachedData)]
        
        if (nrow(fTableAna) > 0) # still work to do?
        {
            if (!is.null(maxMSRtWindow))
            {
                fTableAna[, retmin := max(retmin, ret - maxMSRtWindow), by = seq_len(nrow(fTableAna))]
                fTableAna[, retmax := min(retmax, ret + maxMSRtWindow), by = seq_len(nrow(fTableAna))]
            }
            
            
            fTableAna[, frameIDsMS := list(list(framesMS[Time %between% c(retmin, retmax)]$Id)), by = seq_len(nrow(fTableAna))]
            
            msplMS <- getTIMSPeakLists(fp, fTableAna$frameIDsMS, precursorMZs = fTableAna$mz,
                                       minIntensityPre = avgFeatParamsMS$minIntensityPre,
                                       minIntensityPost = avgFeatParamsMS$minIntensityPost,
                                       minIntensityFinal = avgFeatParamsMS$minIntensityFinal,
                                       onlyWithPrecursor = avgFeatParamsMS$pruneMissingPrecursor,
                                       method = avgFeatParamsMS$method,
                                       topMost = avgFeatParamsMS$topMost,
                                       mzWindow = avgFeatParamsMS$clusterMzWindow,
                                       minAbundance = avgFeatParamsMS$minAbundance)
            names(msplMS) <- fTableAna$group
            msplMS <- lapply(msplMS, setDT)
            msplMS <- Map(msplMS, fTableAna$mz, f = assignPrecursorToMSPeakList)
            mspl <- sapply(msplMS, function(x) list(MS = x), simplify = FALSE)  # move to MS sub fields
            
            if (nrow(pasefInfo) > 0)
            {
                pasefInfo[, c("isomz_min", "isomz_max") := .(IsolationMz - (IsolationWidth/2),
                                                             IsolationMz + (IsolationWidth/2))]
                # reduce frame search range a bit to optimize below
                framesMSMS <- frames[MsMsType != 0] # UNDONE: check how different types need to be supported
                framesMSMS <- framesMSMS[Time %between% c(min(fTableAna$retmin), max(fTableAna$retmax))]
                
                # only consider PASEF frames of interest
                pasefInfo <- pasefInfo[Frame %in% framesMSMS$Id & isomz_min < max(fTableAna$mz) & isomz_max > min(fTableAna$mz)]
                pasefInfo[, retFrame := framesMSMS[match(Frame, Id)]$Time]
                
                # split all relevant PASEF info rows for each feature, which makes further processing much faster
                frameSubTabs <- Map(fTableAna$mz, fTableAna$retmin, fTableAna$retmax, f = function(m, rmin, rmax)
                {
                    return(pasefInfo[retFrame %between% c(rmin, rmax) & numGTE(m, isomz_min) & numLTE(m, isomz_max)])
                })
                
                fTableAna[, frameIDsMSMS := lapply(frameSubTabs, function(tab) unique(tab$Frame))]
                
                # omit without MSMS to speedup
                withMSMS <- lengths(fTableAna$frameIDsMSMS) > 0
                frameSubTabs <- frameSubTabs[withMSMS]
                fTableAnaMSMS <- fTableAna[withMSMS]
                
                fTableAnaMSMS[, scanStarts := lapply(frameSubTabs, function(tab) tab[, list(list(ScanNumBegin)),
                                                                                     by = "Frame"][[2]])]
                fTableAnaMSMS[, scanEnds := lapply(frameSubTabs, function(tab) tab[, list(list(ScanNumEnd)),
                                                                                   by = "Frame"][[2]])]
                
                msplMSMS <- getTIMSPeakLists(fp, fTableAnaMSMS$frameIDsMSMS, precursorMZs = fTableAnaMSMS$mz,
                                             minIntensityPre = avgFeatParamsMSMS$minIntensityPre,
                                             minIntensityPost = avgFeatParamsMSMS$minIntensityPost,
                                             minIntensityFinal = avgFeatParamsMSMS$minIntensityFinal,
                                             onlyWithPrecursor = avgFeatParamsMSMS$pruneMissingPrecursor,
                                             method = avgFeatParamsMSMS$method,
                                             topMost = avgFeatParamsMSMS$topMost,
                                             mzWindow = avgFeatParamsMSMS$clusterMzWindow,
                                             minAbundance = avgFeatParamsMSMS$minAbundance,
                                             scanStartsListN = fTableAnaMSMS$scanStarts,
                                             scanEndsListN = fTableAnaMSMS$scanEnds)
                names(msplMSMS) <- fTableAnaMSMS$group
                msplMSMS <- lapply(msplMSMS, setDT)
                msplMSMS <- Map(msplMSMS, fTableAnaMSMS$mz, f = assignPrecursorToMSPeakList)
                
                # merge with MS data: combine data for overlapping fGroups
                mspl <- sapply(union(names(mspl), names(msplMSMS)), function(fg) list(MS = mspl[[fg]]$MS,
                                                                                      MSMS = msplMSMS[[fg]]),
                               simplify = FALSE)
                # add missing data
                unFGMSMS <- setdiff(names(msplMSMS), names(mspl))
                mspl[unFGMSMS] <- lapply(msplMSMS[unFGMSMS], function(x) list(MSMS = x))
            }
            
            for (fg in names(mspl))
                saveCacheData("MSPeakListsTIMS", mspl[[fg]], fTableAna[group == fg]$hash, cacheDB)
        }
        
        if (length(cachedData) > 0)
            mspl <- c(mspl, cachedData)
        
        mspl <- mspl[intersect(gNames, names(mspl))] # sync fg order
        
        if (is.null(cachedSet))
            saveCacheSet("MSPeakListsTIMS", fTableAna$hash, setHash, cacheDB)
        
        setTxtProgressBar(prog, anai)
        return(mspl)
    })
    names(plists) <- anaInfo$analysis
    
    close(prog)

    return(MSPeakLists(peakLists = plists, metadata = list(), avgPeakListArgs = avgFGroupParams,
                       origFGNames = gNames, algorithm = "tims"))
})

#' @rdname generateMSPeakListsTIMS
#' @export
setMethod("generateMSPeakListsTIMS", "featureGroupsSet", function(fGroups, ...)
{
    generateMSPeakListsSet(fGroups, generateMSPeakListsTIMS, ...)
})
