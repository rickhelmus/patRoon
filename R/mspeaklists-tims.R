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

    # UNDONE: check validity input data
    
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
    setHash <- makeHash(fGroups, maxMSRtWindow, topMost, avgFeatParams)
    cachedSet <- loadCacheSet("MSPeakListsTIMS", setHash, cacheDB)
    resultHashes <- vector("character", anaCount * gCount)
    resultHashCount <- 0

    avgFeatParamsMS <- avgFeatParamsMSMS <-
        avgFeatParams[setdiff(names(avgFeatParams), c("pruneMissingPrecursorMS", "retainPrecursorMSMS"))]
    avgFeatParamsMS$retainPrecursor <- TRUE;
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
        fp <- getBrukerAnalysisPath(ana, anaInfo$path[anai]) # UNDONE: streamline this function mz(X)ML versions (also for DA algos)
        
        TIMSDB <- openTIMSMetaDBScope(f = fp)
        frames <- getTIMSMetaTable(TIMSDB, "Frames", c("Id", "Time", "MsMsType"))
        framesMS <- frames[MsMsType == 0]
        pasefInfo <- getTIMSMetaTable(TIMSDB, "PasefFrameMsMsInfo", c("Frame", "ScanNumBegin", "ScanNumEnd",
                                                                      "IsolationMz", "isolationWidth"))

        # UNDONE: file hash? (also for DA)
        baseHash <- makeHash(ana, maxMSRtWindow, topMost, avgFeatParams)
        
        fTableAna <- copy(fTable[[ana]])
        if (!is.null(maxMSRtWindow))
        {
            fTableAna[, retmin := max(retmin, ret - maxMSRtWindow), by = seq_len(nrow(fTableAna))]
            fTableAna[, retmax := min(retmax, ret + maxMSRtWindow), by = seq_len(nrow(fTableAna))]
        }
        fTableAna[, frameIDsMS := list(list(framesMS[Time %between% c(retmin, retmax)]$Id)), by = seq_len(nrow(fTableAna))]
        
        # UNDONE: more args, currently missing in PListParams
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
        msplMS <- lapply(msplMS, as.data.table)
        msplMS <- Map(msplMS, fTableAna$mz, f = assignPrecursorToMSPeakList)
        mspl <- sapply(msplMS, function(x) list(MS = x), simplify = FALSE)  # move to MS sub fields
        
        if (nrow(pasefInfo) > 0)
        {
            fTableAna[, c("frameIDsMSMS", "scanStarts", "scanEnds") := {
                # UNDONE: window needs to be halved?
                pinf <- pasefInfo[numGTE(mz, IsolationMz - (IsolationWidth/2)) &
                                      numLTE(mz, IsolationMz + (IsolationWidth/2))]
                pinf <- pinf[Frame %in% frames[Time %between% c(retmin, retmax)]$Id]
                
                # NOTE: lots of nested lists to store list data in DTs...
                list(
                    list(pinf$Frame), # relevant MS/MS frame IDs
                    list(pinf[, list(list(ScanNumBegin)), by = Frame][[2]]), # scanStarts/frame
                    list(pinf[, list(list(ScanNumEnd)), by = Frame][[2]]) # scanEnds/frame
                )
            }, by = seq_len(nrow(fTableAna))]
            fTableAnaMSMS <- fTableAna[lengths(frameIDsMSMS) > 0]

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
            msplMSMS <- lapply(msplMSMS, as.data.table)
            msplMSMS <- Map(msplMSMS, fTableAnaMSMS$mz, f = assignPrecursorToMSPeakList)

            # merge with MS data: combine data for overlapping fGroups
            mspl <- sapply(union(names(mspl), names(msplMSMS)), function(fg) list(MS = mspl[[fg]]$MS,
                                                                                  MSMS = msplMSMS[[fg]]),
                           simplify = FALSE)
            # add missing data
            unFGMSMS <- setdiff(names(msplMSMS), names(mspl))
            mspl[unFGMSMS] <- lapply(msplMSMS[unFGMSMS], function(x) list(MSMS = x))
        }
        
        mspl <- mspl[intersect(gNames, names(mspl))] # sync fg order
        
        setTxtProgressBar(prog, anai)
        return(mspl)
    })
    names(plists) <- anaInfo$analysis
    
    close(prog)

    # if (is.null(cachedSet))
    #     saveCacheSet("MSPeakListsTIMS", resultHashes[seq_len(resultHashCount)], setHash, cacheDB)

    return(MSPeakLists(peakLists = plists, metadata = list(), avgPeakListArgs = avgFGroupParams,
                       origFGNames = gNames, algorithm = "tims"))
})

#' @rdname generateMSPeakListsTIMS
#' @export
setMethod("generateMSPeakListsTIMS", "featureGroupsSet", function(fGroups, ...)
{
    generateMSPeakListsSet(fGroups, generateMSPeakListsTIMS, ...)
})
