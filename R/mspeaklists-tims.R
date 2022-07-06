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
                method = "dist",
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

    # structure: [[analysis]][[fGroup]][[MSType]][[MSPeak]]
    plists <- list()

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
    
    for (anai in seq_len(anaCount))
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
        fTableAna[, retmin := max(retmin, ret - maxMSRtWindow), by = seq_len(nrow(fTableAna))]
        fTableAna[, retmax := min(retmax, ret + maxMSRtWindow), by = seq_len(nrow(fTableAna))]
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
        names(mspl) <- fTableAna$group
        
        if (nrow(pasefInfo) > 0)
        {
            fTableAna[, frameIDsMSMS := {
                # UNDONE: window needs to be halved?
                framesPASEF <- pasefInfo[mz %between% ((IsolationMz + c(-IsolationWidth, IsolationWidth))/2)]$Frame
                intersect(framesPASEF, framesMS[[1]])
            }, by = seq_len(nrow(fTableAna))]
            fTableAnaMSMS <- fTableAna[lengths(frameIDsMSMS) > 0]
        }
        
        setTxtProgressBar(prog, anai)
    }
    
    close(prog)

    if (is.null(cachedSet))
        saveCacheSet("MSPeakListsTIMS", resultHashes[seq_len(resultHashCount)], setHash, cacheDB)

    return(MSPeakLists(peakLists = plists, metadata = list(), avgPeakListArgs = avgFGroupParams,
                       origFGNames = gNames, algorithm = "tims"))
})

#' @rdname generateMSPeakListsTIMS
#' @export
setMethod("generateMSPeakListsTIMS", "featureGroupsSet", function(fGroups, ...)
{
    generateMSPeakListsSet(fGroups, generateMSPeakListsTIMS, ...)
})
