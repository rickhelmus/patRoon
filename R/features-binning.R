#' @include features.R
NULL

#' @rdname features-class
#' @export
featuresBinning <- setClass("featuresBinning", contains = "features")

setMethod("initialize", "featuresBinning",
          function(.Object, ...) callNextMethod(.Object, algorithm = "binning", ...))

#' @export
findFeaturesBinning <- function(analysisInfo, retRange = NULL, mzRange = c(50, 400), mzStep = 0.02, mobRange = NULL,
                                mobStep = 0.04, minIntensityIMS = 25, peaksParam, verbose = TRUE)
{
    # UNDONE: add refs to docs, and highlight changes
    # UNDONE: mobRange/mobStep defaults
    # UNDONE: default OK for minIntensityIMS?
    # UNDONE: print messages, stats etc, taking verbose into account
    # UNDONE: support parallel? Can take a lot of memory and currently not supported by applyMSData()
    # UNDONE: use BP intensity?
    
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, add = ac)
    aapply(assertRange, . ~ retRange + mzRange + mobRange, null.ok = c(TRUE, FALSE, TRUE), fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ mzStep + mobStep, lower = 0.000001, finite = TRUE, fixed = list(add = ac))
    assertFindPeaksParam(peaksParam, add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    if (is.null(retRange))
        retRange <- c(0, 0)
    
    cacheDB <- openCacheDBScope()
    baseHash <- makeHash(retRange, mzRange, mzStep, mobRange, mobStep, minIntensityIMS, peaksParam)
    anaHashes <- getMSFileHashesFromAvailBackend(analysisInfo, needIMS = withIMS)
    anaHashes <- sapply(anaHashes, makeHash, baseHash)
    cachedData <- pruneList(loadCacheData("featuresBinning", anaHashes, simplify = FALSE, dbArg = cacheDB))
    if (length(cachedData) > 0)
    {
        names(cachedData) <- names(anaHashes)[match(names(cachedData), anaHashes)]
        anaInfoTBD <- analysisInfo[!analysis %in% names(cachedData)]
    }
    else
        anaInfoTBD <- analysisInfo
    
    getEICInfoList <- function(withIMS)
    {
        # UNDONE: also support other binning approaches?
        binsMZ <- seq(mzRange[1], mzRange[2], by = mzStep * 0.5)
        names(binsMZ) <- paste0("bin_M", binsMZ)
        
        binsIMS <- NULL
        if (withIMS)
        {
            binsIMS <- seq(mobRange[1], mobRange[2], by = mobStep * 0.5)
            names(binsIMS) <- paste0("bin_I", binsIMS)
        }
        
        EICInfoAna <- data.table(mzmin = binsMZ, mzmax = binsMZ + mzStep, retmin = retRange[1], retmax = retRange[2],
                                 EIC_ID_MZ = names(binsMZ))
        if (withIMS)
        {
            tab <- CJ(EIC_ID_MZ = names(binsMZ), EIC_ID_IMS = names(binsIMS), sorted = FALSE)
            tab[, c("mobmin", "mobmax") := .(binsIMS[EIC_ID_IMS], binsIMS[EIC_ID_IMS] + mobStep)]
            EICInfoAna <- merge(EICInfoAna, tab, by = "EIC_ID_MZ", sort = FALSE)
        }
        else
            EICInfoAna[, c("mobmin", "mobmax") := 0]
        
        EICInfoAna[, EIC_ID := paste0("EIC_", .I)]
        
        return(setNames(rep(list(EICInfoAna), nrow(anaInfoTBD)), anaInfoTBD$analysis))
    }
    
    getEICsAna <- function(backend, EICInfo)
    {
        # UNDONE: make post filters configurable
        ret <- doGetEICsForAna(backend, EICInfo$mzmin, EICInfo$mzmax, EICInfo$retmin, EICInfo$retmax,
                               EICInfo$mobmin, EICInfo$mobmax, mzExpIMSWindow = 0, minIntensityIMS = minIntensityIMS,
                               compress = FALSE, showProgress = FALSE, withBP = TRUE, minEICIntensity = 1000,
                               minAdjacentTime = 30, minAdjacentPointIntensity = 250)
        names(ret) <- EICInfo$EIC_ID
        ret <- pruneList(ret, checkEmptyElements = TRUE)
        return(ret)
    }
    
    fList <- list()
    if (nrow(anaInfoTBD) > 0)
    {
        withIMS <- !is.null(mobRange)
        EICInfoListMZ <- getEICInfoList(withIMS = FALSE)
        EICInfoListMob <- if (withIMS) getEICInfoList(withIMS = TRUE) else rep(list(NULL), nrow(anaInfoTBD))
        
        fList <- applyMSData(anaInfoTBD, EICInfoListMZ, EICInfoListMob, needIMS = withIMS,
                             func = function(ana, path, backend, EICInfoMZ, EICInfoMob)
        {
            openMSReadBackend(backend, path)
            
            EICs <- getEICsAna(backend, EICInfoMZ)
            # omit missing
            EICInfoMZ <- EICInfoMZ[EIC_ID %chin% names(EICs)]
            
            if (withIMS)
            {
                # With IMS worksflows the mz EICs are only used as a pre-filter. As the EIC object is potentially large we
                # remove it here.
                rm(EICs) 
                gc()
                
                # remove complete m/z bins that were filtered out before
                temp <- EICInfoMZ[, c("mzmin", "mzmax"), with = FALSE]
                setkeyv(temp, c("mzmin", "mzmax"))
                ov <- foverlaps(EICInfoMob, temp, type = "within", nomatch = NULL, which = TRUE)
                EICInfo <- EICInfoMob[ov$xid]
                
                EICs <- getEICsAna(backend, EICInfo)
            }
            else
                EICInfo <- EICInfoMZ
            
            peaks <- findPeaksInEICs(EICs, peaksParam, withBP = TRUE, withMobility = withIMS, cacheDB = cacheDB)
            
            # only keep those peaks with m/z in the "center" of the analyzed m/z range
            peaks[, binMZStart := EICInfo[match(peaks$EIC_ID, EIC_ID)]$mzmin]
            peaks <- peaks[between(mz, binMZStart + mzStep/4, binMZStart + mzStep/4*3) == TRUE]
            if (withIMS)
            {
                peaks[, binMobStart := EICInfo[match(peaks$EIC_ID, EIC_ID)]$mobmin]
                peaks <- peaks[between(mobility, binMobStart + mobStep/4, binMobStart + mobStep/4*3) == TRUE]
                
                if (nrow(peaks) > 1)
                {
                    # remove duplicates and keep the one with the highest intensity
                    peaks[, keep := TRUE]
                    for (r in seq_len(nrow(peaks) - 1))
                    {
                        if (!peaks$keep[r])
                            next
                        rTab <- peaks[r]
                        nextTab <- peaks[seq(r + 1, nrow(peaks))]
                        # UNDONE: make thresholds configurable
                        nextTab <- nextTab[keep == TRUE & numLTE(abs(ret - rTab$ret), 3) & numLTE(abs(mz - rTab$mz), 0.001) &
                                               numLTE(abs(mobility - rTab$mobility), 0.02)]
                        if (nrow(nextTab) > 0)
                        {
                            maxInt <- max(rTab$intensity, nextTab$intensity)
                            peaks[ID %chin% c(rTab$ID, nextTab$ID), keep := numEQ(intensity, maxInt)]
                        }
                    }
                    peaks <- peaks[keep == TRUE][, keep := NULL]
                }
            }
            peaks <- removeDTColumnsIfPresent(peaks, c("binMZStart", "binMobStart"))
            
            doProgress()
            
            return(peaks)
        })
        
        for (a in anaInfoTBD$analysis)
            saveCacheData("featuresBinning", fList[[a]], anaHashes[[a]], dbArg = cacheDB)
    }
    
    if (length(cachedData) > 0)
    {
        fList <- c(fList, cachedData)
        fList <- fList[analysisInfo$analysis] # put original order
    }
    
    if (verbose)
    {
        printf("Done!\n")
        printFeatStats(fList)
    }
    
    return(featuresBinning(analysisInfo = analysisInfo, features = fList, hasMobilities = withIMS))
}
