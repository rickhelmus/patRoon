#' @include features.R
NULL

#' @rdname features-class
#' @export
featuresBinning <- setClass("featuresBinning", contains = "features")

setMethod("initialize", "featuresBinning",
          function(.Object, ...) callNextMethod(.Object, algorithm = "binning", ...))

#' @export
findFeaturesBinning <- function(analysisInfo, retRange = NULL, mzRange = c(50, 400), mzStep = 0.02, mobRange = NULL,
                                mobStep = 0.04, minIntensityIMS = 25,
                                peaksParam, ..., parallel = TRUE, verbose = TRUE)
{
    # UNDONE: add refs to docs, and highlight changes
    # UNDONE: mobRange/mobStep defaults
    # UNDONE: default OK for minIntensityIMS?
    # UNDONE: print messages, stats etc, taking verbose into account
    
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, add = ac)
    aapply(assertRange, . ~ retRange + mzRange + mobRange, null.ok = c(TRUE, FALSE, TRUE), fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ mzStep + mobStep, lower = 0.000001, finite = TRUE, fixed = list(add = ac))
    assertFindPeaksParam(peaksParam, add = ac)
    aapply(checkmate::assertFlag, . ~ parallel + verbose, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    if (is.null(retRange))
        retRange <- c(0, 0)
    
    cacheDB <- openCacheDBScope()
    baseHash <- makeHash(mzRange, mzStep, mobRange, mobStep, minIntensityIMS, peaksParam)
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
        
        EICInfoAna[, EIC_ID := paste0("EIC_", .I)]
        
        return(setNames(rep(list(EICInfoAna), nrow(anaInfoTBD)), anaInfoTBD$analysis))
    }
    
    getAllEICs <- function(EICInfoList, ...)
    {
        ret <- doGetEICs(anaInfoTBD, EICInfoList, minIntensityIMS = minIntensityIMS, compress = FALSE,
                         showProgress = if (verbose) "ana" else FALSE, withBP = TRUE, cacheDB = cacheDB, ...)
        ret <- lapply(ret, setNames, EICInfoList[[1]]$EIC_ID)
        return(ret)
    }
    
    fList <- list()
    if (nrow(anaInfoTBD) > 0)
    {
        EICInfoList <- getEICInfoList(withIMS = FALSE)
        allEICs <- getAllEICs(EICInfoList, minEICIntensity = 1000, minAdjacentTime = 30,
                              minAdjacentPointIntensity = 250)
        allEICs <- lapply(allEICs, pruneList, checkEmptyElements = TRUE)
        # omit missing
        EICInfoList <- Map(EICInfoList, allEICs, f = function(info, eics) info[EIC_ID %chin% names(eics)])

        withIMS <- !is.null(mobRange)
        if (withIMS)
        {
            # With IMS worksflows the mz EICs are only used as a pre-filter. As the EIC object is potentially large we
            # remove it here.
            rm(allEICs) 
            gc()
            
            mzEICInfoList <- EICInfoList
            EICInfoList <- getEICInfoList(withIMS = TRUE)
            EICInfoList <- Map(EICInfoList, mzEICInfoList, f = function(infoMob, infoMZ)
            {
                # remove complete m/z bins that were filtered out before
                infoMZ <- infoMZ[, c("mzmin", "mzmax"), with = FALSE]
                setkeyv(infoMZ, c("mzmin", "mzmax"))
                ov <- foverlaps(infoMob, infoMZ, type = "within", nomatch = NULL, which = TRUE)
                return(infoMob[ov$xid])
            })
            
            allEICs <- getAllEICs(EICInfoList, minEICIntensity = 1000, minAdjacentTime = 30,
                                  minAdjacentPointIntensity = 250)
            allEICs <- lapply(allEICs, pruneList, checkEmptyElements = TRUE)
        }
        fList <- findPeaksInEICs(allEICs, peaksParam, withBP = TRUE, withMobility = withIMS, parallel = parallel,
                                 cacheDB = cacheDB)
        fList <- Map(fList, EICInfoList, f = function(fTab, EICInfoAna)
        {
            # only keep those peaks with m/z in the "center" of the analyzed m/z range
            fTab[, binMZStart := EICInfoAna[match(fTab$EIC_ID, EIC_ID)]$mzmin]
            fTab <- fTab[between(mz, binMZStart + mzStep/4, binMZStart + mzStep/4*3) == TRUE]
            if (withIMS)
            {
                fTab[, binMobStart := EICInfoAna[match(fTab$EIC_ID, EIC_ID)]$mobmin]
                fTab <- fTab[between(mobility, binMobStart + mobStep/4, binMobStart + mobStep/4*3) == TRUE]
                
                if (nrow(fTab) > 1)
                {
                    # remove duplicates and keep the one with the highest intensity
                    fTab[, keep := TRUE]
                    for (r in seq_len(nrow(fTab) - 1))
                    {
                        if (!fTab$keep[r])
                            next
                        rTab <- fTab[r]
                        nextTab <- fTab[seq(r + 1, nrow(fTab))]
                        # UNDONE: make thresholds configurable
                        nextTab <- nextTab[keep == TRUE & numLTE(abs(ret - rTab$ret), 3) & numLTE(abs(mz - rTab$mz), 0.001) &
                                               numLTE(abs(mobility - rTab$mobility), 0.02)]
                        if (nrow(nextTab) > 0)
                        {
                            maxInt <- max(rTab$intensity, nextTab$intensity)
                            fTab[ID %chin% c(rTab$ID, nextTab$ID), keep := numEQ(intensity, maxInt)]
                        }
                    }
                    fTab <- fTab[keep == TRUE][, keep := NULL]
                }
            }
            fTab <- removeDTColumnsIfPresent(fTab, c("binMZStart", "binMobStart"))
            return(fTab)
        })
        
        # UNDONE: use BP intensity?
        
        for (a in anaInfoTBD$analysis)
            saveCacheData("featuresBinning", fList[[a]], anaHashes[[a]])
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
