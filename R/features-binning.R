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
    
    getEICInfoList <- function(withIMS, wide)
    {
        # UNDONE: make factor configurable
        mzst <- if (wide) mzStep * 1 else mzStep
        mobst <- if (wide) mobStep * 1 else mobStep
        
        # UNDONE: also support other binning approaches?
        binsMZ <- seq(mzRange[1], mzRange[2], by = mzst * 0.5)
        names(binsMZ) <- paste0("bin_M", binsMZ)
        
        binsMob <- NULL
        if (withIMS)
        {
            binsMob <- seq(mobRange[1], mobRange[2], by = mobst * 0.5)
            names(binsMob) <- paste0("bin_M", binsMob)
        }
        
        EICInfoAna <- data.table(mzmin = binsMZ, mzmax = binsMZ + mzst, retmin = retRange[1], retmax = retRange[2],
                                 EIC_ID_MZ = names(binsMZ))
        if (withIMS)
        {
            tab <- CJ(EIC_ID_MZ = names(binsMZ), EIC_ID_mob = names(binsMob), sorted = FALSE)
            tab[, c("mobmin", "mobmax") := .(binsMob[EIC_ID_mob], binsMob[EIC_ID_mob] + mobst)]
            EICInfoAna <- merge(EICInfoAna, tab, by = "EIC_ID_MZ", sort = FALSE)
        }
        
        EICInfoAna[, EIC_ID := paste0("EIC_", .I)] # UNDONE
        
        return(setNames(rep(list(EICInfoAna), nrow(anaInfoTBD)), anaInfoTBD$analysis))
    }
    
    getAllEICs <- function(EICInfoList, ...)
    {
        allEICs <- doGetEICs(anaInfoTBD, EICInfoList, minIntensityIMS = minIntensityIMS, compress = FALSE,
                             showProgress = if (verbose) "ana" else FALSE, withBP = TRUE, cacheDB = cacheDB, ...)
        allEICs <- lapply(allEICs, setNames, EICInfoList[[1]]$EIC_ID)
        return(allEICs)
    }
    
    fList <- list()
    if (nrow(anaInfoTBD) > 0)
    {
        wideEICInfoList <- getEICInfoList(withIMS = FALSE, wide = TRUE)
        wideEICs <- getAllEICs(wideEICInfoList, minEICIntensity = 1000, minAdjacentTime = 30,
                               minAdjacentPointIntensity = 250)
        validWideEICs <- lapply(wideEICs, function(wea) lengths(wea) > 0)
        
        EICInfoList <- getEICInfoList(withIMS = !is.null(mobRange), wide = FALSE)
        EICInfoList <- Map(EICInfoList, wideEICInfoList, validWideEICs, f = function(EICInfoAna, EICInfoAnaWide, validWideEICs)
        {
            wide <- EICInfoAnaWide[validWideEICs == TRUE, c("EIC_ID", "mzmin", "mzmax"), with = FALSE]
            setkeyv(wide, c("mzmin", "mzmax"))
            ov <- foverlaps(EICInfoAna, wide, type = "within", nomatch = NULL, which = TRUE)
            return(EICInfoAna[ov$xid])
        })
        
        allEICs <- getAllEICs(EICInfoList, minEICIntensity = 1000, minAdjacentTime = 30,
                              minAdjacentPointIntensity = 250)
        allEICs <- lapply(allEICs, pruneList, checkEmptyElements = TRUE)

        withIMS <- !is.null(mobRange)
        fList <- findPeaksInEICs(allEICs, peaksParam, withBP = TRUE, withMobility = withIMS, parallel = parallel,
                                 cacheDB = cacheDB)
        fList <- Map(fList, EICInfoList, f = function(fTab, EICInfoAna)
        {
            # only keep those peaks with m/z in the "center" of the analyzed m/z range
            fTab[, EIC_ID_MZ := EICInfoAna[match(fTab$EIC_ID, EIC_ID)]$EIC_ID_MZ]
            fTab <- fTab[between(mz, binsMZ[EIC_ID_MZ] + mzStep/4, binsMZ[EIC_ID_MZ] + mzStep/4*3) == TRUE]
            if (withIMS)
            {
                fTab[, EIC_ID_mob := EICInfoAna[match(fTab$EIC_ID, EIC_ID)]$EIC_ID_mob]
                fTab[, mob_bin := binsMob[EIC_ID_mob]] # UNDONE
                fTab <- fTab[between(mobility, binsMob[EIC_ID_mob] + mobStep/4, binsMob[EIC_ID_mob] + mobStep/4*3) == TRUE]
                
                if (nrow(fTab) > 1)
                {
                    # remove duplicates and keep the one with the highest intensity
                    fTab[, keep := TRUE]
                    for (r in seq_len(nrow(fTab)-1))
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
            fTab <- removeDTColumnsIfPresent(fTab, c("EIC_ID_MZ", "EIC_ID_mob", "EIC_ID"))
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
