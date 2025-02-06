#' @include features.R
NULL

removeDuplicateFeatsSusps <- function(tab, checkRet, selectTopIntens)
{
    # marks duplicates in feature or suspect table, and keeps the feature with the highest intensity

    if (nrow(tab) <= 1)
        return(tab)

    isSame <- function(x, y)
    {
        # UNDONE: make thresholds configurable
        yes <- numLTE(abs(x$mz - y$mz), 0.001)
        if (checkRet)
            yes <- yes & ((is.na(x$ret) & is.na(y$ret)) | numLTE(abs(x$ret - y$ret), 3))
        if (!is.null(x[["mobility"]]))
            yes <- yes & numLTE(abs(x$mobility - y$mobility), 0.02)
        return(yes)
    }
    
    tab <- copy(tab)
    tab[, duplicate := FALSE]
    for (r in seq_len(nrow(tab) - 1))
    {
        if (tab$duplicate[r])
            next
        rTab <- tab[r]
        nextTab <- tab[seq(r + 1, nrow(tab))]
        # UNDONE: make thresholds configurable
        nextTab <- nextTab[duplicate == FALSE & isSame(rTab, nextTab)]
        if (nrow(nextTab) > 0)
        {
            if (selectTopIntens)
            {
                maxInt <- max(rTab$intensity, nextTab$intensity)
                tab[ID %chin% c(rTab$ID, nextTab$ID), duplicate := !numEQ(intensity, maxInt)]
            }
            else
            {
                # otherwise just keep current and mark all from next rows
                # UNDONE: name/ID
                tab[name %chin% nextTab$name, duplicate := TRUE]
            }
        }
    }
    return(tab[duplicate == FALSE][, duplicate := NULL])
}

#' @rdname features-class
#' @export
featuresBinning <- setClass("featuresBinning", contains = "features")

setMethod("initialize", "featuresBinning",
          function(.Object, ...) callNextMethod(.Object, algorithm = "binning", ...))

#' @export
findFeaturesBinning <- function(analysisInfo, peaksParam, retRange = NULL, mzRange = c(50, 400), mzStep = 0.02,
                                mobRange = NULL, mobStep = 0.04, suspects = NULL, suspRTWindow = 30,
                                suspMZWindow = 0.005, suspIMSWindow = 0.02, skipInvalid = TRUE, adduct = NULL,
                                prefCalcChemProps = TRUE, neutralChemProps = FALSE, minIntensityIMS = 25, verbose = TRUE)
{
    # UNDONE: add refs to docs, and highlight changes
    # UNDONE: mobRange/mobStep defaults
    # UNDONE: default OK for minIntensityIMS?
    # UNDONE: print messages, stats etc, taking verbose into account
    # UNDONE: support parallel? Can take a lot of memory and currently not supported by applyMSData()
    # UNDONE: use BP intensity?
    # UNDONE: log/print status messages
    # UNDONE: support IM bins with suspects
    # UNDONE: split function for bins / suspects
    
    checkmate::assertFlag(skipInvalid) # not in assert collection, should fail before assertSuspectList
    
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, add = ac)
    assertFindPeaksParam(peaksParam, add = ac)
    aapply(assertRange, . ~ retRange + mzRange + mobRange, null.ok = c(TRUE, FALSE, TRUE), fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ mzStep + mobStep, lower = 0.000001, finite = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ suspRTWindow + suspMZWindow + suspIMSWindow + minIntensityIMS, lower = 0,
           finite = TRUE, fixed = list(add = ac))
    assertSuspectList(suspects, needsAdduct = is.null(adduct), skipInvalid = skipInvalid, null.ok = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ prefCalcChemProps + neutralChemProps + verbose, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    if (!is.null(adduct))
        adduct <- checkAndToAdduct(adduct)
    
    if (!is.null(suspects))
    {
        suspects <- prepareSuspectList(suspects, adduct = adduct, skipInvalid = skipInvalid, checkDesc = TRUE,
                                       prefCalcChemProps = prefCalcChemProps, neutralChemProps = neutralChemProps)
        suspects <- removeDuplicateFeatsSusps(suspects, FALSE, FALSE)
        suspects <- suspects[order(mz)]
    }
    
    doBinning <- is.null(suspects) # UNDONE
    
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
        if (doBinning)
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
        }
        else
        {
            EICInfoAna <- data.table(mzmin = suspects$mz - suspMZWindow, mzmax = suspects$mz + suspMZWindow,
                                     EIC_ID = paste0("EIC_", seq_len(nrow(suspects))))

            # UNDONE: also put in retmin/retmax? This will affect peak finding (eg noise estimation)...
            EICInfoAna[, c("retmin", "retmax") := .(0, 0)]
            
            if (withIMS && !is.null(suspects[["mobility"]]))
                EICInfoAna[, c("mobmin", "mobmax") := .(suspects$mobility - suspIMSWindow, suspects$mobility + suspIMSWindow)]
            else
                EICInfoAna[, c("mobmin", "mobmax") := .(0, 0)]
        }
        
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

            if (doBinning)
            {
                # only keep those peaks with m/z in the "center" of the analyzed m/z and mobility range
                peaks[, binMZStart := EICInfo[match(peaks$EIC_ID, EIC_ID)]$mzmin]
                peaks <- peaks[between(mz, binMZStart + mzStep/4, binMZStart + mzStep/4*3) == TRUE]
                if (withIMS)
                {
                    peaks[, binMobStart := EICInfo[match(peaks$EIC_ID, EIC_ID)]$mobmin]
                    peaks <- peaks[between(mobility, binMobStart + mobStep/4, binMobStart + mobStep/4*3) == TRUE]
                }
            }

            peaks <- removeDuplicateFeatsSusps(peaks, TRUE, TRUE)
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
