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

getFeatEICsInfo <- function(params, withIMS)
{
    EICInfoAna <- NULL
    if (params$methodMZ == "bins")
    {
        # UNDONE: also support other binning approaches?
        binsMZ <- seq(params$mzRange[1], params$mzRange[2], by = params$mzStep * 0.5)
        names(binsMZ) <- paste0("bin_M", binsMZ)
        
        binsIMS <- NULL
        if (!is.null(params[["methodIMS"]]))
        {
            binsIMS <- seq(params$mobRange[1], params$mobRange[2], by = params$mobStep * 0.5)
            names(binsIMS) <- paste0("bin_I", binsIMS)
        }
        
        EICInfoAna <- data.table(mzmin = binsMZ, mzmax = binsMZ + params$mzStep, retmin = params$retRange[1],
                                 retmax = params$retRange[2], EIC_ID_MZ = names(binsMZ))
        if (withIMS)
        {
            tab <- CJ(EIC_ID_MZ = names(binsMZ), EIC_ID_IMS = names(binsIMS), sorted = FALSE)
            tab[, c("mobmin", "mobmax") := .(binsIMS[EIC_ID_IMS], binsIMS[EIC_ID_IMS] + params$mobStep)]
            EICInfoAna <- merge(EICInfoAna, tab, by = "EIC_ID_MZ", sort = FALSE)
        }
        else
            EICInfoAna[, c("mobmin", "mobmax") := 0]
        
        EICInfoAna[, EIC_ID := paste0("EIC_", .I)]
    }
    else
    {
        EICInfoAna <- data.table(mzmin = params$suspects$mz - params$mzWindow,
                                 mzmax = params$suspects$mz + params$mzWindow,
                                 EIC_ID = paste0("EIC_", seq_len(nrow(params$suspects))))
        
        # UNDONE: also put in retmin/retmax? This will affect peak finding (eg noise estimation)...
        EICInfoAna[, c("retmin", "retmax") := .(0, 0)]
        
        if (withIMS && !is.null(params$suspects[["mobility"]]))
            EICInfoAna[, c("mobmin", "mobmax") := .(params$suspects$mobility - params$IMSWindow,
                                                    params$suspects$mobility + params$IMSWindow)]
        else
            EICInfoAna[, c("mobmin", "mobmax") := .(0, 0)]
    }
    
    return(EICInfoAna)
}

#' @export
getFeaturesEICsParams <- function(methodMZ, methodIMS = NULL, ...)
{
    checkmate::assertChoice(methodMZ, c("bins", "suspects"))
    checkmate::assertChoice(methodIMS, c("bins", "suspects"), null.ok = TRUE)
    
    if (methodMZ != "suspects" && identical(methodIMS, "suspects"))
        stop("methodIMS can only be 'suspects' if methodMZ is also set to 'suspects'", call. = FALSE)
    
    ret <- list(methodMZ = methodMZ, methodIMS = methodIMS, retRange = NULL)
    
    if (methodMZ == "bins")
    {
        ret <- modifyList(ret, list(
            mzRange = c(80, 600),
            mzStep = 0.02
        ))
    }
    else if (methodMZ == "suspects")
    {
        ret <- modifyList(ret, list(
            rtWindow = 30,
            mzWindow = 0.005,
            # UNDONE these to separate param and also use elsewhere?
            skipInvalid = TRUE,
            prefCalcChemProps = TRUE,
            neutralChemProps = FALSE
        ))
    }
    
    if (!is.null(methodIMS))
    {
        if (methodIMS == "bins")
        {
            ret <- modifyList(ret, list(
                mobRange = c(0.5, 1.3),
                mobStep = 0.04
            ))
        }
        else if (methodIMS == "suspects")
        {
            ret <- modifyList(ret, list(
                IMSWindow = 0.02
            ))
        }
    }
    
    return(modifyList(ret, list(...), keep.null = TRUE))
}

#' @rdname features-class
#' @export
featuresBinning <- setClass("featuresBinning", contains = "features")

setMethod("initialize", "featuresBinning",
          function(.Object, ...) callNextMethod(.Object, algorithm = "binning", ...))

#' @export
findFeaturesBinning <- function(analysisInfo, featParams, peakParams, minIntensityIMS = 25, verbose = TRUE)
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
    
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, add = ac)
    assertFeaturesEICsParams(featParams, add = ac)
    assertFindPeakParams(peakParams, add = ac)
    checkmate::assertNumber(minIntensityIMS, lower = 0, finite = TRUE, add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    if (!is.null(featParams[["adduct"]]))
        featParams$adduct <- checkAndToAdduct(featParams$adduct)
    
    if (!is.null(featParams$suspects))
    {
        featParams$suspects <- prepareSuspectList(featParams$suspects, adduct = adduct,
                                                  skipInvalid = featParams$skipInvalid, checkDesc = TRUE,
                                                  prefCalcChemProps = featParams$prefCalcChemProps,
                                                  neutralChemProps = featParams$neutralChemProps)
        featParams$suspects <- removeDuplicateFeatsSusps(featParams$suspects, FALSE, FALSE)
        featParams$suspects <- featParams$suspects[order(mz)]
        # UNDONE: print removed suspects
    }
    
    if (is.null(featParams[["retRange"]]))
        featParams$retRange <- c(0, 0)
    
    withIMS <- !is.null(featParams[["methodIMS"]])
    
    cacheDB <- openCacheDBScope()
    baseHash <- makeHash(featParams, peakParams, minIntensityIMS)
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
        fList <- applyMSData(anaInfoTBD, needIMS = withIMS, func = function(ana, path, backend)
        {
            openMSReadBackend(backend, path)
         
            EICInfoMZ <- getFeatEICsInfo(featParams, withIMS = FALSE)
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
                EICInfoMob <- getFeatEICsInfo(featParams, withIMS = TRUE)
                ov <- foverlaps(EICInfoMob, temp, type = "within", nomatch = NULL, which = TRUE)
                EICInfo <- EICInfoMob[ov$xid]
                
                EICs <- getEICsAna(backend, EICInfo)
            }
            else
                EICInfo <- EICInfoMZ
            
            peaks <- findPeaksInEICs(EICs, peakParams, withBP = TRUE, withMobility = withIMS, cacheDB = cacheDB)

            if (featParams$methodMZ == "bins")
            {
                # only keep those peaks with m/z in the "center" of the analyzed m/z and mobility range
                peaks[, binMZStart := EICInfo[match(peaks$EIC_ID, EIC_ID)]$mzmin]
                peaks <- peaks[between(mz, binMZStart + featParams$mzStep/4, binMZStart + featParams$mzStep/4*3) == TRUE]
                if (identical(featParams[["methodIMS"]], "bins"))
                {
                    peaks[, binMobStart := EICInfo[match(peaks$EIC_ID, EIC_ID)]$mobmin]
                    peaks <- peaks[between(mobility, binMobStart + featParams$mobStep/4, binMobStart + featParams$mobStep/4*3) == TRUE]
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
