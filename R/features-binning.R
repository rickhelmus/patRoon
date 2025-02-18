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

getFeatEICsInfo <- function(params, withIMS, MS2Info)
{
    EICInfo <- NULL
    if (params$methodMZ == "bins")
    {
        # UNDONE: also support other binning approaches?
        binsMZ <- seq(params$mzRange[1], params$mzRange[2], by = params$mzStep * 0.5)
        names(binsMZ) <- paste0("bin_M", binsMZ)
        
        EICInfo <- data.table(mzmin = binsMZ, mzmax = binsMZ + params$mzStep, EIC_ID_MZ = names(binsMZ))
    }
    else if (params$methodMZ == "suspects")
        EICInfo <- data.table(mzmin = params$suspects$mz - params$mzWindow, mzmax = params$suspects$mz + params$mzWindow)
    else if (params$methodMZ == "ms2")
        EICInfo <- data.table(mzmin = MS2Info$mz - params$mzWindow, mzmax = MS2Info$mz + params$mzWindow)
    
    if (withIMS)
    {
        if (params$methodIMS == "bins")
        {
            binsIMS <- NULL
            {
                binsIMS <- seq(params$mobRange[1], params$mobRange[2], by = params$mobStep * 0.5)
                names(binsIMS) <- paste0("bin_I", binsIMS)
            }
            tab <- CJ(EIC_ID_MZ = names(binsMZ), EIC_ID_IMS = names(binsIMS), sorted = FALSE)
            tab[, c("mobmin", "mobmax") := .(binsIMS[EIC_ID_IMS], binsIMS[EIC_ID_IMS] + params$mobStep)]
            EICInfo <- merge(EICInfo, tab, by = "EIC_ID_MZ", sort = FALSE)
        }
        else if (params$methodIMS == "suspects")
        {
            EICInfo <- EICInfo[, c("mobmin", "mobmax") := .(params$suspects$mobility - params$IMSWindow,
                                                            params$suspects$mobility + params$IMSWindow)]
        }
        else if (params$methodIMS == "ms2")
        {
            EICInfo <- EICInfo[, c("mobmin", "mobmax") := .(MS2Info$mobility - params$IMSWindow,
                                                            MS2Info$mobility + params$IMSWindow)]
        }
    }
    else
        EICInfo[, c("mobmin", "mobmax") := 0]
    
    EICInfo[, c("retmin", "retmax") := .(params$retRange[1], params$retRange[2])]
    EICInfo[, EIC_ID := paste0("EIC_", .I)]
    
    return(EICInfo)
}

#' @export
getFeaturesEICsParams <- function(methodMZ, methodIMS = NULL, ...)
{
    checkmate::assertChoice(methodMZ, c("bins", "suspects", "ms2"))
    checkmate::assertChoice(methodIMS, c("bins", "suspects", "ms2"), null.ok = TRUE)
    
    if (methodMZ != "suspects" && identical(methodIMS, "suspects"))
        stop("methodIMS can only be 'suspects' if methodMZ is also set to 'suspects'", call. = FALSE)
    
    ret <- list(methodMZ = methodMZ, methodIMS = methodIMS, retRange = NULL, minEICIntensity = 5000,
                minEICAdjTime = 5, minEICAdjPoints = 5, minEICAdjIntensity = 250)
    
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
    else if (methodMZ == "ms2")
    {
        ret <- modifyList(ret, list(
            rtWindow = 3,
            mzWindow = 0.005,
            minTIC = 10000,
            clusterMethod = "distance"
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
        else if (methodIMS == "ms2")
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
    
    if (featParams$methodMZ == "suspects")
    {
        featParams$suspects <- prepareSuspectList(featParams$suspects, adduct = adduct,
                                                  skipInvalid = featParams$skipInvalid, checkDesc = TRUE,
                                                  prefCalcChemProps = featParams$prefCalcChemProps,
                                                  neutralChemProps = featParams$neutralChemProps)
        featParams$suspectsOrig <- featParams$suspects
        featParams$suspects <- removeDuplicateFeatsSusps(featParams$suspects, FALSE, FALSE)
        featParams$suspects <- featParams$suspects[order(mz)]
        # UNDONE: print removed suspects?
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
        ret <- doGetEICsForAna(backend, EICInfo$mzmin, EICInfo$mzmax, EICInfo$retmin, EICInfo$retmax,
                               EICInfo$mobmin, EICInfo$mobmax, mzExpIMSWindow = 0, minIntensityIMS = minIntensityIMS,
                               compress = FALSE, showProgress = FALSE, withBP = TRUE,
                               minEICIntensity = featParams$minEICIntensity, minEICAdjTime = featParams$minEICAdjTime,
                               minEICAdjPoints = featParams$minEICAdjPoints,
                               minEICAdjIntensity = featParams$minEICAdjIntensity)
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
         
            MS2Info <- NULL
            if (featParams$methodMZ == "ms2")
            {
                MS2Info <- if (identical(featParams$methodIMS, "ms2"))
                    getIsolationMZsAndMobs(backend, featParams$clusterMethod, featParams$mzWindow, featParams$IMSWindow,
                                           featParams$minTIC)
                else
                    getIsolationMZs(backend, featParams$clusterMethod, featParams$mzWindow, featParams$minTIC)
                setDT(MS2Info)
            }
            
            EICInfoMZ <- getFeatEICsInfo(featParams, withIMS = FALSE, MS2Info = MS2Info)
            EICs <- getEICsAna(backend, EICInfoMZ)
            # omit missing
            EICInfoMZ <- EICInfoMZ[EIC_ID %chin% names(EICs)]
            
            EICInfo <- NULL
            if (withIMS)
            {
                # With IMS worksflows the mz EICs are only used as a pre-filter. As the EIC object is potentially large we
                # remove it here.
                rm(EICs) 
                gc()
                
                # remove complete m/z bins that were filtered out before
                temp <- EICInfoMZ[, c("mzmin", "mzmax"), with = FALSE]
                setkeyv(temp, c("mzmin", "mzmax"))
                EICInfoMob <- getFeatEICsInfo(featParams, withIMS = TRUE, MS2Info = MS2Info)
                ov <- foverlaps(EICInfoMob, temp, type = "within", nomatch = NULL, which = TRUE)
                EICInfo <- EICInfoMob[ov$xid]
                
                EICs <- getEICsAna(backend, EICInfo)
            }
            else
                EICInfo <- EICInfoMZ
            
            peaks <- findPeaksInEICs(EICs, peakParams, withBP = TRUE, withMobility = withIMS,
                                     logPath = file.path("log", "featEICs", paste0(ana, ".txt")), cacheDB = cacheDB)

            # only keep those peaks with m/z in the "center" of the analyzed m/z and mobility range
            if (featParams$methodMZ == "bins")
            {
                peaks[, binMZStart := EICInfo[match(peaks$EIC_ID, EIC_ID)]$mzmin]
                peaks <- peaks[between(mz, binMZStart + featParams$mzStep/4, binMZStart + featParams$mzStep/4*3) == TRUE]
            }
            if (identical(featParams[["methodIMS"]], "bins"))
            {
                peaks[, binMobStart := EICInfo[match(peaks$EIC_ID, EIC_ID)]$mobmin]
                peaks <- peaks[between(mobility, binMobStart + featParams$mobStep/4, binMobStart + featParams$mobStep/4*3) == TRUE]
            }
            if (featParams$methodMZ == "suspects" && is.finite(featParams$rtWindow) && !is.null(featParams$suspects[["rt"]]))
            {
                # only keep peaks with at least one closely eluting suspect
                peaks[, keep := {
                    susp <- featParams$suspectsOrig[numLTE(abs(mz - omz), featParams$mzWindow) &
                                                        numLTE(abs(rt - ort), featParams$rtWindow),
                                                    env = list(omz = mz, ort = ret)]
                    if (identical(featParams[["methodIMS"]], "ms2") && !is.null(featParams$suspects[["mobility"]]))
                        susp <- susp[numLTE(abs(mobility - omob), featParams$IMSWindow), env = list(omob = mobility)]
                    nrow(susp) > 0
                }, by = .I]
                peaks <- peaks[keep == TRUE]
            }
            else if (featParams$methodMZ %in% "ms2" && is.finite(featParams$rtWindow))
            {
                # only keep peaks that elute closely to at least one MS2 spectrum
                peaks[, keep := {
                    MS2Peaks <- MS2Info[numLTE(abs(mz - omz), featParams$mzWindow), env = list(omz = mz)]
                    if (identical(featParams[["methodIMS"]], "ms2"))
                        MS2Peaks <- MS2Peaks[numLTE(abs(mobility - omob), featParams$IMSWindow), env = list(omob = mobility)]
                    allRet <- unlist(MS2Peaks$times)
                    any(numLTE(abs(ret - allRet), featParams$rtWindow))
                }, by = .I]
                peaks <- peaks[keep == TRUE]
            }
            
            peaks <- removeDuplicateFeatsSusps(peaks, TRUE, TRUE)
            peaks <- removeDTColumnsIfPresent(peaks, c("binMZStart", "binMobStart", "keep"))
            
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
