#' @include features.R
NULL

removeDuplicateFeatsSusps <- function(tab, checkRet, selectTopIntens, nameCol)
{
    # marks duplicates in feature or suspect table, and keeps the feature with the highest intensity

    if (nrow(tab) <= 1)
        return(tab)

    isSame <- function(x, y)
    {
        yes <- numLTE(abs(x$mz - y$mz), defaultLim("mz", "very_narrow"))
        if (checkRet)
            yes <- yes & ((is.na(x$ret) & is.na(y$ret)) | numLTE(abs(x$ret - y$ret), defaultLim("retention", "very_narrow")))
        if (!is.null(x[["mobility"]]))
            yes <- yes & numLTE(abs(x$mobility - y$mobility), defaultLim("mobility", "very_narrow"))
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
        nextTab <- nextTab[duplicate == FALSE & isSame(rTab, nextTab)]
        if (nrow(nextTab) > 0)
        {
            if (selectTopIntens)
            {
                maxInt <- max(rTab$intensity, nextTab$intensity)
                tab[get(nameCol) %chin% c(rTab[[nameCol]], nextTab[[nameCol]]), duplicate := !numEQ(intensity, maxInt)]
            }
            else
            {
                # otherwise just keep current and mark all from next rows
                tab[get(nameCol) %chin% nextTab[[nameCol]], duplicate := TRUE]
            }
        }
    }
    return(tab[duplicate == FALSE][, duplicate := NULL])
}

getPiekEICsInfo <- function(params, suspects, withIMS, MS2Info)
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
        EICInfo <- data.table(mzmin = suspects$mz - params$mzWindow, mzmax = suspects$mz + params$mzWindow)
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
            EICInfo <- EICInfo[, c("mobmin", "mobmax") := .(suspects$mobility - params$IMSWindow,
                                                            suspects$mobility + params$IMSWindow)]
        }
        else if (params$methodIMS == "ms2")
        {
            EICInfo <- EICInfo[, c("mobmin", "mobmax") := .(MS2Info$mobility - params$IMSWindow,
                                                            MS2Info$mobility + params$IMSWindow)]
        }
    }
    else
        EICInfo[, c("mobmin", "mobmax") := 0]
    
    EICInfo[, EIC_ID := paste0("EIC_", .I)]
    
    return(EICInfo)
}

#' @export
getPiekGenEICParams <- function(methodMZ, methodIMS = NULL, ...)
{
    checkmate::assertChoice(methodMZ, c("bins", "suspects", "ms2"))
    checkmate::assertChoice(methodIMS, c("bins", "suspects", "ms2"), null.ok = TRUE)
    
    if (methodMZ != "suspects" && identical(methodIMS, "suspects"))
        stop("methodIMS can only be 'suspects' if methodMZ is also set to 'suspects'", call. = FALSE)
    
    ret <- list(methodMZ = methodMZ, methodIMS = methodIMS, retRange = NULL, minEICIntensity = 5000,
                minEICAdjTime = 5, minEICAdjPoints = 5, minEICAdjIntensity = 250, topMostEIC = 10000,
                topMostEICPre = 10000)
    
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
            rtWindow = defaultLim("retention", "medium"),
            mzWindow = defaultLim("mz", "medium"),
            # UNDONE these to separate param and also use elsewhere?
            skipInvalid = TRUE,
            prefCalcChemProps = TRUE,
            neutralChemProps = FALSE
        ))
    }
    else if (methodMZ == "ms2")
    {
        ret <- modifyList(ret, list(
            rtWindow = defaultLim("retention", "very_narrow"),
            mzWindow = defaultLim("mz", "narrow"),
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
                IMSWindow = defaultLim("mobility", "medium")
            ))
        }
        else if (methodIMS == "ms2")
        {
            ret <- modifyList(ret, list(
                IMSWindow = defaultLim("mobility", "medium") # UNDONE: or narrow?
            ))
        }
    }
    
    return(modifyList(ret, list(...), keep.null = TRUE))
}

#' @rdname features-class
#' @export
featuresPiek <- setClass("featuresPiek", contains = "features")

setMethod("initialize", "featuresPiek",
          function(.Object, ...) callNextMethod(.Object, algorithm = "piek", ...))

#' @template minIntensityIMS-arg
#' @export
findFeaturesPiek <- function(analysisInfo, genEICParams, peakParams, suspects = NULL, adduct = NULL,
                             minIntensityIMS = 25, verbose = TRUE)
{
    # UNDONE: add refs to docs, and highlight changes
    # UNDONE: use BP intensity?
    # UNDONE: test empties, eg no EICs
    
    if (!is.null(adduct))
        adduct <- checkAndToAdduct(adduct)
    
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, add = ac)
    assertPiekGenEICParams(genEICParams, add = ac)
    assertFindPeakParams(peakParams, add = ac)
    assertSuspectList(suspects, needsAdduct = is.null(adduct), skipInvalid = genEICParams$skipInvalid, null.ok = TRUE,
                      add = ac)
    checkmate::assertNumber(minIntensityIMS, lower = 0, finite = TRUE, add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    maybePrintf <- \(...) if (verbose) printf(...)
    
    if (genEICParams$methodMZ == "suspects")
    {
        suspects <- prepareSuspectList(suspects, adduct = adduct, skipInvalid = genEICParams$skipInvalid,
                                       checkDesc = TRUE, prefCalcChemProps = genEICParams$prefCalcChemProps,
                                       neutralChemProps = genEICParams$neutralChemProps)
        suspectsOrig <- suspects
        suspects <- removeDuplicateFeatsSusps(suspects, FALSE, FALSE, "name")
        suspects <- suspects[order(mz)]
        suspsRM <- setdiff(suspectsOrig$name, suspects$name)
        maybePrintf("The following %d non-unique suspects were removed %s.\n", length(suspsRM),
                    getStrListWithMax(suspsRM, 10, ", "))
    }
    
    if (is.null(genEICParams[["retRange"]]))
        genEICParams$retRange <- c(0, 0)
    
    withIMS <- !is.null(genEICParams[["methodIMS"]])
    
    cacheDB <- openCacheDBScope()
    baseHash <- makeHash(genEICParams, peakParams, minIntensityIMS)
    anaHashes <- getMSFileHashesFromAvailBackend(analysisInfo, needIMS = withIMS)
    anaHashes <- sapply(anaHashes, makeHash, baseHash)
    cachedData <- pruneList(loadCacheData("featuresPiek", anaHashes, simplify = FALSE, dbArg = cacheDB))
    if (length(cachedData) > 0)
    {
        names(cachedData) <- names(anaHashes)[match(names(cachedData), anaHashes)]
        anaInfoTBD <- analysisInfo[!analysis %in% names(cachedData)]
    }
    else
        anaInfoTBD <- analysisInfo
    
    getEICsAna <- function(backend, EICInfo, mode, topMost)
    {
        args <- list(backend, EICInfo$mzmin, EICInfo$mzmax, genEICParams$retRange[1], genEICParams$retRange[2],
                     EICInfo$mobmin, EICInfo$mobmax, mzExpIMSWindow = 0, minIntensityIMS = minIntensityIMS,
                     mode = mode, minEICIntensity = genEICParams$minEICIntensity,
                     minEICAdjTime = genEICParams$minEICAdjTime,
                     minEICAdjPoints = genEICParams$minEICAdjPoints,
                     minEICAdjIntensity = genEICParams$minEICAdjIntensity, topMost = topMost)
        ret <- do.call(if (mode == "test") getEICList else doGetEICsForAna, args)
        names(ret) <- EICInfo$EIC_ID
        if (mode != "test")
            ret <- pruneList(ret, checkZeroRows = TRUE, keepAttr = TRUE)
        return(ret)
    }
    
    fList <- list()
    if (nrow(anaInfoTBD) > 0)
    {
        fList <- applyMSData(anaInfoTBD, needIMS = withIMS, showProgress = FALSE, func = function(ana, path, backend)
        {
            maybePrintf("\n-------\nFinding features for '%s'...\n", ana)
            
            openMSReadBackend(backend, path)
         
            MS2Info <- NULL
            if (genEICParams$methodMZ == "ms2")
            {
                MS2Info <- if (identical(genEICParams$methodIMS, "ms2"))
                    getIsolationMZsAndMobs(backend, genEICParams$clusterMethod, genEICParams$mzWindow,
                                           genEICParams$IMSWindow, genEICParams$minTIC)
                else
                    getIsolationMZs(backend, genEICParams$clusterMethod, genEICParams$mzWindow, genEICParams$minTIC)
                setDT(MS2Info)
            }
            
            EICs <- EICInfo <- NULL
            EICInfoMZ <- getPiekEICsInfo(genEICParams, suspects, withIMS = FALSE, MS2Info = MS2Info)
            if (withIMS)
            {
                maybePrintf("Pre-checking %d m/z EICs... ", nrow(EICInfoMZ))
                testEICs <- getEICsAna(backend, EICInfoMZ, "test", genEICParams$topMostEICPre)
                testEICs <- unlist(testEICs)
                EICInfoMZ <- EICInfoMZ[EIC_ID %chin% names(testEICs)[testEICs]]
                maybePrintf("Done! Eliminated %d (%.2f%%) EICs\n", sum(!testEICs),
                            (sum(!testEICs) / length(testEICs)) * 100)
                
                # remove complete m/z bins that were filtered out before
                temp <- EICInfoMZ[, c("mzmin", "mzmax"), with = FALSE]
                setkeyv(temp, c("mzmin", "mzmax"))
                EICInfoMob <- getPiekEICsInfo(genEICParams, suspects, withIMS = TRUE, MS2Info = MS2Info)
                ov <- foverlaps(EICInfoMob, temp, type = "within", nomatch = NULL, which = TRUE)
                EICInfo <- EICInfoMob[ov$xid]
                
                maybePrintf("Loading %d m/z+mobility EICs... ", nrow(EICInfo))
                EICs <- getEICsAna(backend, EICInfo, "full", genEICParams$topMostEIC)
                maybePrintf("Done!\n")
            }
            else
            {
                maybePrintf("Loading %d m/z EICs... ", nrow(EICInfoMZ))
                EICs <- getEICsAna(backend, EICInfoMZ, "full_mz", genEICParams$topMostEIC)
                EICInfo <- EICInfoMZ[EIC_ID %chin% names(EICs)] # omit missing
                maybePrintf("Done!\n")
            }
            
            maybePrintf("Finding peaks in %d remaining EICs... ", length(EICs))
            peaks <- findPeaksInEICs(EICs, peakParams, withMobility = withIMS,
                                     logPath = file.path("log", "featEICs", paste0(ana, ".txt")), cacheDB = cacheDB)
            maybePrintf("Done! Found %d peaks.\n", nrow(peaks))

            # only keep those peaks with m/z in the "center" of the analyzed m/z and mobility range
            if (genEICParams$methodMZ == "bins")
            {
                peaks[, binMZStart := EICInfo[match(peaks$EIC_ID, EIC_ID)]$mzmin]
                peaks <- peaks[between(mz, binMZStart + genEICParams$mzStep/4, binMZStart + genEICParams$mzStep/4*3) == TRUE]
            }
            if (identical(genEICParams[["methodIMS"]], "bins"))
            {
                peaks[, binMobStart := EICInfo[match(peaks$EIC_ID, EIC_ID)]$mobmin]
                peaks <- peaks[between(mobility, binMobStart + genEICParams$mobStep/4, binMobStart + genEICParams$mobStep/4*3) == TRUE]
            }
            if (genEICParams$methodMZ == "suspects" && is.finite(genEICParams$rtWindow) && !is.null(suspects[["rt"]]))
            {
                # only keep peaks with at least one closely eluting suspect
                peaks[, keep := {
                    susp <- suspectsOrig[numLTE(abs(mz - omz), genEICParams$mzWindow) &
                                             numLTE(abs(rt - ort), genEICParams$rtWindow),
                                         env = list(omz = mz, ort = ret)]
                    if (identical(genEICParams[["methodIMS"]], "ms2") && !is.null(suspects[["mobility"]]))
                        susp <- susp[numLTE(abs(mobility - omob), genEICParams$IMSWindow), env = list(omob = mobility)]
                    nrow(susp) > 0
                }, by = .I]
                peaks <- peaks[keep == TRUE]
            }
            else if (genEICParams$methodMZ %in% "ms2" && is.finite(genEICParams$rtWindow))
            {
                # only keep peaks that elute closely to at least one MS2 spectrum
                peaks[, keep := {
                    MS2Peaks <- MS2Info[numLTE(abs(mz - omz), genEICParams$mzWindow), env = list(omz = mz)]
                    if (identical(genEICParams[["methodIMS"]], "ms2"))
                        MS2Peaks <- MS2Peaks[numLTE(abs(mobility - omob), genEICParams$IMSWindow), env = list(omob = mobility)]
                    allRet <- unlist(MS2Peaks$times)
                    any(numLTE(abs(ret - allRet), genEICParams$rtWindow))
                }, by = .I]
                peaks <- peaks[keep == TRUE]
            }
            
            peaks <- removeDuplicateFeatsSusps(peaks, TRUE, TRUE, "ID")
            peaks <- removeDTColumnsIfPresent(peaks, c("binMZStart", "binMobStart", "keep"))

            maybePrintf("%d peaks remain after filtering.\n", nrow(peaks))            
            
            return(peaks)
        })
        
        for (a in anaInfoTBD$analysis)
            saveCacheData("featuresPiek", fList[[a]], anaHashes[[a]], dbArg = cacheDB)
    }
    
    if (length(cachedData) > 0)
    {
        fList <- c(fList, cachedData)
        fList <- fList[analysisInfo$analysis] # put original order
    }
    
    if (verbose)
    {
        printf("\n-------\nDone!\n")
        printFeatStats(fList)
    }
    
    return(featuresPiek(analysisInfo = analysisInfo, features = fList, hasMobilities = withIMS))
}
