#' @include features.R
NULL

#' @rdname features-class
#' @export
featuresSuspects <- setClass("featuresSuspects", slots = c(suspects = "data.table"), contains = "features")

#' @export
findFeaturesSuspects <- function(analysisInfo, suspects, peaksParam, rtWindow = 12, mzWindow = 0.005,
                                 mobWindow = NULL, adduct = NULL, skipInvalid = TRUE, prefCalcChemProps = TRUE,
                                 neutralChemProps = FALSE, parallel = TRUE, verbose = TRUE)
{
    # UNDONE: doc that feature RT is used for checking, instead of group RT for screenSuspects()
    # UNDONE: test with large suspect lists
    # UNDONE: use adduct to set annotations? (then need to support that for sets)
    # UNDONE: force mob window option, like findMobilities()
    # UNDONE: check mobilities in prepareSuspectList()

    checkmate::assertFlag(skipInvalid) # not in assert collection, should fail before assertSuspectList
    
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, add = ac)
    assertSuspectList(suspects, needsAdduct = is.null(adduct), skipInvalid = skipInvalid, add = ac)
    assertFindPeaksParam(peaksParam, add = ac)
    aapply(checkmate::assertNumber, . ~ rtWindow + mzWindow, lower = 0, finite = TRUE, fixed = list(add = ac))
    checkmate::assertNumber(mobWindow, lower = 0, finite = TRUE, null.ok = TRUE)
    aapply(checkmate::assertFlag, . ~ skipInvalid + prefCalcChemProps + neutralChemProps + parallel + verbose,
           fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    if (!is.null(adduct))
        adduct <- checkAndToAdduct(adduct)
    
    # do this before checking cache to ensure proper errors/warnings are thrown!
    suspects <- prepareSuspectList(suspects, adduct, skipInvalid, checkDesc = TRUE,
                                   prefCalcChemProps = prefCalcChemProps, neutralChemProps = neutralChemProps)
    
    needIMS <- !is.null(mobWindow)
    if (needIMS && is.null(suspects[["mobility"]]))
    {
        warning("Suspect list does not contain mobilities: mobility data will be ignored!", call. = FALSE)
        needIMS <- FALSE
    }
    
    # UNDONE: do something like this in findMobilities()
    # if (needIMS)
    # {
    #     suspects[, order := seq_len(.N)]
    #     suspects[!is.na(mobility), ims_parent_ID := name]
    #     # expand suspect list so we can also include 'non-IMS' features
    #     suspsNoMob <- copy(suspects[!is.na(mobility)])
    #     suspsNoMob[, mobility := NA_real_]
    #     suspects[!is.na(mobility), name := appendMobToName(name, mobility)]
    #     suspects <- rbind(suspects, suspsNoMob)
    #     setorderv(suspects, c("order", "mobility"), na.last = FALSE)
    #     suspects[, order := NULL]
    # }
    
    cacheDB <- openCacheDBScope()
    baseHash <- makeHash(suspects, peaksParam, rtWindow, mzWindow, mobWindow, adduct, skipInvalid, prefCalcChemProps,
                         neutralChemProps)
    anaHashes <- getMSFileHashesFromAvailBackend(analysisInfo, needIMS = needIMS)
    anaHashes <- sapply(anaHashes, makeHash, baseHash)
    cachedData <- loadCacheData("featuresSuspects", anaHashes, simplify = FALSE, dbArg = cacheDB)
    if (!is.null(cachedData))
    {
        cachedData <- pruneList(setNames(cachedData, names(anaHashes)[match(names(cachedData), anaHashes)]))
        anaInfoTBD <- analysisInfo[!analysis %in% names(cachedData)]
    }
    else
        anaInfoTBD <- analysisInfo

    EICInfo <- data.table(mzmin = suspects$mz - mzWindow, mzmax = suspects$mz + mzWindow, retmin = 0, retmax = 0)
    if (needIMS)
        EICInfo[, c("mobmin", "mobmax") := .(suspects$mobility - mobWindow, suspects$mobility + mobWindow)]
    
    fList <- list()
    if (nrow(anaInfoTBD) > 0)
    {
        if (verbose)
            printf("Finding features of %d suspects in %d analyses ...\n", nrow(suspects), nrow(anaInfoTBD))

        # UNDONE: somehow limit RT range if suspect rt is given?
        # UNDONE: make withBP configurable?
        printf("Loading EICs...\n")
        allEICs <- doGetEICs(anaInfoTBD, rep(list(EICInfo), nrow(anaInfoTBD)), compress = FALSE, cacheDB = cacheDB)
        allEICs <- lapply(allEICs, setNames, suspects$name)
        
        fList <- findPeaksInEICs(allEICs, peaksParam, withBP = FALSE, parallel = parallel, cacheDB = cacheDB)
        fList <- lapply(fList, function(peaks)
        {
            peaks <- copy(peaks)
            setnames(peaks, "EIC_ID", "suspect")
            
            if (!is.null(suspects[["rt"]]))
            {
                peaks[, susp_rt := suspects[match(suspect, name)]$rt]
                peaks <- peaks[is.na(susp_rt) | numLTE(abs(ret - susp_rt), rtWindow)][, susp_rt := NULL]
            }
            
            if (needIMS)
                peaks[, ims_parent_ID := NA_character_]
            else
                peaks <- removeDTColumnsIfPresent(peaks, c("mobmin", "mobmax", "mobility"))
            
            return(peaks[])
        })
        
        for (a in anaInfoTBD$analysis)
            saveCacheData("featuresSuspects", fList[[a]], anaHashes[[a]])
    }
    
    if (!is.null(cachedData))
    {
        fList <- c(fList, cachedData)
        fList <- fList[analysisInfo$analysis] # put original order
    }
    
    if (verbose)
    {
        printf("Done!\n")
        printFeatStats(fList)
        printf("Found %d out of %d suspects\n", uniqueN(unlist(lapply(fList, "[[", "suspect"))), nrow(suspects))
    }
    
    return(featuresSuspects(analysisInfo = analysisInfo, features = fList, suspects = suspects,
                            algorithm = paste0("suspects-", peaksParam$algorithm)))
}
