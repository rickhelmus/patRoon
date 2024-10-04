#' @include features.R
NULL

#' @rdname features-class
#' @export
featuresSuspects <- setClass("featuresSuspects", slots = c(suspects = "data.table"), contains = "features")

#' @export
findfeaturesSuspects <- function(analysisInfo, suspects, findPeaksAlgo, rtWindow = 12, mzWindow = 0.005, adduct = NULL,
                                 skipInvalid = TRUE, prefCalcChemProps = TRUE, neutralChemProps = FALSE, ...,
                                 parallel = TRUE, verbose = TRUE)
{
    # UNDONE: doc that feature RT is used for checking, instead of group RT for screenSuspects()
    # UNDONE: test with large suspect lists
    # UNDONE: use adduct to set annotations?

    checkmate::assertFlag(skipInvalid) # not in assert collection, should fail before assertSuspectList
    
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, add = ac)
    assertSuspectList(suspects, needsAdduct = is.null(adduct), skipInvalid = skipInvalid, add = ac)
    checkmate::assertString(findPeaksAlgo, min.chars = 1, add = ac) # UNDONE check algo choice if findPeaks() will not be exported
    aapply(checkmate::assertNumber, . ~ rtWindow + mzWindow, lower = 0, finite = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertFlag, . ~ skipInvalid + prefCalcChemProps + neutralChemProps + parallel + verbose,
           fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    cacheDB <- openCacheDBScope()
    anaCount <- nrow(analysisInfo)
    
    if (!is.null(adduct))
        adduct <- checkAndToAdduct(adduct)
    
    # do this before checking cache to ensure proper errors/warnings are thrown!
    suspects <- prepareSuspectList(suspects, adduct, skipInvalid, checkDesc = TRUE,
                                   prefCalcChemProps = prefCalcChemProps, neutralChemProps = neutralChemProps)
    
    if (verbose)
        printf("Finding features from suspects for %d analyses ...\n", anaCount)

    # UNDONE: fill in IMS ranges (if available)
    EICInfoList <- rep(list(data.table(mzmin = suspects$mz - mzWindow, mzmax = suspects$mz + mzWindow,
                                       retmin = 0, retmax = 0)), nrow(analysisInfo))
    printf("Loading EICs...\n")
    allEICs <- doGetEICs(analysisInfo, EICInfoList, compress = FALSE, cacheDB = cacheDB)

    # NOTE: EICs are cached in doGetEICs(), so here we mainly do caching for peaks
    baseHash <- makeHash(suspects, findPeaksAlgo, rtWindow, mzWindow, adduct, skipInvalid, prefCalcChemProps,
                         neutralChemProps, list(...))
    
    printf("Finding peaks...\n")
    fList <- doApply("sapply", parallel, allEICs, function(EICs)
    {
        hash <- makeHash(baseHash, EICs)
        cd <- loadCacheData("featuresSuspects", hash, cacheDB)
        if (!is.null(cd))
        {
            doProgress()
            return(cd)
        }
        
        names(EICs) <- suspects$name
        EICs <- lapply(EICs, setDT)
        
        # UNDONE: limit RT range if suspect rt is given?
        
        peaks <- findPeaks(EICs, findPeaksAlgo, ..., verbose = FALSE)
        peaks <- rbindlist(peaks, idcol = "suspect")
        
        if (!is.null(suspects[["rt"]]))
        {
            peaks[, susp_rt := suspects[match(suspect, name)]$rt]
            peaks <- peaks[is.na(susp_rt) | numLTE(abs(ret - susp_rt), rtWindow)][, susp_rt := NULL]
        }
        
        peaks[, c("mzmin", "mzmax", "mz", "mobmin", "mobmax", "mobility") := {
            eic <- EICs[[suspect]][intensity != 0 & time %between% c(retmin, retmax)]
            if (is.null(eic[["mobility"]]))
                eic[, mobility := NA_real_]
            list(min(eic$mz), max(eic$mz), weighted.mean(eic$mz, eic$intensity),
                 min(eic$mobility), max(eic$mobility), weighted.mean(eic$mobility, eic$intensity))
        }, by = seq_len(nrow(peaks))]

        # NOTE: we could also set mobilities after checking if data is available, but then we need to repeat the EIC subsetting above
        if (is.null(EICs$mobility))
            peaks[, c("mobmin", "mobmax", "mobility") := NULL]

        # make unique IDs
        peaks[, ID := make.unique(suspect)]
        
        saveCacheData("featuresSuspects", peaks, hash, cacheDB)
        
        doProgress()
        
        return(peaks)
    }, simplify = FALSE)
    
    if (verbose)
    {
        printf("Done!\n")
        printFeatStats(fList)
    }
    
    return(featuresSuspects(analysisInfo = analysisInfo, features = fList, suspects = suspects,
                            algorithm = paste0("suspects-", findPeaksAlgo)))
}
