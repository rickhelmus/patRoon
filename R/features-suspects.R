#' @include features.R
NULL

#' @rdname features-class
#' @export
featuresSuspects <- setClass("featuresSuspects", slots = c(suspects = "data.table"), contains = "features")

#' @export
findFeaturesSuspects <- function(analysisInfo, suspects, findPeaksAlgo, rtWindow = 12, mzWindow = 0.005, adduct = NULL,
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
    assertFindPeaksAlgo(findPeaksAlgo, add = ac)
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
    allEICs <- lapply(allEICs, setNames, suspects$name)
    
    # UNDONE: somehow limit RT range if suspect rt is given?
    
    fList <- findPeaksInEICs(allEICs, findPeaksAlgo, ..., parallel = parallel, cacheDB = cacheDB)
    fList <- lapply(fList, function(peaks)
    {
        peaks <- copy(peaks)
        setnames(peaks, "EIC_ID", "suspect")
        
        if (!is.null(suspects[["rt"]]))
        {
            peaks[, susp_rt := suspects[match(suspect, name)]$rt]
            peaks <- peaks[is.na(susp_rt) | numLTE(abs(ret - susp_rt), rtWindow)][, susp_rt := NULL]
        }
        
        return(peaks)
    })
    
    if (verbose)
    {
        printf("Done!\n")
        printFeatStats(fList)
    }
    
    return(featuresSuspects(analysisInfo = analysisInfo, features = fList, suspects = suspects,
                            algorithm = paste0("suspects-", findPeaksAlgo)))
}
