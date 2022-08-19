#' @include features.R
NULL

#' @rdname features-class
#' @export
featuresSuspects <- setClass("featuresSuspects", contains = "features")

#' @export
findfeaturesSuspects <- function(analysisInfo, suspects, findPeaksAlgo, rtWindow = 12, mzWindow = 0.005, adduct = NULL,
                                 skipInvalid = TRUE, prefCalcChemProps = TRUE, ..., parallel = TRUE, verbose = TRUE)
{
    # UNDONE: doc that feature RT is used for checking, instead of group RT for screenSuspects()
    # UNDONE: test with large suspect lists
    # UNDONE: export? or just internal function?
    # UNDONE: store suspect list in slot? Probably yes, to make screenInfo
    
    checkmate::assertFlag(skipInvalid) # not in assert collection, should fail before assertSuspectList
    
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, verifyCentroided = TRUE, add = ac)
    assertSuspectList(suspects, needsAdduct = is.null(adduct), skipInvalid = skipInvalid, add = ac)
    checkmate::assertString(findPeaksAlgo, min.chars = 1, add = ac) # UNDONE check algo choice if findPeaks() will not be exported
    aapply(checkmate::assertNumber, . ~ rtWindow + mzWindow, lower = 0, finite = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertFlag, . ~ skipInvalid + prefCalcChemProps + parallel + verbose, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    anaCount <- nrow(analysisInfo)
    
    if (!is.null(adduct))
        adduct <- checkAndToAdduct(adduct)
    
    # do this before checking cache to ensure proper errors/warnings are thrown!
    suspects <- prepareSuspectList(suspects, adduct, skipInvalid, checkDesc = TRUE,
                                   prefCalcChemProps = prefCalcChemProps)
    
    if (verbose)
        printf("Finding features from suspects for %d analyses ...\n", anaCount)
    
    doFindFeats <- function(ana, path)
    {
        # UNDONE: also support regular mzR interface
        fp <- getBrukerAnalysisPath(ana, path)
        
        TIMSDB <- openTIMSMetaDBScope(f = fp)
        frames <- getTIMSMetaTable(TIMSDB, "Frames", c("Id", "Time", "MsMsType"))
        frames <- frames[MsMsType == 0]
        
        # UNDONE: fill in IMS ranges (if available)
        EICs <- getTIMSEICs(fp, frames$Id, suspects$mz - mzWindow, suspects$mz + mzWindow,
                            rep(0, nrow(suspects)), rep(0, nrow(suspects)), FALSE)
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
        
        peaks[, mz := suspects[match(suspect, name)]$mz]
        peaks[, c("mzmin", "mzmax") := .(mz - mzWindow, mz + mzWindow)]
        
        # UNDONE: also assign mobility
        
        # make unique IDs
        peaks[, ID := make.unique(suspect)]
        
        return(peaks)
    }
    
    fList <- if (parallel)
        withProg(nrow(analysisInfo), TRUE, future.apply::future_Map(analysisInfo$analysis, analysisInfo$path,
                                                                    f = doFindFeats))
    else
        withProg(nrow(analysisInfo), FALSE, Map(analysisInfo$analysis, analysisInfo$path, f = doFindFeats))
    
    
    if (verbose)
    {
        printf("Done!\n")
        printFeatStats(fList)
    }
    
    return(featuresSuspects(analysisInfo = analysisInfo, features = fList, algorithm = paste0("suspects-", findPeaksAlgo)))
}
