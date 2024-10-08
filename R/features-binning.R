#' @include features.R
NULL

#' @rdname features-class
#' @export
featuresBinning <- setClass("featuresBinning", contains = "features")

setMethod("initialize", "featuresBinning",
          function(.Object, ...) callNextMethod(.Object, algorithm = "binning", ...))

#' @export
findFeaturesBinning <- function(analysisInfo, mzRange = c(50, 400), mzStep = 0.02, findPeaksAlgo, ..., parallel = TRUE,
                                verbose = TRUE)
{
    # UNDONE: add refs to docs, and highlight changes
    
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, add = ac)
    assertRange(mzRange, add = ac)
    checkmate::assertNumber(mzStep, lower = 0.000001, finite = TRUE, add = ac)
    assertFindPeaksAlgo(findPeaksAlgo, add = ac)
    aapply(checkmate::assertFlag, . ~ parallel + verbose, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    cacheDB <- openCacheDBScope()
    
    bins <- seq(mzRange[1], mzRange[2], by = mzStep * 0.5)
    names(bins) <- paste0("bin_M", bins)
    
    if (verbose)
        printf("Finding features in %d bins for %d analyses ...\n", length(bins), nrow(analysisInfo))

    printf("Loading EICs...\n")
    EICInfoList <- rep(list(data.table(mzmin = bins, mzmax = bins + mzStep, retmin = 0, retmax = 0)), nrow(analysisInfo))
    allEICs <- doGetEICs(analysisInfo, EICInfoList, compress = FALSE, withBP = TRUE, cacheDB = cacheDB)
    allEICs <- lapply(allEICs, setNames, names(bins))
    
    fList <- findPeaksInEICs(allEICs, findPeaksAlgo, withBP = TRUE, ..., parallel = parallel, cacheDB = cacheDB)
    
    fList <- lapply(fList, function(fTab)
    {
        # only keep those peaks with m/z in the "center" of the analyzed m/z range
        fTab[between(mz, bins[EIC_ID] + mzStep/4, bins[EIC_ID] + mzStep/4*3) == TRUE][, EIC_ID := NULL]
    })
    
    # UNDONE: use BP intensity?
        
    if (verbose)
    {
        printf("Done!\n")
        printFeatStats(fList)
    }
    
    return(featuresBinning(analysisInfo = analysisInfo, features = fList))
}
