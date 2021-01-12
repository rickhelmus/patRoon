#' @include main.R
#' @include components-features.R
NULL

#' @export
componentsCliqueMS <- setClass("componentsCliqueMS", contains = "componentsFeatures")

setMethod("initialize", "componentsCliqueMS",
          function(.Object, ...) callNextMethod(.Object, algorithm = "openms", ...))

#' @rdname component-generation
#' @export
setMethod("generateComponentsCliqueMS", "featureGroups", function(fGroups, ionization,
                                                                  mzWindow = 0.005, minSize = 2,
                                                                  relMinAdductAbundance = 0.75,
                                                                  extraOpts = NULL)
{
    # UNDONE: all features are currently annotated (ie including not in a group), should be fine once featng is merged
    # UNDONE: more parameters and proper defaults
    # UNDONE: parallel?
    
    checkPackage("cliqueMS", "https://github.com/rickhelmus/cliqueMS") # UNDONE
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertChoice(ionization, c("positive", "negative"), add = ac)
    checkmate::assertCount(minSize, positive = TRUE, add = ac)
    aapply(checkmate::assertNumber, . ~ mzWindow + relMinAdductAbundance, finite = TRUE, lower = 0,
           fixed = list(add = ac))
    checkmate::assertList(extraOpts, any.missing = FALSE, names = "unique", null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    if (length(fGroups) == 0)
        return(componentsCliqueMS())
        
    anas <- analyses(fGroups)
    fList <- getFeatures(fGroups)
    
    printf("Annotating all features with CliqueMS for %d analyses ...\n", nrow(anaInfo))
    
    hash <- makeHash(fList, ionization, mzWindow, minSize, relMinAdductAbundance, extraOpts)
    
    # UNDONE: save/load from cache
    
    if (ionization == "positive")
        data(positive.adinfo, package = "cliqueMS", envir = environment())
    else
        data(negative.adinfo, package = "cliqueMS", envir = environment())
    
    prog <- openProgBar(0, length(anas))
    
    featComponents <- setNames(lapply(seq_along(anas), function(i)
    {
        xdata <- getXCMSnExp(fList[i], verbose = FALSE)
        
        cliques <- cliqueMS::getCliques(xdata, filter = TRUE)
        cliques <- cliqueMS::getIsotopes(cliques, ppm = 10)
        cliques <- cliqueMS::getAnnotation(cliques, ppm = 10,
                                           adinfo = if (ionization == "positive") positive.adinfo else negative.adinfo,
                                           polarity = ionization, normalizeScore = TRUE)

        peakTab <- as.data.table(cliqueMS::getPeaklistanClique(cliques))
        setnames(peakTab, "rt", "ret")
        peakTab[, ID := seq_len(.N)]
        peakTab <- peakTab[, c("ID", "ret", "mz", "cliqueGroup", "isotope", paste0("an", seq_len(5)),
                               paste0("score", seq_len(5)), paste0("mass", seq_len(5)))]
        setTxtProgressBar(prog, i)
        
        return(split(peakTab, by = "cliqueGroup", keep.by = FALSE))
    }), anas)
    
    close(prog)
    
    return(componentsCliqueMS(fGroups = fGroups, mzWindow = mzWindow, minSize = minSize,
                              relMinAdductAbundance = relMinAdductAbundance, featureComponents = featComponents))
})
