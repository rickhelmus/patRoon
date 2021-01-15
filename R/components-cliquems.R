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
    
    printf("Annotating all features with CliqueMS for %d analyses ...\n", length(anas))
    
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

        # For now we just take the highest ranking annotation. To further simplify, each clique is further separated per
        # neutral mass.
        
        peakTab <- as.data.table(cliqueMS::getPeaklistanClique(cliques))
        setnames(peakTab,
                 c("rt", "an1", "score1", "mass1"),
                 c("ret", "adduct", "score", "neutralMass"))
        peakTab[, ID := seq_len(.N)]
        
        # split isotope annotations in groups and clusters (similar to CAMERA/RAMClustR), remove defaulted annotations
        
        
        peakTab <- peakTab[, c("ID", "ret", "mz", "cliqueGroup", "isotope", "adduct", "score", "neutralMass"),
                           with = FALSE]
        
        # UNDONE: split works fine with (equal) numerics?
        ret <- split(peakTab, by = c("cliqueGroup", "neutralMass"), keep.by = FALSE)
        
        browser()
        
        setTxtProgressBar(prog, i)
        return(ret)
    }), anas)
    
    close(prog)
    
    return(componentsCliqueMS(fGroups = fGroups, mzWindow = mzWindow, minSize = minSize,
                              relMinAdductAbundance = relMinAdductAbundance, featureComponents = featComponents))
})
