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
                                                                  relMinAdductAbundance = 0.75)
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
    checkmate::reportAssertions(ac)

    if (length(fGroups) == 0)
        return(componentsCliqueMS())
        
    anas <- analyses(fGroups)
    fList <- getFeatures(fGroups)
    
    printf("Annotating all features with CliqueMS for %d analyses ...\n", length(anas))
    
    if (ionization == "positive")
        data(positive.adinfo, package = "cliqueMS", envir = environment())
    else
        data(negative.adinfo, package = "cliqueMS", envir = environment())

    db <- openCacheDBScope()
    baseHash <- makeHash(ionization, mzWindow, minSize, relMinAdductAbundance)
    
    prog <- openProgBar(0, length(anas))
    
    featComponents <- setNames(lapply(seq_along(anas), function(i)
    {
        hash <- makeHash(fList[[i]], baseHash)
        ret <- loadCacheData("componentsCliqueMS", hash, db)
        if (is.null(ret))
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
            peakTab[!is.na(adduct) & nzchar(adduct),
                    adduct := sapply(adduct, function(a) as.character(as.adduct(a, format = "cliquems")))]
            
            # remove defaulted isotope annotations (eg w/out cluster)
            peakTab[!grepl("[", isotope, fixed = TRUE ), isotope := NA_character_]
            
            # split isotope annotations in groups and clusters (similar to CAMERA/RAMClustR)
            peakTab[!is.na(isotope), isogroup := sub("M[[:digit:]]+ \\[([[:digit:]]+)\\]", "\\1", isotope)]
            peakTab[!is.na(isotope), isonr := sub("M([[:digit:]]+).*", "\\1", isotope)]
            
            # get charges
            isoTab <- as.data.table(cliqueMS::getIsolistanClique(cliques))
            # only keep isotope annotations with equal charge as M0, see https://github.com/osenan/cliqueMS/issues/7
            isoTab[, keep := charge == .SD[grade == 0]$charge, by = "cluster"]
            isoTab <- isoTab[keep == TRUE]
            # remove any clusters with size <2 that now may have appeared
            isoTab[, keep := .N >= 2, by = "cluster"]
            isoTab <- isoTab[keep == TRUE]
            peakTab[!is.na(isotope), charge := isoTab[match(ID, feature)]$charge]
            
            peakTab <- peakTab[, c("ID", "ret", "mz", "cliqueGroup", "isogroup", "isonr", "charge", "adduct",
                                   "score", "neutralMass"),
                               with = FALSE]
            
            # UNDONE: split works fine with (equal) numerics?
            ret <- split(peakTab, by = c("cliqueGroup", "neutralMass"), keep.by = TRUE)
        }
        
        saveCacheData("componentsCliqueMS", ret, hash, db)
        
        setTxtProgressBar(prog, i)
        return(ret)
    }), anas)
    
    close(prog)
    
    return(componentsCliqueMS(fGroups = fGroups, mzWindow = mzWindow, minSize = minSize,
                              relMinAdductAbundance = relMinAdductAbundance, featureComponents = featComponents))
})
