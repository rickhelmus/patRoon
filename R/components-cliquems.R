# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
#' @include components-features.R
NULL

#' @rdname components-class
#' @export
componentsCliqueMS <- setClass("componentsCliqueMS", slots = c(cliques = "list"),
                               contains = "componentsFeatures")

setMethod("initialize", "componentsCliqueMS",
          function(.Object, ...) callNextMethod(.Object, algorithm = "cliquems", ...))

#' Componentization of adducts, isotopes etc. with cliqueMS
#'
#' Uses \href{https://github.com/osenan/cliqueMS}{cliqueMS} to generate components using the
#' \code{\link[cliqueMS:getCliques]{cliqueMS::getCliques}} function.
#'
#' @templateVar algo cliqueMS
#' @templateVar do generate components
#' @templateVar generic generateComponents
#' @templateVar algoParam cliquems
#' @template algo_generator
#'
#' @details The grouping of features in each component ('clique') is based on high similarity of chromatographic elution
#'   profiles. All features in each component are then annotated with the
#'   \code{\link[cliqueMS:getIsotopes]{cliqueMS::getIsotopes}} and
#'   \code{\link[cliqueMS:getAnnotation]{cliqueMS::getAnnotation}} functions.
#'
#' @param maxCharge,maxGrade,ppm Arguments passed to \code{\link[cliqueMS:getIsotopes]{cliqueMS::getIsotopes}} and/or
#'   \code{\link[cliqueMS:getAnnotation]{cliqueMS::getAnnotation}}.
#' @param adductInfo Sets the \code{adinfo} argument to \code{\link[cliqueMS:getAnnotation]{cliqueMS::getAnnotation}}.
#'   If \code{NULL} then the default adduct information from \pkg{cliqueMS} is used (\emph{i.e.} the
#'   \code{positive.adinfo}/\code{negative.adinfo} package datasets).
#' @param absMzDev Maximum absolute \emph{m/z} deviation.
#' @param extraOptsCli,extraOptsIso,extraOptsAnn Named \code{list} with further arguments to be passed to
#'   \code{\link[cliqueMS:getCliques]{cliqueMS::getCliques}}, \code{\link[cliqueMS:getIsotopes]{cliqueMS::getIsotopes}}
#'   and \code{\link[cliqueMS:getAnnotation]{cliqueMS::getAnnotation}}, respectively. Set to \code{NULL} to ignore.
#'
#' @templateVar ion TRUE
#' @templateVar minSize TRUE
#' @template compon_algo-args
#' @template parallel-arg
#'
#' @inheritParams generateComponents
#'
#' @return A \code{\link{componentsFeatures}} derived object.
#'
#' @template compon_gen-feat
#'
#' @templateVar class componentsSet
#' @template compon_gen-sets-merged
#'
#' @references \insertRef{Senan2019}{patRoon}
#'
#' @templateVar what generateComponentsCliqueMS
#' @template main-rd-method
#' @export
setMethod("generateComponentsCliqueMS", "featureGroups", function(fGroups, ionization = NULL, maxCharge = 1,
                                                                  maxGrade = 2, ppm = 10,
                                                                  adductInfo = NULL,
                                                                  absMzDev = 0.005, minSize = 2,
                                                                  relMinAdductAbundance = 0.75,
                                                                  adductConflictsUsePref = TRUE,
                                                                  NMConflicts = c("preferential", "mostAbundant", "mostIntense"),
                                                                  prefAdducts = c("[M+H]+", "[M-H]-"),
                                                                  extraOptsCli = NULL, extraOptsIso = NULL,
                                                                  extraOptsAnn = NULL, parallel = TRUE)
{
    checkPackage("cliqueMS", "rickhelmus/cliqueMS") # UNDONE
    
    ac <- checkmate::makeAssertCollection()
    ionization <- checkAndGetIonization(ionization, fGroups, add = ac)
    aapply(checkmate::assertCount, . ~ maxCharge + maxGrade + minSize, positive = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ ppm + absMzDev + relMinAdductAbundance, finite = TRUE, lower = 0,
           fixed = list(add = ac))
    checkmate::assertDataFrame(adductInfo, any.missing = FALSE, col.names = "unique", null.ok = TRUE, add = ac)
    checkmate::assertFlag(adductConflictsUsePref, add = ac)
    checkmate::assertSubset(NMConflicts, c("preferential", "mostAbundant", "mostIntense"), empty.ok = FALSE, add = ac)
    checkmate::assertCharacter(prefAdducts, min.chars = 1, any.missing = FALSE, unique = TRUE, add = ac)
    aapply(checkmate::assertList, . ~ extraOptsCli + extraOptsIso + extraOptsAnn, any.missing = FALSE,
           names = "unique", null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertFlag(parallel, add = ac)
    checkmate::reportAssertions(ac)

    if (length(fGroups) == 0)
        return(componentsCliqueMS(fGroups = fGroups, absMzDev = absMzDev, minSize = minSize,
                                  relMinAdductAbundance = relMinAdductAbundance))
        
    anas <- analyses(fGroups)
    fList <- getFeatures(fGroups)
    
    printf("Annotating all features with CliqueMS for %d analyses ...\n", length(anas))

    if (is.null(adductInfo))
    {
        if (ionization == "positive")
        {
            data(positive.adinfo, package = "cliqueMS", envir = environment())
            adductInfo <- positive.adinfo
        }
        else
        {
            data(negative.adinfo, package = "cliqueMS", envir = environment())
            adductInfo <- negative.adinfo
        }
    }

    db <- openCacheDBScope()
    baseHash <- makeHash(ionization, maxCharge, maxGrade, ppm, adductInfo, absMzDev, minSize, relMinAdductAbundance,
                         extraOptsCli, extraOptsIso, extraOptsAnn)
    hashes <- sapply(featureTable(fGroups), makeHash, baseHash)
    cachedCliques <- sapply(hashes, loadCacheData, category = "componentsCliqueMS", dbArg = db, simplify = FALSE)
    cachedCliques <- pruneList(cachedCliques)
    
    getCliques <- function(xdata)
    {
        suppressMessages(invisible(utils::capture.output({
            cliques <- do.call(cliqueMS::getCliques, c(list(xdata), extraOptsCli))
            cliques <- do.call(cliqueMS::getIsotopes,
                               c(list(cliques, maxCharge = maxCharge, maxGrade = maxGrade, ppm = ppm), extraOptsIso))
            cliques <- do.call(cliqueMS::getAnnotation,
                               c(list(cliques, ppm = ppm, adinfo = adductInfo, polarity = ionization,
                                      normalizeScore = TRUE), extraOptsAnn))
        })))
        patRoon:::doProgress()
        return(cliques)
    }

    anasTBD <- setdiff(anas, names(cachedCliques))
    if (length(anasTBD) > 0)
    {
        cat("Exporting to XCMS features... ")
        xds <- sapply(anasTBD, function(a) getXCMSnExp(fList[a], verbose = FALSE, loadRawData = TRUE),
                      simplify = FALSE)
        cat("Done!\n")
        
        allCliques <- doApply("sapply", parallel, xds, getCliques, simplify = FALSE)

        for (a in anasTBD)
            saveCacheData("componentsCliqueMS", allCliques[[a]], hashes[[a]], db)
        
        if (length(cachedCliques) > 0)
            allCliques <- c(allCliques, cachedCliques)[anas] # merge and re-order
    }
    else
        allCliques <- cachedCliques
        
    featComponents <- list()
    if (length(allCliques) > 0)
    {
        # UNDONE also cache?
        
        printf("Processing %d cliques... ",
               sum(sapply(allCliques, function(cl) length(cliqueMS::getlistofCliques(cl)))))
    
        featComponents <- sapply(allCliques, function(cl)
        {
            # For now we just take the highest ranking annotation. To further simplify, each clique is further separated per
            # neutral mass.
            peakTab <- as.data.table(cliqueMS::getPeaklistanClique(cl))
            setnames(peakTab,
                     c("rt", "an1", "score1", "mass1"),
                     c("ret", "adduct_ion", "score", "neutralMass"))
            peakTab[, ID := seq_len(.N)]
            peakTab[!nzchar(adduct_ion), adduct_ion := NA_character_] # BUG: sometimes adducts are empty instead of NA?
            peakTab[!is.na(adduct_ion),
                    adduct_ion := sapply(adduct_ion, function(a) as.character(as.adduct(a, format = "cliquems")))]
            
            # remove defaulted isotope annotations (eg w/out cluster)
            peakTab[!grepl("[", isotope, fixed = TRUE ), isotope := NA_character_]
            
            # split isotope annotations in groups and clusters (similar to CAMERA/RAMClustR)
            peakTab[!is.na(isotope), isogroup := sub("M[[:digit:]]+ \\[([[:digit:]]+)\\]", "\\1", isotope)]
            peakTab[!is.na(isotope), isonr := sub("M([[:digit:]]+).*", "\\1", isotope)]
            
            # get charges
            isoTab <- as.data.table(cliqueMS::getIsolistanClique(cl))
            # only keep isotope annotations with equal charge as M0, see https://github.com/osenan/cliqueMS/issues/7
            isoTab[, keep := charge == .SD[grade == 0]$charge, by = "cluster"]
            isoTab <- isoTab[keep == TRUE]
            # remove any clusters with size <2 that now may have appeared
            isoTab[, keep := .N >= 2, by = "cluster"]
            isoTab <- isoTab[keep == TRUE]
            peakTab[!is.na(isotope), charge := isoTab[match(ID, feature)]$charge]
            
            # unassign removed clusters
            peakTab[!is.na(isotope) & !isogroup %in% isoTab$cluster, c("isogroup", "isonr") := NA]
            
            peakTab <- peakTab[, c("ID", "ret", "mz", "cliqueGroup", "isogroup", "isonr", "charge", "adduct_ion",
                                   "score", "neutralMass"),
                               with = FALSE]
            
            # UNDONE: split works fine with (equal) numerics?
            ret <- unname(split(peakTab, by = c("cliqueGroup", "neutralMass"), keep.by = TRUE))
            
            return(ret)
        }, simplify = FALSE)
        
        printf("Done!\n")
    }
    
    return(componentsCliqueMS(cliques = allCliques, fGroups = fGroups, absMzDev = absMzDev, minSize = minSize,
                              relMinAdductAbundance = relMinAdductAbundance,
                              adductConflictsUsePref = adductConflictsUsePref, NMConflicts = NMConflicts,
                              prefAdducts = prefAdducts, featureComponents = featComponents))
})

#' @rdname generateComponentsCliqueMS
#' @export
setMethod("generateComponentsCliqueMS", "featureGroupsSet", function(fGroups, ionization = NULL, ...)
{
    generateComponentsSet(fGroups, ionization, generateComponentsCliqueMS, setIonization = TRUE, ...)
})

