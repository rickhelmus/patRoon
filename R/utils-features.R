makeFGroupName <- function(id, ret, mz) sprintf("M%.0f_R%.0f_%d", mz, ret, id)

showAnaInfo <- function(anaInfo)
{
    rGroups <- unique(anaInfo$group)
    blGroups <- unique(anaInfo$blank)
    printf("Analyses: %s (%d total)\n", getStrListWithMax(anaInfo$analysis, 6, ", "), nrow(anaInfo))
    printf("Replicate groups: %s (%d total)\n", getStrListWithMax(rGroups, 8, ", "), length(rGroups))
    printf("Replicate groups used as blank: %s (%d total)\n", getStrListWithMax(blGroups, 8, ", "), length(blGroups))
}

reGenerateFTIndex <- function(fGroups)
{
    gNames <- names(fGroups)
    fGroups@ftindex <- setnames(rbindlist(lapply(featureTable(fGroups),
                                                 function(ft) as.list(chmatch(gNames, ft$group, 0)))), gNames)
    return(fGroups)
}
isFGSet <- function(fGroups) inherits(fGroups, "featureGroupsSet")

featureQualities <- function()
{
    checkPackage("MetaClean")
    list(ApexBoundaryRatio = list(func = MetaClean::calculateApexMaxBoundaryRatio, HQ = "LV", range = c(0, 1)),
         FWHM2Base = list(func = MetaClean::calculateFWHM, HQ = "HV", range = c(0, 1)),
         Jaggedness = list(func = MetaClean::calculateJaggedness, HQ = "LV", range = Inf),
         Modality = list(func = MetaClean::calculateModality, HQ = "LV", range = Inf),
         Symmetry = list(func = MetaClean::calculateSymmetry, HQ = "HV", range = c(-1, 1)),
         GaussianSimilarity = list(func = MetaClean::calculateGaussianSimilarity, HQ = "HV", range = c(0, 1)),
         Sharpness = list(func = MetaClean::calculateSharpness, HQ = "HV", range = Inf),
         TPASR = list(func = MetaClean::calculateTPASR, HQ = "LV", range = Inf),
         ZigZag = list(func = MetaClean::calculateZigZagIndex, HQ = "LV", range = Inf))
}

featureGroupQualities <- function()
{
    checkPackage("MetaClean")
    list(
        ElutionShift = list(func = MetaClean::calculateElutionShift, HQ = "LV", range = Inf),
        RetentionTimeCorrelation = list(func = MetaClean::calculateRetentionTimeConsistency, HQ = "LV", range = Inf)
    )
}

# normalize, invert if necessary to get low (worst) to high (best) order
scoreFeatQuality <- function(quality, values)
{
    if (all(is.finite(quality$range)))
    {
        if (!isTRUE(all.equal(quality$range, c(0, 1)))) # no need to normalize 0-1
            values <- normalize(values, minMax = quality$range[1] < 0, xrange = quality$range)
    }
    else
        values <- normalize(values, TRUE)
    
    if (quality$HQ == "LV")
        values <- 1 - values
    
    return(values)
}

hasFGroupScores <- function(fGroups) nrow(groupScores(fGroups)) > 0

doFGroupsFilter <- function(fGroups, what, hashParam, func, cacheCateg = what, verbose = TRUE)
{
    if (verbose)
    {
        printf("Applying %s filter... ", what)
        oldFCount <- length(getFeatures(fGroups)); oldGCount <- length(fGroups)
    }
    
    cacheName <- sprintf("filterFGroups_%s", cacheCateg)
    hash <- makeHash(fGroups, what, hashParam)
    ret <- loadCacheData(cacheName, hash)
    if (is.null(ret))
    {
        ret <- if (length(fGroups) > 0) func(fGroups) else fGroups
        saveCacheData(cacheName, ret, hash)
    }
    
    if (verbose)
    {
        newFCount <- length(getFeatures(ret)); newGCount <- length(ret)
        newn <- length(ret)
        printf("Done! Filtered %d (%.2f%%) features and %d (%.2f%%) feature groups. Remaining: %d features in %d groups.\n",
               oldFCount - newFCount, if (oldFCount > 0) (1 - (newFCount / oldFCount)) * 100,
               oldGCount - newGCount, if (oldGCount > 0) (1 - (newGCount / oldGCount)) * 100,
               newFCount, newGCount)
    }
    
    return(ret)
}

# used by adducts()<- methods
updateAnnAdducts <- function(annTable, gInfo, adducts)
{
    if (nrow(annTable) > 0)
    {
        # only consider changed
        adducts <- adducts[adducts != annTable$adduct]
        
        if (length(adducts) == 0)
            return(annTable) # nothing changed
    }
    
    adducts <- sapply(adducts, checkAndToAdduct, .var.name = "value", simplify = FALSE)
    adductsChr <- sapply(adducts, as.character) # re-make characters: standardize format
    nm <- calculateMasses(gInfo[names(adducts), "mzs"], adducts, type = "neutral")
    
    if (nrow(annTable) > 0)
    {
        # update table
        annTable <- copy(annTable)
        annTable[match(names(adducts), group), c("adduct", "neutralMass") := .(adductsChr, nm)][]
    }
    else # initialize table
        annTable <- data.table(group = names(adducts), adduct = adductsChr, neutralMass = nm)
    
    return(annTable)
}

unSetAnaInfo <- function(anaInfo) anaInfo[, setdiff(names(anaInfo), "set")]

maybeAutoNormalizeFGroups <- function(fGroups)
{
    if (length(fGroups) == 0 || !is.null(featureTable(fGroups)[[1]][["intensity_rel"]]))
        return(fGroups) # no features or already normalized
    printf("Automatically normalizing feature groups, see ?normInts for more options.\n")
    return(normInts(fGroups, featNorm = "none", groupNorm = TRUE))
}

getAnnotationsFromSetFeatures <- function(fGroups)
{
    if (length(fGroups) > 0)
    {
        anaInfo <- analysisInfo(fGroups)
        ftind <- groupFeatIndex(fGroups)
        fTable <- featureTable(fGroups)
        
        ret <- rbindlist(sapply(sets(fGroups), function(s)
        {
            anaInds <- which(anaInfo$set == s)
            anas <- anaInfo[anaInds, "analysis"]
            grps <- names(fGroups)[sapply(ftind[anaInds], function(x) any(x != 0))]
            firstFeats <- rbindlist(lapply(ftind[anaInds, grps, with = FALSE], function(x)
            {
                firstAna <- which(x != 0)[1]
                return(fTable[[anas[firstAna]]][x[firstAna]])
            }))
            
            return(data.table(group = grps, adduct = firstFeats$adduct))
        }, simplify = FALSE), idcol = "set", fill = TRUE) # set fill for empty objects
        ret[, neutralMass := groupInfo(fGroups)[ret$group, "mzs"]]
    }
    else
        ret <- data.table()
}

filterEICs <- function(EICs, fGroups, analysis = NULL, groupName = NULL, topMost = NULL, topMostByRGroup = FALSE,
                       onlyPresent = FALSE)
{
    if (!is.null(analysis))
        EICs <- EICs[names(EICs) %chin% analysis]
    else
        analysis <- analyses(fGroups)
    if (!is.null(groupName))
        EICs <- lapply(EICs, function(e) e[names(e) %chin% groupName])

    if (onlyPresent)
    {
        gTable <- groupTable(fGroups)
        EICs <- Map(names(EICs), EICs, f = function(ana, aeic)
        {
            anaInd <- match(ana, analyses(fGroups))
            absentFGs <- names(gTable)[gTable[anaInd] == 0]
            return(aeic[!names(aeic) %chin% absentFGs])
        })
    }
    
    if (!is.null(topMost))
    {
        gTable <- copy(groupTable(fGroups))
        gTable[, c("analysis", "rGroup") := analysisInfo(fGroups)[, c("analysis", "group")]]
        for (fg in names(fGroups))
        {
            anasWithFG <- Map(names(EICs), EICs, f = function(ana, aeic) if (fg %chin% names(aeic)) ana else character())
            anasWithFG <- pruneList(anasWithFG, checkEmptyElements = TRUE)
            anasWithFG <- unlist(anasWithFG)
            tab <- gTable[analysis %chin% anasWithFG, c(fg, "analysis", "rGroup"), with = FALSE]
            if (nrow(tab) > 0)
            {
                rmAnas <- if (topMostByRGroup)
                {
                    tab[, rank := frank(-get(fg), ties.method = "first"), by = "rGroup"]
                    tab[rank > topMost]$analysis
                }
                else if (nrow(tab) > topMost)
                {
                    setorderv(tab, fg, order = -1L)
                    tab[seq(topMost + 1, nrow(tab))]$analysis
                }
                else
                    character()
                EICs[rmAnas] <- lapply(EICs[rmAnas], "[[<-", fg, NULL)
            }
        }
    }

    return(pruneList(EICs, checkEmptyElements = TRUE))
}

massConcUnitBase <- function(mcu)
{
    return(switch(mcu,
                  gL = 1,
                  mgL = 1E-3,
                  ugL = 1E-6,
                  ngL = 1E-9))
}

convertConc <- function(conc, unitFrom, unitTo, MW)
{
    parseUnit <- function(unit)
    {
        logBase <- if (startsWith(unit, "log10 "))
            10
        else if (startsWith(unit, "log2 "))
            2
        else if (startsWith(unit, "log "))
            exp(1)
        else
            0
        if (logBase != 0)
            unit <- sub("^log(2|10)? ", "", unit)
        base <- switch(substring(unit, 1, 1),
                       n = 1E-9,
                       u = 1E-6,
                       m = 1E-3,
                       1)
        molar <- endsWith(unit, "M")
        return(list(logBase = logBase, base = base, molar = molar))
    }
    
    unitFromP <- parseUnit(unitFrom); unitToP <- parseUnit(unitTo)
    
    if (unitFromP$logBase != 0 && unitToP$logBase == 0)
        conc <- unitFromP$logBase^conc
    else if (unitFromP$logBase == 0 && unitToP$logBase)
        conc <- log(conc, unitToP$logBase)
    else if (unitFromP$logBase != unitToP$logBase)
        conc <- log(unitFromP$logBase^conc, unitToP$logBase)
    
    conc <- conc * (unitFromP$base / unitToP$base)
    
    if (unitFromP$molar && !unitToP$molar)
        conc <- conc * MW
    else if (!unitFromP$molar && unitToP$molar)
        conc <- conc / MW
    
    return(conc)
}

calcFeatureConcs <- function(fGroups, resp, areas)
{
    # get concentration data from response factors
    
    gt <- transpose(groupTable(fGroups, areas = areas)[, resp$group, with = FALSE])
    setnames(gt, analyses(fGroups))
    
    concs <- copy(resp)
    concs[, (analyses(fGroups)) := lapply(gt, function(ints) ints / RF)]
    concs[, RF := NULL]

    return(concs[])
}

finalizeFeaturePredTab <- function(pred)
{
    co <- intersect(c("group", "type", "candidate", "candidate_name", "sets"), names(pred))
    setcolorder(pred, co)
    setorderv(pred, c("group", "type", "candidate"))
    return(pred)
}

doCalcConcSets <- function(fGroups, featureAnn, areas)
{
    if (is.null(featureAnn))
    {
        # just from screening info
        usFGroups <- sapply(sets(fGroups), unset, obj = fGroups, simplify = FALSE)
        usFGroups <- sapply(usFGroups, calculateConcs, areas = areas, simplify = FALSE)
    }
    else
    {
        usFeatAnns <- checkAndUnSetOther(sets(fGroups), featureAnn, "featureAnn")
        usFeatAnns <- usFeatAnns[lengths(usFeatAnns) > 0]
        
        if (length(usFeatAnns) == 0)
        {
            cat("No feature annotations, nothing to do...\n")
            return(fGroups)
        }
        
        usFGroups <- sapply(names(usFeatAnns), unset, obj = fGroups, simplify = FALSE)
        usFGroups <- Map(usFGroups, usFeatAnns, f = calculateConcs, MoreArgs = list(areas = areas))
    }
    
    usFGroups <- usFGroups[sapply(usFGroups, function(ufg) nrow(concentrations(ufg)) > 0)]
    if (length(usFGroups) == 0)
        return(fGroups)
    
    if (length(usFGroups) == 1)
        fGroups@concentrations <- copy(concentrations(usFGroups[[1]]))
    else
    {
        fGroups@concentrations <- Reduce(lapply(usFGroups, concentrations), f = function(left, right)
        {
            cols <- intersect(c("group", "type", "candidate", "candidate_name"), c(names(left), names(right)))
            return(merge(left, right, by = cols, all = TRUE))
        })
        fGroups@concentrations <- finalizeFeaturePredTab(fGroups@concentrations)
    }
    
    # add missing analyses
    missingAnaCols <- setdiff(analyses(fGroups), names(fGroups@concentrations))
    if (length(missingAnaCols) > 0)
        fGroups@concentrations[, (missingAnaCols) := NA_real_][]
    
    return(fGroups)
}

doCalcToxSets <- function(fGroups, featureAnn)
{
    if (is.null(featureAnn))
    {
        # just from screening info
        usFGroups <- sapply(sets(fGroups), unset, obj = fGroups, simplify = FALSE)
        usFGroups <- sapply(usFGroups, calculateTox, simplify = FALSE)
    }
    else
    {
        usFeatAnns <- checkAndUnSetOther(sets(fGroups), featureAnn, "featureAnn")
        usFeatAnns <- usFeatAnns[lengths(usFeatAnns) > 0]
        
        if (length(usFeatAnns) == 0)
        {
            cat("No feature annotations, nothing to do...\n")
            return(fGroups)
        }
        
        usFGroups <- sapply(names(usFeatAnns), unset, obj = fGroups, simplify = FALSE)
        usFGroups <- Map(usFGroups, usFeatAnns, f = calculateTox)
    }
    
    usFGroups <- usFGroups[sapply(usFGroups, function(ufg) nrow(toxicities(ufg)) > 0)]
    if (length(usFGroups) == 0)
        return(fGroups)
    
    if (length(usFGroups) == 1)
    {
        fGroups@toxicities <- copy(toxicities(usFGroups[[1]]))
        fGroups@toxicities[, sets := names(usFGroups)]
    }
    else
    {
        allTox <- sapply(usFGroups, toxicities, simplify = FALSE)
        
        allToxTab <- rbindlist(allTox, idcol = "sets")
        
        # NOTE: only for SIRIUS_FP it makes sense to split sets data, the rest we can merge
        allToxTab[type != "SIRIUS_FP", sets := paste0(unique(sets), collapse = ","), by = c("group", "candidate", "type")]
        allToxTab <- unique(allToxTab, by = c("group", "type", "candidate", "sets"))
        
        fGroups@toxicities <- finalizeFeaturePredTab(allToxTab)
    }
    
    return(fGroups)
}
