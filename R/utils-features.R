makeFGroupName <- function(id, ret, mz) sprintf("M%.0f_R%.0f_%d", mz, ret, id)

showAnaInfo <- function(anaInfo)
{
    rGroups <- unique(anaInfo$group)
    blGroups <- unique(anaInfo$blank)
    printf("Analyses: %s (%d total)\n", getStrListWithMax(anaInfo$analysis, 6, ", "), nrow(anaInfo))
    printf("Replicate groups: %s (%d total)\n", getStrListWithMax(rGroups, 8, ", "), length(rGroups))
    printf("Replicate groups used as blank: %s (%d total)\n", getStrListWithMax(blGroups, 8, ", "), length(blGroups))
}

checkAnaInfoAggrGrouping <- function(anaInfo, what, aggrBy, groupBy)
{
    for (ag in unique(anaInfo[[aggrBy]]))
    {
        ai <- anaInfo[get(aggrBy) == ag]
        if (uniqueN(ai[[groupBy]]) > 1)
        {
            stop(sprintf("The following analyses will be %s but do not have the same value set for the '%s' column in the analysis information:\n%s",
                         what, groupBy, paste0(sprintf("%s: %s", ai$analysis, ai[[groupBy]]), collapse = "\n")), call. = FALSE)
        }
    }
}

doSubsetFeaturesByAna <- function(obj, i, ni, reorder)
{
    if (!missing(ni))
    {
        if (!missing(i))
            stop("Cannot simulatenously specify i and ni arguments.", call. = FALSE)
        anaInfo <- analysisInfo(obj)
        # From https://stackoverflow.com/a/62043225
        # NOTE: set env=parent.frame() here, as ni is passed from another function (`[`) 
        anaInfo <- eval(substitute(anaInfo[.i], list(.i = substitute(ni, env = parent.frame()))))
        i <- anaInfo$analysis
    }
    
    if (!missing(i))
    {
        i <- assertSubsetArgAndToChr(i, analyses(obj))
        obj <- delete(obj, setdiff(analyses(obj), i))
        if (reorder && !isTRUE(all.equal(i, analyses(obj), check.attributes = FALSE)))
            obj <- reorderAnalyses(obj, i)
    }
    
    return(obj)
}

reGenerateFTIndex <- function(fGroups)
{
    gNames <- names(fGroups)
    fGroups@ftindex <- setnames(rbindlist(lapply(featureTable(fGroups),
                                                 function(ft) as.list(chmatch(gNames, ft$group, 0)))), gNames)
    return(fGroups)
}
isFGSet <- function(fGroups) inherits(fGroups, "featureGroupsSet")

calcFeatureRegression <- function(xvec, ints)
{
    NARet <- list(RSQ = NA_real_, intercept = NA_real_, slope = NA_real_, p = NA_real_, lm = NULL)
    
    notna <- !is.na(xvec)
    if (!any(notna) || length(xvec) == 1)
        return(NARet)
    
    ints <- ints[notna]
    ints[ints == 0] <- NA
    if (all(is.na(ints)))
        return(NARet)
    
    rlm <- lm(ints ~ xvec[notna])
    suppressWarnings(reg <- summary(rlm))
    slope <- pv <- NA_real_
    if (nrow(reg[["coefficients"]]) > 1)
    {
        slope <- reg[["coefficients"]][2, 1]
        pv <- reg[["coefficients"]][2, 4]
    }
    
    return(list(RSQ = reg[["r.squared"]], intercept = reg[["coefficients"]][1, 1], slope = slope, p = pv, lm = rlm))
}

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
            anas <- anaInfo$analysis[anaInds]
            grps <- names(fGroups)[sapply(ftind[anaInds], function(x) any(x != 0))]
            firstFeats <- rbindlist(lapply(ftind[anaInds, grps, with = FALSE], function(x)
            {
                firstAna <- which(x != 0)[1]
                return(fTable[[anas[firstAna]]][x[firstAna]])
            }))
            
            return(data.table(group = grps, adduct = firstFeats$adduct))
        }, simplify = FALSE), idcol = "set", fill = TRUE) # set fill for empty objects
        ret[, neutralMass := groupInfo(fGroups)[ret$group, "mzs"]]
        adducts <- sapply(unique(ret$adduct), as.adduct)
        ret[, ion_mz := calculateMasses(neutralMass, adducts[adduct], type = "mz")]
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
        gTable[, c("analysis", "rGroup") := analysisInfo(fGroups)[, c("analysis", "group"), with = FALSE]]
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

aggregateConcs <- function(concs, anaInfo, aggrParams, splitSuspects = FALSE)
{
    concs <- copy(concs)
    concs <- subsetDTColumnsIfPresent(concs, c("group", "type", "candidate", "candidate_name", anaInfo$analysis))
    
    if (splitSuspects && !any(concs$type == "suspect"))
        splitSuspects <- FALSE # no suspects, nothing to split
    
    if (aggrParams$preferType != "none")
    {
        concs[, keep := !aggrParams$preferType %in% type | type == aggrParams$preferType, by = "group"]
        concs <- concs[keep == TRUE][, keep := NULL]
    }
    
    ignoreFGs <- character()
    if (splitSuspects)
    {
        # fgs with suspect results will be reported for each suspect per row and results from other types will be
        # ignored
        fgsWithSuspect <- unique(concs[type == "suspect"]$group)
        concs <- concs[!group %chin% fgsWithSuspect | type == "suspect"] # remove all others
        ignoreFGs <- fgsWithSuspect # no further aggregation needed
    }
    else if (!is.null(concs[["candidate_name"]]))
        concs[, candidate_name := NULL]

    doAggr <- function(func, by, combineTypes = FALSE)
    {
        concs[!group %chin% ignoreFGs, (anaInfo$analysis) := lapply(.SD, aggrVec, func),
              .SDcols = anaInfo$analysis, by = by]
        dups <- duplicated(concs, by = by)
        if (combineTypes)
            concs[, type := paste0(unique(type), collapse = ","), by = by]            
        return(concs[group %chin% ignoreFGs | dups == FALSE])
    }
    concs <- doAggr(aggrParams$candidateFunc, c("group", "type", "candidate"))
    concs <- doAggr(aggrParams$typeFunc, c("group", "type"))
    concs <- doAggr(aggrParams$groupFunc, "group", combineTypes = TRUE)
    
    mby <- if (splitSuspects) c("group", "candidate_name") else "group"
    concs <- unique(concs, by = mby)
    
    concs[, candidate := NULL]
    
    return(concs[])
}

aggregateTox <- function(tox, aggrParams, splitSuspects = FALSE)
{
    tox <- copy(tox)
    tox <- subsetDTColumnsIfPresent(tox, c("group", "type", "candidate", "candidate_name", "LC50"))
    
    if (splitSuspects && !any(tox$type == "suspect"))
        splitSuspects <- FALSE # no suspects, nothing to split

    if (aggrParams$preferType != "none")
    {
        tox[, keep := !aggrParams$preferType %in% type | type == aggrParams$preferType, by = "group"]
        tox <- tox[keep == TRUE][, keep := NULL]
    }
    
    ignoreFGs <- character()
    if (splitSuspects)
    {
        # fgs with suspect results will be reported for each suspect per row and results from other types will be
        # ignored
        fgsWithSuspect <- unique(tox[type == "suspect"]$group)
        tox <- tox[!group %chin% fgsWithSuspect | type == "suspect"] # remove all others
        ignoreFGs <- fgsWithSuspect # no further aggregation needed
    }
    else if (!is.null(tox[["candidate_name"]]))
        tox[, candidate_name := NULL]

    doAggr <- function(func, by, combineTypes = FALSE)
    {
        tox[!group %chin% ignoreFGs, LC50 := aggrVec(LC50, func), by = by]
        dups <- duplicated(tox, by = by)
        if (combineTypes)
            tox[, type := paste0(unique(type), collapse = ","), by = by]
        return(tox[group %chin% ignoreFGs | dups == FALSE])
    }
    tox <- doAggr(aggrParams$candidateFunc, c("group", "type", "candidate"))
    tox <- doAggr(aggrParams$typeFunc, c("group", "type"))
    tox <- doAggr(aggrParams$groupFunc, "group", combineTypes = TRUE)

    mby <- if (splitSuspects) c("group", "candidate_name") else "group"
    tox <- unique(tox, by = mby)
    
    tox[, candidate := NULL]
    
    return(tox[])
}
