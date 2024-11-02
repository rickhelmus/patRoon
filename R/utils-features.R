# SPDX-FileCopyrightText: 2016 - 2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

makeFGroupName <- function(id, ret, mz) sprintf("M%.0f_R%.0f_%d", mz, ret, id)

appendMobToName <- function(n, mob) make.unique(sprintf("%s_I%.2f", n, mob)) # UNDONE: does this scale well for non-Bruker data?

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

getFilteredFGroups <- function(fGroups, retFilter) 
{
    if (!is.null(retFilter$subset))
    {
        if (is.null(retFilter$subset$i) && !is.null(retFilter$subset$j))
            fGroups <- fGroups[, retFilter$subset$j]
        else if (!is.null(retFilter$subset$i) && is.null(retFilter$subset$j))
            fGroups <- fGroups[retFilter$subset$i]
        else
            fGroups <- fGroups[retFilter$subset$i, retFilter$subset$j]  
    }
    if (!is.null(retFilter$delete))
        fGroups <- delete(fGroups, i = retFilter$delete$i, j = retFilter$delete$j)
    return(fGroups)
}

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
        ret <- if (length(fGroups) > 0) func(fGroups) else NULL
        if (!is.null(ret))
            saveCacheData(cacheName, ret, hash)
    }

    fGroups <- getFilteredFGroups(fGroups, ret)
    
    if (verbose)
    {
        newFCount <- length(getFeatures(fGroups)); newGCount <- length(fGroups)
        newn <- length(ret)
        printf("Done! Filtered %d (%.2f%%) features and %d (%.2f%%) feature groups. Remaining: %d features in %d groups.\n",
               oldFCount - newFCount, if (oldFCount > 0) (1 - (newFCount / oldFCount)) * 100,
               oldGCount - newGCount, if (oldGCount > 0) (1 - (newGCount / oldGCount)) * 100,
               newFCount, newGCount)
    }
    
    return(fGroups)
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
    nm <- calculateMasses(gInfo[match(names(adducts), group)]$mz, adducts, type = "neutral")
    
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
        ret[, neutralMass := groupInfo(fGroups)$mz[match(ret$group, groupInfo(fGroups)$group)]]
        adducts <- sapply(unique(ret$adduct), as.adduct)
        ret[, ion_mz := calculateMasses(neutralMass, adducts[adduct], type = "mz")]
    }
    else
        ret <- data.table()
}

finishFGroupsForSets <- function(fGroups, ..., verbose)
{
    otherArgs <- list(...)
    feat <- getFeatures(fGroups) # may have been changed (eg in initialize())
    ret <- featureGroupsSet(groupAlgo = algorithm(fGroups), groupArgs = otherArgs, groupVerbose = verbose,
                            groups = groupTable(fGroups), groupInfo = groupInfo(fGroups), features = feat,
                            ftindex = groupFeatIndex(fGroups), algorithm = makeSetAlgorithm(list(fGroups)))
    ret@annotations <- getAnnotationsFromSetFeatures(ret)
    return(ret)
}

groupFeaturesSets <- function(feat, grouper, ..., verbose)
{
    fGroups <- selectMethod(grouper, "features")(feat = feat, ...)
    return(finishFGroupsForSets(fGroups, ..., verbose = verbose))
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

setMethod("getEICFGroupInfo", "featureGroups", function(fGroups, analysis, groupName, EICParams)
{
    takeAnalysis <- analysis # copy name to workaround DT access below
    
    anaInfo <- analysisInfo(fGroups)[analysis %chin% takeAnalysis]
    featTab <- as.data.table(getFeatures(fGroups))
    
    topMost <- if (!is.null(EICParams$topMost))
        min(EICParams$topMost, nrow(anaInfo))
    else
        NULL
    
    # subset relevant things in advance
    featTab <- subsetDTColumnsIfPresent(featTab, c("group", "analysis", "intensity", "retmin", "retmax", "mzmin",
                                                   "mzmax", "mobmin", "mobmax"))
    
    # NOTE: we subset and split here in advance, as doing it in the loop below gets quite slow with many fGroups
    featTabAnaSub <- featTab[analysis %chin% takeAnalysis]
    featTabSplitGrp <- split(featTab, by = "group", keep.by = FALSE)
    featTabAnaSubSplitGrp <- split(featTabAnaSub, by = "group", keep.by = FALSE)
    
    return(sapply(groupName, function(fg)
    {
        ret <- featTabAnaSubSplitGrp[[fg]]
        if (is.null(ret))
            ret <- featTabSplitGrp[[fg]][0] # not present for this analysis, take full table to get all columns
        ret <- copy(ret)
        
        # add missing analyses if needed
        if (!EICParams$onlyPresent)
        {
            if (any(!analysis %chin% ret$analysis))
            {
                ftAllAna <- featTabSplitGrp[[fg]]
                ret <- rbind(ret, data.table(analysis = setdiff(analysis, ret$analysis), intensity = 0,
                                             retmin = min(ftAllAna$retmin), retmax = max(ftAllAna$retmax),
                                             mzmin = min(ftAllAna$mzmin) - EICParams$mzExpWindow,
                                             mzmax = max(ftAllAna$mzmax) + EICParams$mzExpWindow))
            }
        }
        
        if (!is.null(topMost))
        {
            if (EICParams$topMostByRGroup)
            {
                ret[, rGroup := anaInfo$group[match(analysis, anaInfo$analysis)]]
                ret[, rank := frank(-intensity, ties.method = "first"), by = "rGroup"]
                ret <- ret[rank <= topMost]
            }
            else
            {
                setorderv(ret, "intensity", order = -1L)
                ret <- ret[seq_len(topMost)]
            }
        }
        
        ret[, c("retmin", "retmax") := .(retmin - EICParams$rtWindow, retmax + EICParams$rtWindow)]
        return(ret)
    }, simplify = FALSE))
})

setMethod("getEICFGroupInfo", "featureGroupsSet", function(fGroups, analysis, groupName, EICParams)
{
    ret <- callNextMethod()
    
    anaInfo <- analysisInfo(fGroups)
    featTab <- as.data.table(getFeatures(fGroups))
    
    # HACK: since feature tables store the character form, it's easier to keep it all the same
    EICParams$setsAdductPos <- as.character(EICParams$setsAdductPos)
    EICParams$setsAdductNeg <- as.character(EICParams$setsAdductNeg)
    
    # 'ionize' m/zs
    return(Map(names(ret), ret, f = function(grp, ranges)
    {
        featTabGrp <- featTab[group == grp]
        ranges[, adduct := featTabGrp[match(ranges$analysis, analysis)]$adduct]
        
        if (!EICParams$onlyPresent && any(is.na(ranges$adduct))) # adduct will be NA for 'missing' features
        {
            # First try to get adduct from other features in the same set: assume that adduct per set for a single
            # feature group is always the same
            adductSets <- unique(featTabGrp[, c("adduct", "set"), with = FALSE])
            ranges[is.na(adduct), adduct := {
                s <- anaInfo$set[match(analysis, anaInfo$analysis)]
                adductSets[set == s]$adduct
            }, by = "analysis"]
            
            # Then fallback to default adducts for a polarity. For this we get the polarity from another feature in the
            # same analysis (if present)
            ranges[is.na(adduct), adduct := sapply(analysis, function(ana)
            {
                t <- featTab[analysis == ana]
                if (nrow(t) == 0)
                    NA # all features were removed
                else if (as.adduct(t$adduct[1])@charge > 0)
                    EICParams$setsAdductPos
                else
                    EICParams$setsAdductNeg
            })]
            
            if (any(is.na(ranges$adduct)))
            {
                # If all failed then simply omit
                warning(sprintf("Cannot get adduct information for group '%s' for features in analyses %s", grp,
                                paste0(ranges[is.na(adduct)]$analysis, collapse = ", ")), call. = FALSE)
                ranges <- ranges[!is.na(adduct)]
            }
        }
        
        # NOTE: this is mostly copied from unset features method        
        allAdducts <- sapply(unique(ranges$adduct), as.adduct)
        mzmins <- calculateMasses(ranges$mzmin, allAdducts[ranges$adduct], type = "mz")
        nmd <- mzmins - ranges$mzmin
        set(ranges, j = c("mzmin", "mzmax"), value = list(mzmins, ranges$mzmax + nmd))
        
        return(ranges)
    }))
})

setMethod("getEICsForFGroups", "featureGroups", function(fGroups, analysis, groupName, EICParams)
{
    if (length(fGroups) == 0 || length(analysis) == 0 || length(groupName) == 0)
        return(list())
    
    takeAnalysis <- analysis # for DT subset below
    anaInfo <- analysisInfo(fGroups)[analysis %chin% takeAnalysis]
    
    EICInfoTab <- getEICFGroupInfo(fGroups, analysis, groupName, EICParams)
    EICInfo <- split(rbindlist(EICInfoTab, idcol = "group"), by = "analysis")
    EICInfo <- EICInfo[intersect(anaInfo$analysis, names(EICInfo))] # sync order
    anaInfoEICs <- anaInfo[analysis %in% names(EICInfo)]
    
    EICs <- doGetEICs(anaInfoEICs, EICInfo)
    EICs <- Map(EICs, lapply(EICInfo, "[[", "group"), f = setNames)
    
    return(pruneList(EICs))
})

setMethod("getEICsForFeatures", "features", function(features, RTWindow = NULL, onlyMob = FALSE, ...)
{
    if (length(features) == 0)
        return(list())
    
    EICInfoList <- featureTable(features)
    if (!is.null(RTWindow) || onlyMob)
    {
        EICInfoList <- lapply(featureTable(features), function(ft)
        {
            if (!is.null(RTWindow))
            {
                ft <- copy(ft)
                ft[, c("retmin", "retmax") := .(retmin - RTWindow, retmax + RTWindow)]
            }
            if (onlyMob)
                ft <- ft[!is.null(ft[["mobility"]]) & !is.na(mobility)]
            return(ft)
        })
    }
    
    ret <- doGetEICs(analysisInfo(features), EICInfoList, ...)
    ret <- Map(ret, featureTable(features), f = function(eics, ft)
    {
        wh <- !onlyMob | (!is.null(ft[["mobility"]]) & !is.na(ft$mobility))
        names(eics) <- ft[wh]$ID
        return(eics)
    })
    return(pruneList(ret))
})

setMethod("getEICsForFeatures", "featuresSet", function(features, ...)
{
    unsetFeatList <- sapply(sets(features), unset, obj = features, simplify = FALSE)
    EICList <- sapply(unsetFeatList, getEICsForFeatures, simplify = FALSE)
    EICs <- unlist(EICList, recursive = FALSE, use.names = FALSE) # use.names gives combined set/ana name, we just want ana
    names(EICs) <- unlist(lapply(EICList, names))
    EICs <- EICs[intersect(analyses(features), names(EICs))] # sync order
    return(EICs)
})

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

findPeaksInEICs <- function(allEICs, peaksParam, withBP, parallel, cacheDB = NULL)
{
    baseHash <- makeHash(peaksParam)
    
    doApply("sapply", parallel, allEICs, function(EICs)
    {
        # NOTE: EICs must be named
        
        hash <- makeHash(baseHash, EICs)
        cd <- loadCacheData("peaksEIC", hash, cacheDB)
        if (!is.null(cd))
        {
            doProgress()
            return(cd)
        }

        # convert EICs to data.tables: this is necessary for findPeaks()
        # NOTE: we don't store (or hash) the EICs as DTs, as this makes things slower
        
        peaks <- findPeaks(EICs, peaksParam, verbose = FALSE)
        peaks <- rbindlist(peaks, idcol = "EIC_ID")
        
        if (nrow(peaks) == 0)
        {
            peaks[, EIC_ID := character()]
            peaks[, c("retmin", "retmax", "ret", "area", "intensity", "mzmin", "mzmax", "mz", "mobmin", "mobmax", "mobility") := numeric()]
        }
        else
        {
            peaks[, c("mzmin", "mzmax", "mz", "mobmin", "mobmax", "mobility") := {
                eic <- EICs[[EIC_ID]][EICs[[EIC_ID]]$intensity != 0 & EICs[[EIC_ID]]$time %between% c(retmin, retmax), ]
                if (nrow(eic) == 0)
                    numeric(1)
                else
                {
                    if (is.null(eic[["mobility"]]))
                        eic$mobility <- NA_real_
                    # UNDONE: also use mobility BP data?
                    list(min(eic$mzmin), max(eic$mzmax), weighted.mean(if (withBP) eic$mzBP else eic$mz, eic$intensity),
                         min(eic$mobility), max(eic$mobility), weighted.mean(eic$mobility, eic$intensity))
                }
            }, by = seq_len(nrow(peaks))]
        }
        
        # NOTE: we could also set mobilities after checking if data is available, but then we need to repeat the EIC subsetting above
        if (length(EICs) == 0 || is.null(EICs[[1]][["mobility"]]))
            peaks[, c("mobmin", "mobmax", "mobility") := NULL]
        
        # make unique IDs
        peaks[, ID := make.unique(EIC_ID)]
        
        saveCacheData("peaksEIC", peaks, hash, cacheDB)
        
        doProgress()
        
        return(peaks)
    }, simplify = FALSE)
}

assignFeatureMobilities <- function(features, peaksParam, mzWindow, IMSWindow, clusterMethod, minIntensityIMS,
                                    maxMSRTWindow, assignedMobilities)
{
    printf("Finding mobilities for all features...\n")
    
    anaInfo <- analysisInfo(features)
    
    features@features <- applyMSData(anaInfo, features@features, needIMS = TRUE, func = function(ana, path, backend, fTable)
    {
        openMSReadBackend(backend, path)
        
        fTable <- copy(fTable)
        
        fTableForPeaks <- if (!is.null(assignedMobilities))
            copy(fTable[!group %chin% assignedMobilities$group])
        else
            copy(fTable)
        if (!is.null(maxMSRTWindow))
        {
            fTableForPeaks[, retmin := pmax(retmin, ret - maxMSRTWindow)]
            fTableForPeaks[, retmax := pmin(retmax, ret + maxMSRTWindow)]
        }
        
        # NOTE: mzmin/mzmax may be too narrow here, hence use a user specified mz range
        EIMs <- getMobilograms(backend, fTableForPeaks$mz - mzWindow, fTableForPeaks$mz + mzWindow,
                               fTableForPeaks$retmin, fTableForPeaks$retmax, clusterMethod, IMSWindow,
                               minIntensityIMS, FALSE)
        names(EIMs) <- fTableForPeaks$ID
        EIMs <- lapply(EIMs, setDT)
        
        # pretend we have EICs so we can find peaks
        EICs <- lapply(EIMs, copy)
        EICs <- lapply(EICs, setnames, old = "mobility", new = "time")
        peaksList <- findPeaks(EICs, peaksParam, verbose = FALSE)
        
        peaksTable <- rbindlist(peaksList, idcol = "ims_parent_ID")
        
        if (!is.null(assignedMobilities))
        {
            am <- copy(assignedMobilities)
            am <- am[group %chin% fTable$group]
            am[, ims_parent_ID := fTable[match(am$group, group, nomatch = 0)]$ID][, group := NULL]
            peaksTable <- rbind(peaksTable, am, fill = TRUE)
        }
        
        peaksTable[, mobOrd := seq_len(.N), by = "ims_parent_ID"]
        
        mobNumCols <- c("mobility", "mobmin", "mobmax", "mob_area", "mob_intensity")
        if (nrow(peaksTable) == 0)
        {
            fTable[, ims_parent_ID := NA_character_]
            fTable[, (mobNumCols) := NA_real_]
        }
        else
        {
            fTable[, ord := seq_len(.N)]
            
            setnames(peaksTable, c("ret", "retmin", "retmax", "area", "intensity"), mobNumCols, skip_absent = TRUE)
            # NOTE: we subset columns here to remove any algo specific columns that may also be present in the feature
            # table (UNDONE?)
            peaksTable <- subsetDTColumnsIfPresent(peaksTable, c(mobNumCols, "mobOrd", "ims_parent_ID"))
            
            peaksTable[, ID := appendMobToName(ims_parent_ID, mobility)]
            
            # add feature data
            peaksTable <- merge(peaksTable, fTable, by.x = "ims_parent_ID", by.y = "ID", sort = FALSE)
            setcolorder(peaksTable, names(fTable))
            
            # merge mobility features
            fTable <- rbind(fTable, peaksTable, fill = TRUE)
            
            setorderv(fTable, c("ord", "mobOrd"), na.last = FALSE)
            fTable[, c("ord", "mobOrd") := NULL]
        }
        
        doProgress()
        
        return(fTable)
    })
    
    mobFeatsN <- sum(sapply(features@features, function(ft) sum(!is.na(ft$mobility))))
    printf("Assigned %d mobility features.\n", mobFeatsN)
    
    return(features)
}

# UNDONE: make this an exported method?
reintegrateFeatures <- function(features, RTWindow, calcArea, peaksParam, fallbackEIC, onlyMob, parallel)
{
    anaInfo <- analysisInfo(features)
    cacheDB <- openCacheDBScope()
    
    printf("Loading EICs...\n")
    allEICs <- getEICsForFeatures(features, RTWindow, onlyMob, compress = FALSE, cacheDB = cacheDB)
    
    if (!is.null(peaksParam))
    {
        # UNDONE: make withBP configurable?
        peaksList <- findPeaksInEICs(allEICs, peaksParam, withBP = FALSE, parallel = parallel, cacheDB = cacheDB)
        peaksList <- Map(peaksList, featureTable(features), f = function(anaPLs, ft)
        {
            # filter out peaks outside original retmin/retmax
            anaPLs <- anaPLs[numGTE(ret, ft[match(EIC_ID, ID)]$retmin) & numLTE(ret, ft[match(EIC_ID, ID)]$retmax)]
            # filter out all peaks for EICs with >1 result
            anaPLs[, N := .N, by = "EIC_ID"]
            return(anaPLs[N == 1][, N := NULL])
        })
    }
    else
        peaksList <- vector("list", nrow(anaInfo))
    
    updatedFeatsFromEICs <- updatedFeatsFromPeaks <- notAssigned <- 0
    features@features <- Map(featureTable(features), peaksList, allEICs, f = function(ft, pt, eics)
    {
        ft <- copy(ft)
        if (!is.null(pt))
        {
            cols <- c("ret", "retmin", "retmax", "area", "intensity")
            ft[pt, (cols) := mget(paste0("i.", cols)), on = c(ID = "EIC_ID")]
            peakIDs <- pt$EIC_ID
            updatedFeatsFromPeaks <<- updatedFeatsFromPeaks + sum(ft$ID %in% pt$EIC_ID)
        }
        else
            peakIDs <- character()
        
        if (fallbackEIC)
        {
            # update those not assigned by a peak from EICs
            doIDs <- ft$ID[(!onlyMob | (!is.null(ft[["mobility"]]) & !is.na(ft$mobility))) & !ft$ID %chin% peakIDs]
        
            ft[ID %chin% doIDs, intensity := mapply(ret, eics[doIDs], FUN = function(r, eic)
            {
                eic <- eic[eic$intensity > 0, ]
                if (nrow(eic) > 0) eic[which.min(abs(eic$time - r)), "intensity"] else 0
            })]
            ft[ID %chin% doIDs, area := mapply(retmin, retmax, eics[doIDs], FUN = function(rmin, rmax, eic)
            {
                wh <- eic$time %between% c(rmin, rmax)
                a <- sum(eic$intensity[wh])
                if (calcArea == "integrate")
                    a <- a * ((rmax - rmin) / sum(wh))
                return(a)
            })]
            
            updatedFeatsFromEICs <<- updatedFeatsFromEICs + length(doIDs)
        }
        else
        {
            # only keep those updated from a new peak or ignored in case onlyMob==T
            wh <- (onlyMob & (is.null(ft[["mobility"]]) | is.na(ft$mobility))) | ft$ID %chin% peakIDs
            notAssigned <<- notAssigned + sum(!wh)
            # UNDONE: use delete()? would give issues for eg XCMS objects
            ft <- ft[wh == TRUE]
        }
        return(ft)
    })
    
    if (fallbackEIC)
        printf("Re-integrated %d features (%d from newly found peaks and %d from EICs)\n",
               updatedFeatsFromEICs + updatedFeatsFromPeaks, updatedFeatsFromPeaks, updatedFeatsFromEICs)
    else
        printf("Re-integrated %d features and removed %d unassigned\n", updatedFeatsFromPeaks, notAssigned)
    
    return(features)
}

doFindMobilities <- function(fGroups, mobPeaksParam, mzWindow, IMSWindow, clusterMethod, minIntensityIMS,
                             maxMSRTWindow, chromPeaksParam, RTWindow, calcArea, fallbackEIC, assignedMobilities,
                             parallel)
{
    ac <- checkmate::makeAssertCollection()
    assertFindPeaksParam(mobPeaksParam, add = ac)
    aapply(checkmate::assertNumber, . ~ mzWindow + IMSWindow + minIntensityIMS + RTWindow, finite = TRUE,
           fixed = list(add = ac))
    checkmate::assertChoice(clusterMethod, c("bin", "distance", "hclust"), add = ac)
    checkmate::assertNumber(maxMSRTWindow, lower = 1, finite = TRUE, null.ok = TRUE, add = ac)
    assertFindPeaksParam(chromPeaksParam, null.ok = TRUE, add = ac)
    checkmate::assertChoice(calcArea, c("integrate", "sum"), add = ac)
    aapply(checkmate::assertFlag, . ~ fallbackEIC + parallel, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    if (length(fGroups) == 0)
        return(fGroups) # nothing to do...
    
    fTable <- featureTable(fGroups)
    hash <- makeHash(fGroups, mobPeaksParam, mzWindow, IMSWindow, clusterMethod, minIntensityIMS, maxMSRTWindow,
                     chromPeaksParam, RTWindow, calcArea, fallbackEIC)
    cd <- loadCacheData("findMobilities", hash)
    if (!is.null(cd))
        return(cd)
    
    fGroups@features <- assignFeatureMobilities(fGroups@features, mobPeaksParam, mzWindow, IMSWindow, clusterMethod,
                                                minIntensityIMS, maxMSRTWindow, assignedMobilities)
    fGroups@features <- reintegrateFeatures(fGroups@features, RTWindow, calcArea, chromPeaksParam, fallbackEIC, TRUE,
                                            parallel)
    
    fTable <- featureTable(fGroups)
    
    printf("Clustering mobilities... ")
    
    # cluster features within original fGroups with similar mobilities together    
    fTableAll <- rbindlist(fTable, idcol = "analysis")
    fTableAll[!is.na(mobility), gClust := {
        if (.N == 0)
            numeric()
        else if (.N == 1)
            1
        else
        {
            hc <- fastcluster::hclust(dist(mobility))
            cutree(hc, h = IMSWindow)
        }
    }, by = "group"]
    printf("Done!\n")
    
    printf("Updating feature group data... ")
    
    # prepare group info
    gMobInfo <- fTableAll[, .(mobility = mean(mobility)), by = c("group", "gClust")]
    setnames(gMobInfo, "group", "ims_parent_group")
    gMobInfo[, group := fifelse(!is.na(mobility), appendMobToName(ims_parent_group, mobility), ims_parent_group)]
    
    # update features
    setnames(fTableAll, "group", "ims_parent_group") # UNDONE: better colname
    fTableAll[gMobInfo, group := i.group, on = c("ims_parent_group", "gClust")]
    featureTable(fGroups) <- split(fTableAll[, -"gClust"], by = "analysis", keep.by = FALSE)
    
    # update gInfo
    gInfo <- copy(groupInfo(fGroups))
    setnames(gInfo, "group", "ims_parent_group")
    gInfo <- merge(gInfo, gMobInfo, by = "ims_parent_group", sort = FALSE)
    setcolorder(gInfo, c("group", "ret", "mz", "mobility", "ims_parent_group"))
    fGroups@groupInfo <- gInfo[, -"gClust"]
    
    # re-fill group table
    fTablePerGroup <- split(fTableAll, by = "group")
    fTablePerGroup <- fTablePerGroup[gInfo$group]
    anaInfo <- analysisInfo(fGroups)
    gTable <- data.table()
    gTable[, (gInfo$group) := lapply(fTablePerGroup, function(ft)
    {
        ints <- numeric(nrow(anaInfo))
        ints[match(ft$analysis, anaInfo$analysis)] <- ft$intensity
        return(ints)
    })]
    fGroups@groups <- gTable
    
    # update ftindex
    fGroups <- reGenerateFTIndex(fGroups)
    
    for (sl in c("groupQualities", "groupScores", "ISTDs", "ISTDAssignments", "annotations", "concentrations",
                 "toxicities"))
    {
        d <- slot(fGroups, sl)
        if (length(d) > 0)
        {
            warning("Clearing all data from ", sl, call. = FALSE)
            slot(fGroups, sl) <- if (is.data.table(d)) data.table() else list()
        }
    }
    
    printf("Done!\n")
    
    saveCacheData("findMobilities", fGroups, hash)
    
    return(fGroups)
}
