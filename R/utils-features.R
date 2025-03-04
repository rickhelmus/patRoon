# SPDX-FileCopyrightText: 2016 - 2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

makeFGroupName <- function(id, ret, mz) sprintf("M%.0f_R%.0f_%d", mz, ret, id)
makeIMSFGroupName <- function(id, ret, mz, mob) sprintf("M%.0f_R%.0f_M%.2f_%d", mz, ret, mob, id) # UNDONE: does this scale well for non-Bruker data?
appendMobToName <- function(n, mob) make.unique(sprintf("%s_I%.2f", n, mob)) # UNDONE: does this scale well for non-Bruker data?

showAnaInfo <- function(anaInfo)
{
    replicates <- unique(anaInfo$replicate)
    blGroups <- unique(anaInfo$blank)
    printf("Analyses: %s (%d total)\n", getStrListWithMax(anaInfo$analysis, 6, ", "), nrow(anaInfo))
    printf("Replicates: %s (%d total)\n", getStrListWithMax(replicates, 8, ", "), length(replicates))
    printf("Replicates used as blank: %s (%d total)\n", getStrListWithMax(blGroups, 8, ", "), length(blGroups))
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

doSubsetFeaturesByAna <- function(obj, i, ..., reorder, env)
{
    # NOTE: ni is passed as dots here, as otherwise NSE will _not_ work!
    if (!missing(..1) > 0) # check ..1: https://stackoverflow.com/a/25396892
    {
        if (!missing(i))
            stop("Cannot simulatenously specify i and ni arguments.", call. = FALSE)
        anaInfo <- analysisInfo(obj)
        # From https://stackoverflow.com/a/62043225
        anaInfo <- eval(substitute(anaInfo[.i], list(.i = substitute(..., env = env))))
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

doFGroupsFilter <- function(fGroups, what, hashParam, func, cacheCateg = what, applyIMS = "both", verbose = TRUE)
{
    if (verbose)
    {
        printf("Applying %s filter... ", what)
        oldFCount <- length(getFeatures(fGroups)); oldGCount <- length(fGroups)
    }
    
    cacheName <- sprintf("filterFGroups_%s", cacheCateg)
    hash <- makeHash(fGroups, what, hashParam, applyIMS)
    ret <- loadCacheData(cacheName, hash)
    if (is.null(ret))
    {
        ret <- if (length(fGroups) > 0) func(fGroups) else NULL
        if (!is.null(ret))
        {
            if (applyIMS != "both" && hasMobilities(fGroups))
            {
                # only remove mobility features/feature groups or their parents
                gInfoKeep <- groupInfo(fGroups)
                gInfoKeep <- if (applyIMS) gInfoKeep[is.na(mobility)] else gInfoKeep[!is.na(mobility)]
                
                if (!is.null(ret[["delete"]]) && !is.null(ret[["delete"]][["j"]]))
                {
                    if (is.character(ret$delete$j))
                        ret$delete$j <- setdiff(ret$delete$j, gInfoKeep$group)
                    else if (checkmate::testAtomicVector(ret$delete$j))
                        ret$delete$j <- setdiff(names(fGroups)[ret$delete$j], gInfoKeep$group)
                    else # DT
                    {
                        ret$delete$j <- copy(ret$delete$j)
                        ret$delete$j[, (gInfoKeep$group) := FALSE]
                    }
                }
                if (!is.null(ret[["subset"]]) && !is.null(ret[["subset"]][["j"]]))
                {
                    if (is.character(ret$subset$j))
                        ret$subset$j <- intersect(ret$subset$j, gInfoKeep$group)
                    else # logical/inds
                        ret$subset$j <- intersect(names(fGroups)[ret$subset$j], gInfoKeep$group)
                }
            }
            
            saveCacheData(cacheName, ret, hash)
        }
    }

    fGroupsFiltered <- getFilteredFGroups(fGroups, ret)
    
    if (verbose)
    {
        newFCount <- length(getFeatures(fGroupsFiltered)); newGCount <- length(fGroupsFiltered)
        printf("Done! Filtered %d (%.2f%%) features and %d (%.2f%%) feature groups. Remaining: %d features in %d groups.\n",
               oldFCount - newFCount, if (oldFCount > 0) (1 - (newFCount / oldFCount)) * 100,
               oldGCount - newGCount, if (oldGCount > 0) (1 - (newGCount / oldGCount)) * 100,
               newFCount, newGCount)
    }
    
    return(fGroupsFiltered)
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
    return(ret)
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

getDefEIXParams <- function()
{
    list(
        topMost = NULL,
        topMostByReplicate = FALSE,
        onlyPresent = TRUE,
        mzExpWindow = defaultLim("mz", "very_narrow"),
        mobExpWindow = defaultLim("mobility", "very_narrow"),
        mzExpIMSWindow = defaultLim("mz", "medium"),
        minIntensityIMS = 25,
        setsAdductPos = "[M+H]+",
        setsAdductNeg = "[M-H]-"
    )
}

filterEIXs <- function(EIXs, fGroups, analysis = NULL, groupName = NULL, topMost = NULL, topMostByReplicate = FALSE,
                       onlyPresent = FALSE)
{
    if (!is.null(analysis))
        EIXs <- EIXs[names(EIXs) %chin% analysis]
    else
        analysis <- analyses(fGroups)
    if (!is.null(groupName))
        EIXs <- lapply(EIXs, function(e) subListAttr(e, names(e) %chin% groupName))
    else
        groupName <- names(fGroups)

    if (onlyPresent)
    {
        gTable <- groupTable(fGroups)
        EIXs <- Map(names(EIXs), EIXs, f = function(ana, aeix)
        {
            anaInd <- match(ana, analyses(fGroups))
            absentFGs <- names(gTable)[gTable[anaInd] == 0]
            return(subListAttr(aeix, !names(aeix) %chin% absentFGs))
        })
    }
    
    if (!is.null(topMost))
    {
        gTable <- copy(groupTable(fGroups))
        gTable[, c("analysis", "replicate") := analysisInfo(fGroups)[, c("analysis", "replicate"), with = FALSE]]
        for (fg in groupName)
        {
            anasWithFG <- Map(names(EIXs), EIXs, f = function(ana, aeic) if (fg %chin% names(aeic)) ana else character())
            anasWithFG <- pruneList(anasWithFG, checkEmptyElements = TRUE)
            anasWithFG <- unlist(anasWithFG)
            tab <- gTable[analysis %chin% anasWithFG, c(fg, "analysis", "replicate"), with = FALSE]
            if (nrow(tab) > 0)
            {
                rmAnas <- if (topMostByReplicate)
                {
                    tab[, rank := frank(-get(fg), ties.method = "first"), by = "replicate"]
                    tab[rank > topMost]$analysis
                }
                else if (nrow(tab) > topMost)
                {
                    setorderv(tab, fg, order = -1L)
                    tab[seq(topMost + 1, nrow(tab))]$analysis
                }
                else
                    character()
                EIXs[rmAnas] <- lapply(EIXs[rmAnas], "[[<-", fg, NULL)
            }
        }
    }

    return(pruneList(EIXs, checkEmptyElements = TRUE, keepAttr = TRUE))
}

extendEIXInputTab <- function(tab, type, EIXParams)
{
    tab <- copy(tab)
    if (type == "EIC")
        tab[, c("retmin", "retmax") := .(retmin - EIXParams$window, retmax + EIXParams$window)]
    else # "EIM"
        tab[, c("mobmin", "mobmax") := .(mobmin - EIXParams$window, mobmax + EIXParams$window)]
    return(tab)
}

getEICsOREIMs <- function(obj, type, inputTab, EIXParams, ...)
{
    if (type == "EIC")
        doGetEICs(analysisInfo(obj), inputTab, EIXParams$mzExpIMSWindow, EIXParams$minIntensityIMS, ...)
    else # EIM
        doGetEIMs(analysisInfo(obj), inputTab, EIXParams$IMSWindow, EIXParams$clusterMethod, EIXParams$minIntensityIMS,
                  EIXParams$mzExpIMSWindow, ...)
}

setMethod("getFeatureEIXInputTab", "features", function(obj, type, EIXParams, selectFunc)
{
    return(lapply(featureTable(obj), function(tab)
    {
        if (!is.null(selectFunc))
            tab <- selectFunc(tab)
        
        tab <- copy(tab)
        # HACK: we keep group column for featureGroups method
        tab <- subsetDTColumnsIfPresent(tab, c("group", "analysis", "intensity", "retmin", "retmax", "mzmin",
                                               "mzmax", "mobmin", "mobmax"))

        if (type == "EIM" && is.null(tab[["mobmin"]]))
            tab[, c("mobmin", "mobmax") := NA_real_]

        tab <- extendEIXInputTab(tab, type, EIXParams)

        return(tab)
    }))
})

setMethod("getFeatureEIXInputTab", "featureGroups", function(obj, type, analysis, groupName, EIXParams)
{
    takeAnalysis <- analysis # copy name to workaround DT access below
    
    anaInfo <- analysisInfo(obj)[analysis %chin% takeAnalysis]
    
    topMost <- if (!is.null(EIXParams$topMost))
        min(EIXParams$topMost, nrow(anaInfo))
    else
        NULL
    
    featTab <- as.data.table(getFeatures(obj))
    
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
        if (!EIXParams$onlyPresent)
        {
            if (any(!analysis %chin% ret$analysis))
            {
                ftAllAna <- featTabSplitGrp[[fg]]
                dummy <- data.table(analysis = setdiff(analysis, ret$analysis), intensity = 0,
                                    retmin = min(ftAllAna$retmin), retmax = max(ftAllAna$retmax),
                                    mzmin = min(ftAllAna$mzmin) - EIXParams$mzExpWindow,
                                    mzmax = max(ftAllAna$mzmax) + EIXParams$mzExpWindow)
                if (!is.null(ftAllAna[["mobmin"]]))
                    dummy[, c("mobmin", "mobmax") := .(min(ftAllAna$mobmin) - EIXParams$mobExpWindow,
                                                       max(ftAllAna$mobmax) + EIXParams$mobExpWindow)]
                ret <- rbind(ret, dummy)
            }
        }
        
        if (!is.null(topMost))
        {
            if (EIXParams$topMostByReplicate)
            {
                ret[, replicate := anaInfo$replicate[match(analysis, anaInfo$analysis)]]
                ret[, rank := frank(-intensity, ties.method = "first"), by = "replicate"]
                ret <- ret[rank <= topMost]
            }
            else
            {
                setorderv(ret, "intensity", order = -1L)
                ret <- ret[seq_len(topMost)]
            }
        }
        
        if (type == "EIM" && is.null(ret[["mobmin"]]))
            ret[, c("mobmin", "mobmax") := NA_real_]
        ret <- extendEIXInputTab(ret, type, EIXParams)
        
        return(ret)
    }, simplify = FALSE))
})

setMethod("getFeatureEIXInputTab", "featureGroupsSet", function(obj, type, analysis, groupName, EIXParams)
{
    ret <- callNextMethod()
    
    anaInfo <- analysisInfo(obj)
    featTab <- as.data.table(getFeatures(obj))
    
    # HACK: since feature tables store the character form, it's easier to keep it all the same
    EIXParams$setsAdductPos <- as.character(EIXParams$setsAdductPos)
    EIXParams$setsAdductNeg <- as.character(EIXParams$setsAdductNeg)
    
    # 'ionize' m/zs
    return(Map(names(ret), ret, f = function(grp, ranges)
    {
        featTabGrp <- featTab[group == grp]
        ranges[, adduct := featTabGrp[match(ranges$analysis, analysis)]$adduct]
        
        if (!EIXParams$onlyPresent && any(is.na(ranges$adduct))) # adduct will be NA for 'missing' features
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
                    EIXParams$setsAdductPos
                else
                    EIXParams$setsAdductNeg
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

setMethod("getFeatureEIXs", "features", function(obj, type, analysis = analyses(obj), EIXParams, selectFunc = NULL, ...)
{
    if (length(obj) == 0)
        return(list())
    inputTab <- getFeatureEIXInputTab(obj, type, EIXParams, selectFunc)
    EIXs <- getEICsOREIMs(obj, type, inputTab, EIXParams, ...)
    EIXs <- Map(EIXs, featureTable(obj), f = function(eics, ft)
    {
        names(eics) <- if (!is.null(selectFunc)) selectFunc(ft)$ID else ft$ID
        return(eics)
    })
    return(pruneList(EIXs))
})

setMethod("getFeatureEIXs", "featuresSet", function(obj, type, ...)
{
    unsetFeatList <- sapply(sets(obj), unset, obj = obj, simplify = FALSE)
    EIXList <- sapply(unsetFeatList, getFeatureEIXs, type = type, ..., simplify = FALSE)
    EIXs <- unlist(EIXList, recursive = FALSE, use.names = FALSE) # use.names gives combined set/ana name, we just want ana
    names(EIXs) <- unlist(lapply(EIXList, names))
    EIXs <- EIXs[intersect(analyses(obj), names(EIXs))] # sync order
    return(EIXs)
})

setMethod("getFeatureEIXs", "featureGroups", function(obj, type, analysis = analyses(obj), groupName = names(obj),
                                                      EIXParams, ...)
{
    if (length(obj) == 0 || length(analysis) == 0 || length(groupName) == 0)
        return(list())
    
    takeAnalysis <- analysis # for DT subset below
    anaInfo <- analysisInfo(obj)[analysis %chin% takeAnalysis]
    
    inputTab <- getFeatureEIXInputTab(obj, type, analysis, groupName, EIXParams)
    inputTab <- split(rbindlist(inputTab, idcol = "group"), by = "analysis")
    inputTab <- inputTab[intersect(anaInfo$analysis, names(inputTab))] # sync order
    anaInfoEIXs <- anaInfo[analysis %in% names(inputTab)]
    
    EIXs <- getEICsOREIMs(obj, type, inputTab, EIXParams, ...)
    EIXs <- Map(EIXs, lapply(inputTab, "[[", "group"), f = setNames)
    
    return(pruneList(EIXs))
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
    
    gt <- groupTable(fGroups, areas = areas)[, resp$group, with = FALSE]
    if (hasMobilities(fGroups))
    {
        # for mobility features, we take parent intensities (if available) as mobility filtering may lead to reduced
        # intensities that currently are not taken into account when calculating response factors
        # UNDONE: make this configurable?
        gInfoChange <- groupInfo(fGroups)[!is.na(mobility) & group %chin% resp$group & ims_parent_group %chin% group]
        if (nrow(gInfoChange) > 0)
        {
            for (grpi in seq_along(gInfoChange))
                set(gt, j = gInfoChange$group[grpi], value = fGroups[[gInfoChange$ims_parent_group[grpi]]])
        }
    }
    
    gt <- transpose(gt)
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

findPeaksInEICs <- function(EICs, peakParams, withMobility, logPath, cacheDB = NULL)
{
    baseHash <- makeHash(peakParams)
    
    # NOTE: EICs must be named
    
    hash <- makeHash(baseHash, EICs)
    cd <- loadCacheData("peaksEIC", hash, cacheDB)
    if (!is.null(cd))
        return(cd)
    
    peaks <- findPeaks(EICs, TRUE, peakParams, logPath)
    peaks <- rbindlist(peaks, idcol = "EIC_ID")
    
    if (!is.null(EICs[[1]][["mz"]])) # EICs made with mode == "full"?
    {
        if (nrow(peaks) == 0)
        {
            peaks[, EIC_ID := character()]
            peaks[, c("retmin", "retmax", "ret", "area", "intensity", "mzmin", "mzmax", "mz", "mobmin", "mobmax", "mobility") := numeric()]
        }
        else
        {
            peaks[, c("mzmin", "mzmax", "mz", "mobmin", "mobmax", "mobility") := {
                eic <- EICs[[EIC_ID]]
                eic <- eic[numGTETol(eic$time, retmin) & numLTETol(eic$time, retmax), ]
                if (nrow(eic) == 0)
                    numeric(1)
                else
                {
                    if (is.null(eic[["mobility"]]))
                        eic$mobility <- eic$mobilityBP <- eic$mobmin <- eic$mobmax <- NA_real_
                    # UNDONE: also use mobility BP data?
                    list(min(eic$mzmin), max(eic$mzmax), weighted.mean(eic$mzBP, eic$intensity),
                         min(eic$mobmin), max(eic$mobmax),
                         weighted.mean(eic$mobilityBP, eic$intensity))
                }
            }, by = seq_len(nrow(peaks))]
        }
        
        # NOTE: we could also set mobilities after checking if data is available, but then we need to repeat the EIC subsetting above
        if (!withMobility || length(EICs) == 0 || is.null(EICs[[1]][["mobility"]]))
            peaks[, c("mobmin", "mobmax", "mobility") := NULL]
    }
    
    # make unique IDs
    peaks[, ID := make.unique(EIC_ID)]
    
    saveCacheData("peaksEIC", peaks, hash, cacheDB)
        
    return(peaks)
}

getMobilityCols <- function() c("mobility", "mobmin", "mobmax", "mob_area", "mob_intensity")
countMobilityFeatures <- function(feat) sum(sapply(featureTable(feat), function(ft) sum(!is.null(ft[["mobility"]]) & !is.na(ft$mobility))))

warnAndClearAssignedMobilities <- function(fGroups)
{
    if (hasMobilities(fGroups))
    {
        warning("Mobility features already have been assigned, these will be cleared now!", call. = FALSE)
        fGroups <- selectIMSFilter(fGroups, IMS = FALSE, verbose = FALSE)
    }
    return(fGroups)
}

doAssignFeatureMobilities <- function(fTable, mobTable)
{
    fTable <- copy(fTable)
    fTable[, ord := seq_len(.N)]
    
    mobTable <- copy(mobTable)
    mobTable[, mobOrd := seq_len(.N), by = "ims_parent_ID"]
    mobTable[, ID := appendMobToName(ims_parent_ID, mobility)]
    
    if (nrow(mobTable) == 0)
    {
        if (is.null(fTable[["mobility"]]))
        {
            # at least initialize columns
            mobNumCols <- intersect(getMobilityCols(), names(mobTable))
            fTable[, ims_parent_ID := NA_character_]
            fTable[, (mobNumCols) := NA_real_]
            fTable[, mob_assign_method := NA_character_]
        }
    }
    else
    {
        # add feature data: merge fTable, while making sure that no columns overlap
        mobTable <- merge(mobTable, fTable[, c("ID", setdiff(names(fTable), names(mobTable))), with = FALSE],
                          by.x = "ims_parent_ID", by.y = "ID", sort = FALSE)
        setcolorder(mobTable, names(fTable))
        
        # merge mobility features
        fTable <- rbind(fTable, mobTable, fill = TRUE)
        fTable[is.na(mobility), ims_parent_ID := NA_character_]
        
        setorderv(fTable, c("ord", "mobOrd"), na.last = FALSE)
    }
    
    fTable <- removeDTColumnsIfPresent(fTable, c("ord", "mobOrd"))
    
    return(fTable)
}

assignFeatureMobilitiesPeaks <- function(features, peakParams, IMSWindow, clusterMethod, minIntensityIMS,
                                         maxMSRTWindow)
{
    hash <- makeHash(features, peakParams, IMSWindow, clusterMethod, minIntensityIMS, maxMSRTWindow)
    cd <- loadCacheData("assignFeatureMobilitiesPeaks", hash)
    if (!is.null(cd))
        return(cd)
    
    printf("Finding mobilities for all features from mobilogram peaks...\n")
    oldCount <- countMobilityFeatures(features)
    
    EIMParams <- getDefEIMParams(window = NULLToZero(maxMSRTWindow), clusterMethod = clusterMethod,
                                 IMSWindow = IMSWindow, minIntensityIMS = minIntensityIMS)
    # skip any mobility features and IMS parents
    EIMSelFunc <- \(tab) if (is.null(tab[["mobility"]])) tab else tab[is.na(mobility) & !ID %chin% ims_parent_ID]
    allEIMs <- getFeatureEIXs(features, "EIM", EIXParams = EIMParams, selectFunc = EIMSelFunc, compress = FALSE)
    
    mobNumCols <- getMobilityCols()
    
    features@features <- Map(features@features, allEIMs, analyses(features), f = function(fTable, EIMs, ana)
    {
        peaksTable <-  data.table()
        if (length(EIMs) > 0)
        {
            # pretend we have EICs so we can find peaks
            EIMs <- lapply(EIMs, setnames, old = "mobility", new = "time")
            peaksList <- findPeaks(EIMs, FALSE, peakParams, file.path("log", "assignMobilities", paste0("mobilogram_peaks-", ana, ".txt")))
            peaksTable <- rbindlist(peaksList, idcol = "ims_parent_ID")                
            setnames(peaksTable, c("ret", "retmin", "retmax", "area", "intensity"), mobNumCols, skip_absent = TRUE)
            # NOTE: we subset columns here to remove any algo specific columns that may also be present in the feature
            # table (UNDONE?)
            peaksTable <- subsetDTColumnsIfPresent(peaksTable, c(mobNumCols, "ims_parent_ID"))
        }
        if (length(peaksTable) == 0)
            peaksTable <- data.table()[, (mobNumCols) := numeric()][, ims_parent_ID := character()]

        peaksTable[, mob_assign_method := "peak"]
        return(doAssignFeatureMobilities(fTable, peaksTable))
    })

    printf("Assigned %d mobility features.\n", countMobilityFeatures(features) - oldCount)
    
    features@hasMobilities <- TRUE
    
    saveCacheData("assignFeatureMobilitiesPeaks", features, hash)
    
    return(features)
}

# UNDONE: make this an exported method?
reintegrateMobilityFeatures <- function(features, EICRTWindow, peakRTWindow, calcArea, minIntensityIMS, peakParams,
                                        fallbackEIC, parallel)
{
    cacheDB <- openCacheDBScope()
    hash <- makeHash(features, EICRTWindow, peakRTWindow, calcArea, minIntensityIMS, peakParams, fallbackEIC)
    cd <- loadCacheData("reintegrateMobilityFeatures", hash, cacheDB)
    if (!is.null(cd))
        return(cd)
    
    anaInfo <- analysisInfo(features)
    
    printf("Loading EICs...\n")
    EICSelFunc <- \(tab) tab[!is.null(tab[["mobility"]]) & !is.na(mobility)]
    allEICs <- getFeatureEIXs(features, type = "EIC", EIXParams = getDefEICParams(window = EICRTWindow,
                                                                                  minIntensityIMS = minIntensityIMS),
                              selectFunc = EICSelFunc, cacheDB = cacheDB)
    
    if (!is.null(peakParams))
    {
        peaksList <- doApply("Map", parallel, allEICs, featureTable(features), analyses(features), f = function(EICs, ft, ana)
        {
            peaks <- findPeaksInEICs(EICs, peakParams, withMobility = FALSE,
                                     logPath = file.path("log", "assignMobilities", paste0("reintegrate-", ana, ".txt")),
                                     cacheDB = cacheDB)
            # filter out peaks outside original retmin/retmax and with high RT deviation
            parFT <- ft[match(peaks$EIC_ID, ID)]
            peaks <- peaks[numGTE(ret, parFT$retmin) & numLTE(ret, parFT$retmax) & numLTE(abs(ret - parFT$ret), peakRTWindow)]
            # filter out all peaks for EICs with >1 result
            peaks[, N := .N, by = "EIC_ID"]
            
            doProgress()
            
            return(peaks[N == 1][, N := NULL])
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
            ft[pt, (c(cols, "mob_reintegr_method")) := c(mget(paste0("i.", cols)), list("peak")), on = c(ID = "EIC_ID")]
            peakIDs <- pt$EIC_ID
            updatedFeatsFromPeaks <<- updatedFeatsFromPeaks + sum(ft$ID %in% pt$EIC_ID)
        }
        else
            peakIDs <- character()
        
        ft[, keep := TRUE]
        
        if (fallbackEIC)
        {
            # update those not assigned by a peak from EICs
            doRows <- which(ft$ID %chin% names(eics) & !ft$ID %chin% peakIDs)

            ft[doRows, c("intensity", "area") := {
                eic <- eics[[ID]]
                eic <- eic[numGTETol(eic$time, retmin) & numLTETol(eic$time, retmax), ]
                
                eicnz <- eic[eic$intensity > 0, ]
                i <- if (nrow(eicnz) > 0) eicnz[which.min(abs(eicnz$time - ret)), "intensity"] else 0
                
                a <- 0
                if (nrow(eic) > 0)
                {
                    a <- sum(eic$intensity)
                    if (calcArea == "integrate")
                        a <- a * ((retmax - retmin) / nrow(eic))
                }
                
                list(i, a)
            }, by = seq_along(doRows)]
            ft[doRows, mob_reintegr_method := "EIC"]

            ft[doRows, keep := intensity > 0]
            updatedFeatsFromEICs <<- updatedFeatsFromEICs + ft[doRows][keep == TRUE, .N]
        }
        else
        {
            # only keep those updated from a new peak
            ft[, keep := ((is.null(ft[["mobility"]]) | is.na(mobility))) | ID %chin% peakIDs]
        }

        notAssigned <<- notAssigned + sum(!ft$keep)
        # UNDONE: use delete()? would give issues for eg XCMS objects
        ft <- ft[keep == TRUE, -"keep"]
        
        return(ft)
    })
    
    printf("Re-integrated %d features (%d from newly found peaks, %d from EICs and removed %d unassigned)\n",
           updatedFeatsFromEICs + updatedFeatsFromPeaks, updatedFeatsFromPeaks, updatedFeatsFromEICs, notAssigned)
    
    saveCacheData("reintegrateMobilityFeatures", features, hash)
    
    return(features)
}

clusterFGroupMobilities <- function(fGroups, IMSWindow, sets)
{
    # UNDONE: do something more polymorphic than hackish sets arg?
    
    hash <- makeHash(fGroups, IMSWindow)
    cd <- loadCacheData("clusterFGroupMobilities", hash)
    if (!is.null(cd))
        return(cd)
    
    fTableAll <- rbindlist(featureTable(fGroups), idcol = "analysis")
    if (sets)
        fTableAll[, set := analysisInfo(fGroups)$set[match(analysis, analyses(fGroups))]]
    
    printf("Clustering mobilities... ")
    # cluster features within original fGroups with similar mobilities together    
    byCols <- if (sets) c("group", "set") else "group"
    fTableAll[!is.na(mobility), IMSClust := {
        cl <- if (.N == 0)
            integer()
        else if (.N == 1)
            1
        else
        {
            hc <- fastcluster::hclust(dist(mobility))
            cutree(hc, h = IMSWindow)
        }
        paste0(paste0(unlist(.BY), collapse = "_"), "-", cl)
    }, by = byCols]
    printf("Done!\n")

    printf("Updating feature group data... ")
    
    # prepare group info
    gMobInfo <- fTableAll[, .(mobility = mean(mobility)), by = c("group", "IMSClust")]
    setnames(gMobInfo, "group", "ims_parent_group")
    gMobInfo[is.na(mobility), group := ims_parent_group]
    gMobInfo[!is.na(mobility), group := appendMobToName(ims_parent_group, mobility)]
    
    # update features
    setnames(fTableAll, "group", "ims_parent_group") # UNDONE: better colname
    fTableAll[gMobInfo, group := i.group, on = c("ims_parent_group", "IMSClust")]
    fTableAllClean <- removeDTColumnsIfPresent(fTableAll, c("IMSClust", "set", "ims_parent_group"))
    fTable <- split(fTableAllClean, by = "analysis", keep.by = FALSE)
    # NOTE: the above will not restore any empty feature tables
    missingAna <- setdiff(analyses(fGroups), names(fTable))
    fTable[missingAna] <- rep(list(fTableAllClean[0]), length(missingAna)) # get empty table with all cols
    fTable <- fTable[analyses(fGroups)] # restore order
    featureTable(fGroups) <- fTable
    
    # update gInfo
    gInfo <- copy(groupInfo(fGroups))
    setnames(gInfo, "group", "ims_parent_group")
    gInfo <- merge(gInfo, gMobInfo, by = "ims_parent_group", sort = FALSE)
    setcolorder(gInfo, c("group", "ret", "mz", "mobility", "ims_parent_group"))
    gInfo[is.na(mobility), ims_parent_group := NA_character_]
    fGroups@groupInfo <- gInfo[, -"IMSClust"]
    
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
    
    for (sl in c("groupQualities", "groupScores", "annotations", "concentrations", "toxicities"))
    {
        d <- slot(fGroups, sl)
        if (length(d) > 0)
        {
            printf("NOTE: copying parent data from %s\n", sl)
            
            # slots are all data.tables with group column
            d <- rbindlist(lapply(split(d, seq_len(nrow(d))), function(r)
            {
                og <- gInfo[!is.na(mobility) & ims_parent_group == r$group]$group
                rbind(r, data.table(group = og, r[, -"group"]))
            }))
            
            # HACK: for sets workflows, ensure that we only have keep the relevant mobility fGroup annotation
            if (sets && sl == "annotations")
            {
                fGroupsPerSet <- fTableAll[, .(group = unique(group)), by = "set"]
                keep <- mapply(d$group, d$set, FUN = function(g, s) g %chin% fGroupsPerSet[set == s]$group)
                d <- d[keep == TRUE]
            }
            
            slot(fGroups, sl) <- d
        }
    }
    
    if (length(internalStandardAssignments(fGroups)) > 0)
    {
        # copy parent IS assignments
        
        addMobISTDAssigns <- function(ISA)
        {
            gi <- gInfo[group %chin% names(ISA) | ims_parent_group %chin% names(ISA)]
            mobFGs <- gi[!is.na(mobility)]$group
            mobFGParents <- gi[!is.na(mobility)]$ims_parent_group
            setNames(ISA[mobFGParents], mobFGs)
        }

        if (sets)
            fGroups@ISTDAssignments <- lapply(internalStandardAssignments(fGroups),
                                              function(ISA) c(ISA, addMobISTDAssigns(ISA)))
        else
            fGroups@ISTDAssignments <- c(fGroups@ISTDAssignments, addMobISTDAssigns(internalStandardAssignments(fGroups)))
    }
    
    saveCacheData("clusterFGroupMobilities", fGroups, hash)
    
    printf("Done!\n")
    
    return(fGroups)
}

assignFGroupsCCS <- function(fGroups, CCSParams)
{
    if (!hasMobilities(fGroups))
        stop("Cannot calculate CCS values: feature groups are without mobility assignments", call. = FALSE)
    
    charges <- if (nrow(annotations(fGroups)) > 0)
    {
        # NOTE: there is at most 1 mobility fGroup per set, so we don't have to worry about multiple annotations in sets workflows
        ann <- annotations(fGroups)[group %chin% groupInfo(fGroups)[!is.na(mobility)]$group]
        unAdd <- unique(ann$adduct)
        addCharges <- sapply(unAdd, function(a) as.adduct(a)@charge)
        setNames(addCharges[ann$adduct], ann$group)
    }
    else
        setNames(rep(CCSParams$defaultCharge, length(fGroups)), names(fGroups))
    
    fGroups@features@features <- lapply(featureTable(fGroups), function(ft)
    {
        ft <- copy(ft)
        ft[!is.na(mobility), CCS := convertMobilityToCCS(mobility, mz, CCSParams, charges[group])]
        return(ft[])
    })
    
    fGroups@groupInfo <- copy(groupInfo(fGroups))
    fGroups@groupInfo[!is.na(mobility), CCS := convertMobilityToCCS(mobility, mz, CCSParams, charges[group])]
    
    return(fGroups)
}

# similar to selectIMSFilter(), but for features
selectIMSFilterFeatures <- function(features, IMS)
{
    if (IMS == "both" || !hasMobilities(features))
        return(features)
    
    features <- delete(features, j = function(fl, ...)
    {
        if (isTRUE(IMS))
            is.na(fl$mobility)
        else if (isFALSE(IMS))
            !is.na(fl$mobility)
        else # "maybe"
            !is.na(fl$mobility) & fl$ims_parent_ID %chin% fl$ID
    })
    
    if (isFALSE(IMS))
        features <- clearMobilities(features)
    
    return(features)
}

expandTableForIMSFGroups <- function(tab, gInfo)
{
    # map all fGroups to copy: these are all mobility fGroups that are not orphans
    gInfo <- gInfo[group %chin% tab$group | ims_parent_group %chin% tab$group]
    tab <- copy(tab)
    tab[, orderOrig := seq_len(.N)]
    imspars <- fifelse(is.na(gInfo$ims_parent_group), gInfo$group, gInfo$ims_parent_group)
    tab <- tab[match(imspars, group)] # expand & copy
    tab[, group := gInfo$group]
    setorderv(tab, "orderOrig")
    tab[, orderOrig := NULL]
    return(tab[])
}

assignTabIMSDeviations <- function(tab, gInfo)
{
    # used by screening / compounds
    
    tab <- copy(tab)
    tab[, d_mob := gInfo$mobility[match(group, gInfo$group)] - mobility]
    tab[, d_mob_rel := d_mob / mobility]
    if (is.null(gInfo[["CCS"]]))
        tab[, d_CCS := NA_real_]
    else
        tab[, d_CCS := gInfo$CCS[match(group, gInfo$group)] - CCS]
    tab[, d_CCS_rel := d_CCS / CCS]
    
    return(tab[])
}
