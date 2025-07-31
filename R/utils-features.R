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
    if (!missing(..1)) # check ..1: https://stackoverflow.com/a/25396892
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

setMethod("doGroupFeatures", "features", function(feat, grouper, groupAlgo, ..., verbose)
{
    # UNDONE: disabled: groupFeaturesIMS() doesn't really work...
    # if (!hasMobilities(feat))
    #     return(grouper(feat, ..., verbose = verbose))
    # return(groupFeaturesIMS(feat, grouper, groupAlgo, ..., IMSWindow = IMSWindow, verbose = verbose))
    return(grouper(feat, ..., verbose = verbose))
})

reGenerateFTIndex <- function(fGroups)
{
    gNames <- names(fGroups)
    fGroups@ftindex <- setnames(rbindlist(lapply(featureTable(fGroups),
                                                 function(ft) as.list(chmatch(gNames, ft$group, 0)))), gNames)
    return(fGroups)
}

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
        fTable <- featureTable(fGroups)
        EIXs <- Map(names(EIXs), EIXs, f = \(ana, aeix) subListAttr(aeix, names(aeix) %chin% fTable[[ana]]$group))
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
        doGetEICs(analysisInfo(obj), inputTab, EIXParams$gapFactor, EIXParams$mzExpIMSWindow, EIXParams$minIntensityIMS,
                  ...)
    else # EIM
        doGetEIMs(analysisInfo(obj), inputTab, EIXParams$minIntensityIMS, EIXParams$sgOrder, EIXParams$sgLength,
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
                                               "mzmax", "mobmin", "mobmax", "ret"))

        if (type == "EIM")
        {
            if (is.null(tab[["mobmin"]]))
                tab[, c("mobmin", "mobmax") := NA_real_]
            
            tab[, c("retmin", "retmax") := .(pmax(ret - EIXParams$maxRTWindow, retmin),
                                             pmin(ret + EIXParams$maxRTWindow, retmax))]
        }
        
        tab[, ret := NULL]
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
                                                   "mzmax", "mobmin", "mobmax", "ret"))

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
        
        if (type == "EIM")
        {
            if (is.null(ret[["mobmin"]]))
                ret[, c("mobmin", "mobmax") := NA_real_]
            
            ret[, c("retmin", "retmax") := .(pmax(ret - EIXParams$maxRTWindow, retmin),
                                             pmin(ret + EIXParams$maxRTWindow, retmax))]
        }
        
        ret[, ret := NULL]
        ret <- extendEIXInputTab(ret, type, EIXParams)
        
        return(ret)
    }, simplify = FALSE))
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
            for (grpi in seq_len(nrow(gInfoChange)))
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

findPeaksInEICs <- function(EICs, peakParams, withMobility, calcStats, logPath, cacheDB = NULL)
{
    baseHash <- makeHash(peakParams)
    
    # NOTE: EICs must be named
    
    hash <- makeHash(baseHash, EICs)
    cd <- loadCacheData("peaksEIC", hash, cacheDB)
    if (!is.null(cd))
        return(cd)
    
    peaks <- findPeaks(EICs, TRUE, peakParams, logPath)
    peaks <- rbindlist(peaks, idcol = "EIC_ID")

    if (nrow(peaks) == 0)
    {
        peaks[, EIC_ID := character()]
        peaks[, c("retmin", "retmax", "ret", "area", "intensity") := numeric()]
        if (calcStats)
            peaks[, c("mzmin", "mzmax", "mz") := numeric()]
        if (withMobility)
        {
            if (calcStats)
                peaks[, c("mobmin", "mobmax") := numeric()]
            peaks[, mobility := numeric()]
        }
    }
    else if (calcStats)
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
                     weighted.mean(eic$mobility, eic$intensity))
            }
        }, by = seq_len(nrow(peaks))]
        
        # NOTE: we could also set mobilities after checking if data is available, but then we need to repeat the EIC subsetting above
        if (!withMobility || length(EICs) == 0 || is.null(EICs[[1]][["mobility"]]))
            peaks[, c("mobmin", "mobmax", "mobility") := NULL]
    }
    
    # make unique IDs
    peaks[, ID := make.unique(EIC_ID)]
    
    saveCacheData("peaksEIC", peaks, hash, cacheDB)
        
    return(peaks)
}
