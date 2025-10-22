# SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

getMobilityCols <- function() c("mobility", "mobmin", "mobmax", "mob_area", "mob_intensity")
countMobilityFeatures <- function(feat) sum(sapply(featureTable(feat), function(ft) sum(!is.null(ft[["mobility"]]) & !is.na(ft$mobility))))

checkUnsupportedIMS <- function(feat, algorithm)
{
    if (hasMobilities(feat))
        stop(sprintf("The '%s' algorithm does not support ion mobility data", algorithm), call. = FALSE)
}

checkAssignedMobilityFGroups <- function(fGroups)
{
    if (hasMobilities(fGroups))
    {
        if (all(!is.na(groupInfo(fGroups)$mobility)))
            stop("There are no feature groups without mobility assignments available for which mobility features can be assigned.", call. = FALSE)
        
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

# split into separate function for future parallelization
doFindPeaksForMobilities <- function(EIMs, ana, peakParams)
{
    mobNumCols <- getMobilityCols()
    
    peaksTable <- data.table()
    if (length(EIMs) > 0)
    {
        # pretend we have EICs so we can find peaks
        EIMs <- lapply(EIMs, \(e) { colnames(e)[1] <- "time"; e })
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
    
    return(peaksTable)
    # return(doAssignFeatureMobilities(fTable, peaksTable))
}

assignFeatureMobilitiesPeaks <- function(features, peakParams, EIMParams, parallel)
{
    hash <- makeHash(features, peakParams, EIMParams)
    cd <- loadCacheData("assignFeatureMobilitiesPeaks", hash)
    if (!is.null(cd))
        return(cd)
    
    EIMParams$topMost <- NULL; EIMParams$onlyPresent <- TRUE
    
    printf("Finding mobilities for all features from mobilogram peaks...\n")
    oldCount <- countMobilityFeatures(features)
    
    # skip any mobility features and IMS parents
    EIMSelFunc <- \(tab) if (is.null(tab[["mobility"]])) tab else tab[is.na(mobility) & !ID %chin% ims_parent_ID]
    allEIMs <- getFeatureEIXs(features, "EIM", EIXParams = EIMParams, selectFunc = EIMSelFunc, compress = FALSE)
    
    if (identical(parallel, "maybe") && peakParams$algorithm == "piek")
        parallel <- FALSE # piek is already parallelized internally
    peaksList <- doMap(parallel, allEIMs, analyses(features), f = patRoon:::doFindPeaksForMobilities,
                       MoreArgs = list(peakParams = peakParams))
    features@features <- Map(featureTable(features), peaksList, f = doAssignFeatureMobilities)
    printf("Assigned %d mobility features.\n", countMobilityFeatures(features) - oldCount)
    
    features@hasMobilities <- TRUE
    
    saveCacheData("assignFeatureMobilitiesPeaks", features, hash)
    
    return(features)
}

# split into separate function for future parallelization
doFindPeaksForReintegration <- function(EICs, peakParams, peakRTWindow, ft, ana, cacheDB)
{
    peaks <- findPeaksInEICs(EICs, peakParams, withMobility = FALSE, calcStats = FALSE, assignRTWindow = 0,
                             logPath = file.path("log", "assignMobilities", paste0("reintegrate-", ana, ".txt")),
                             cacheDB = cacheDB)
    # filter out peaks outside original retmin/retmax and with high RT deviation
    parFT <- ft[match(peaks$EIC_ID, ID)]
    peaks <- peaks[numGTE(ret, parFT$retmin) & numLTE(ret, parFT$retmax) & numLTE(abs(ret - parFT$ret), peakRTWindow)]
    # filter out all peaks for EICs with >1 result
    peaks[, N := .N, by = "EIC_ID"]
    
    return(peaks[N == 1][, N := NULL])
}


# UNDONE: make this an exported method?
reintegrateMobilityFeatures <- function(features, peakParams, EICParams, peakRTWindow, fallbackEIC, calcArea, parallel)
{
    cacheDB <- openCacheDBScope()
    hash <- makeHash(features, peakParams, EICParams, peakRTWindow, fallbackEIC, calcArea)
    cd <- loadCacheData("reintegrateMobilityFeatures", hash, cacheDB)
    if (!is.null(cd))
        return(cd)
    
    anaInfo <- analysisInfo(features)
    
    EICParams$topMost <- NULL; EICParams$onlyPresent <- TRUE
    
    printf("Loading EICs...\n")
    EICSelFunc <- \(tab) tab[!is.null(tab[["mobility"]]) & !is.na(mobility)]
    allEICs <- getFeatureEIXs(features, type = "EIC", EIXParams = EICParams, selectFunc = EICSelFunc, cacheDB = cacheDB)
    
    if (!is.null(peakParams))
    {
        if (identical(parallel, "maybe") && peakParams$algorithm == "piek")
            parallel <- FALSE # piek is already parallelized internally
        peaksList <- doMap(parallel, allEICs, featureTable(features), analyses(features),
                           f = patRoon:::doFindPeaksForReintegration,
                           MoreArgs = list(peakParams = peakParams, peakRTWindow = peakRTWindow,
                                           cacheDB = if (!parallel) cacheDB))
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
        
        # only keep IMS parents and those updated from a new peak
        ft[, keep := ((is.null(ft[["mobility"]]) | is.na(mobility))) | ID %chin% peakIDs]
        
        if (fallbackEIC)
        {
            # update those not assigned by a peak from EICs, only consider those with mobility assignment from EIMs
            # (i.e. not suspects, as they need some kind of verification)
            doRows <- which(ft$ID %chin% names(eics) & !ft$ID %chin% peakIDs & ft$mob_assign_method == "peak")
            ft[doRows, c("intensity", "area") := {
                eic <- eics[[ID]]
                eic <- eic[numGTETol(eic[ , "time"], retmin) & numLTETol(eic[, "time"], retmax), , drop = FALSE]
                
                eicnz <- eic[eic[, "intensity"] > 0, , drop = FALSE]
                i <- if (nrow(eicnz) > 0) eicnz[which.min(abs(eicnz[, "time"] - ret)), "intensity"] else 0
                
                a <- 0
                if (nrow(eic) > 0)
                {
                    a <- sum(eic[, "intensity"])
                    if (calcArea == "integrate")
                        a <- a * ((retmax - retmin) / nrow(eic))
                }
                
                list(i, a)
            }, by = seq_along(doRows)]
            ft[doRows, mob_reintegr_method := "EIC"]
            
            ft[doRows, keep := intensity > 0] # also keep updated features
            updatedFeatsFromEICs <<- updatedFeatsFromEICs + ft[doRows][keep == TRUE, .N]
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

groupFTableMobilities <- function(feat, IMSWindow, byGroup)
{
    fTableAll <- as.data.table(feat)
    anaInfo <- analysisInfo(feat)
    reps <- replicates(feat)
    
    byCols <- if (byGroup) "group" else character()
    # HACK
    if (!is.null(fTableAll[["set"]]))
        byCols <- c(byCols, "set")
    
    fTableAll[, IMSGroup := NA_character_]
    fTableAll[, replicate := anaInfo$replicate[match(analysis, anaInfo$analysis)]]
    fTableAll[!is.na(mobility), IMSGroup := {
        cl <- if (.N == 0)
            integer()
        else if (.N == 1)
            1L
        else
        {
            # do a greedy grouping with just mobilities, setting dummy values for RT and mz (we don't want to regroup these)
            getGroupIDs(rep(100, .N), rep(100, .N), mobility, intensity, match(analysis, anaInfo$analysis),
                        match(replicate, reps), 5, 1, IMSWindow, c(retention = 1, mz = 1, mobility = 1, intensity = 1))
        }
        paste0(paste0(unlist(.BY), collapse = "_"), "-", cl)
    }, by = byCols]
    fTableAll[, replicate := NULL]
    
    return(fTableAll)
}

updateFGroupsForMobilities <- function(fGroups, IMSWindow, sets)
{
    # UNDONE: do something more polymorphic than hackish sets arg?
    
    hash <- makeHash(fGroups, IMSWindow)
    cd <- loadCacheData("updateFGroupsForMobilities", hash)
    if (!is.null(cd))
        return(cd)
    
    # cluster features within original fGroups with similar mobilities together    
    printf("Grouping mobilities... ")
    fTableAll <- groupFTableMobilities(getFeatures(fGroups), IMSWindow, byGroup = TRUE)
    printf("Done!\n")
    
    printf("Updating feature group data... ")
    
    # prepare group info
    gMobInfo <- fTableAll[, .(mobility = mean(mobility)), by = c("group", "IMSGroup")]
    setnames(gMobInfo, "group", "ims_parent_group")
    gMobInfo[is.na(mobility), group := ims_parent_group]
    gMobInfo[!is.na(mobility), group := appendMobToName(ims_parent_group, mobility)]
    
    # update features
    setnames(fTableAll, "group", "ims_parent_group") # UNDONE: better colname
    fTableAll[gMobInfo, group := i.group, on = c("ims_parent_group", "IMSGroup")]
    fTableAllClean <- removeDTColumnsIfPresent(fTableAll, c("IMSGroup", "set", "ims_parent_group"))
    fTable <- split(fTableAllClean, by = "analysis", keep.by = FALSE)
    # NOTE: the above will not restore any empty feature tables
    missingAna <- setdiff(analyses(fGroups), names(fTable))
    fTable[missingAna] <- rep(list(fTableAllClean[0, -"analysis"]), length(missingAna)) # get empty table with all cols
    fTable <- fTable[analyses(fGroups)] # restore order
    featureTable(fGroups) <- fTable
    
    # update gInfo
    gInfo <- copy(groupInfo(fGroups))
    setnames(gInfo, "group", "ims_parent_group")
    gInfo <- merge(gInfo, gMobInfo, by = "ims_parent_group", sort = FALSE)
    setcolorder(gInfo, c("group", "ret", "mz", "mobility", "ims_parent_group"))
    gInfo[is.na(mobility), ims_parent_group := NA_character_]
    fGroups@groupInfo <- gInfo[, -"IMSGroup"]
    
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
            if (nrow(d) == 0)
                next # skip, otherwise we end up with a numm DT (ie columns removed)
            
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
    
    saveCacheData("updateFGroupsForMobilities", fGroups, hash)
    
    printf("Done!\n")
    
    return(fGroups)
}

assignFGroupsCCS <- function(fGroups, CCSParams)
{
    if (!hasMobilities(fGroups))
        stop("Cannot calculate CCS values: feature groups are without mobility assignments", call. = FALSE)
    
    # HACK: may be do some polymorhic things some day...
    hasSets <- isFGSet(fGroups)
    
    # NOTE: there is at most 1 mobility fGroup per set, so we don't have to worry about multiple annotations in sets workflows
    ann <- annotations(fGroups)[group %chin% groupInfo(fGroups)[!is.na(mobility)]$group]
    
    grpCharges <- if (nrow(ann) > 0)
    {
        unAdd <- unique(ann$adduct)
        addCharges <- sapply(unAdd, function(a) as.adduct(a)@charge)
        setNames(addCharges[ann$adduct], ann$group)
    }
    else
        grpCharges <- setNames(rep(CCSParams$defaultCharge, length(fGroups)), names(fGroups))
    
    grpMZs <- if (hasSets)
        setNames(ann$ion_mz, ann$group)
    else
        grpMZs <- setNames(groupInfo(fGroups)$mz, names(fGroups))

    featMZCol <- if (hasSets) "ion_mz" else "mz"
    fGroups@features@features <- lapply(featureTable(fGroups), function(ft)
    {
        ft <- copy(ft)
        ft[!is.na(mobility), CCS := convertMobilityToCCS(mobility, get(featMZCol), CCSParams, grpCharges[group])]
        return(ft[])
    })
    
    fGroups@groupInfo <- copy(groupInfo(fGroups))
    fGroups@groupInfo[!is.na(mobility), CCS := convertMobilityToCCS(mobility, grpMZs[group], CCSParams, grpCharges[group])][]
    
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
    expanded <- rbindlist(lapply(seq_len(nrow(tab)), function(row)
    {
        rowTab <- tab[row]
        ipInds <- which(rowTab$group == imspars)
        expTab <- rowTab[rep(1, length(ipInds))]
        expTab[, group := gInfo$group[ipInds]]
        return(expTab)
    }))
    setorderv(expanded, "orderOrig")
    expanded[, orderOrig := NULL]
    return(expanded[])
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
