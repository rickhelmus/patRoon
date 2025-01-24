# SPDX-FileCopyrightText: 2016 - 2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include utils.R
#' @include utils-compounds.R
NULL

isScreening <- function(fGroups) inherits(fGroups, c("featureGroupsScreening", "featureGroupsScreeningSet"))
isSuspAnnotated <- function(fGroups) isScreening(fGroups) && !is.null(screenInfo(fGroups)[["estIDLevel"]])

suspMetaDataCols <- function() c("name", "rt", "name_orig", "mz", "mobility", "CCS", "SMILES", "InChI", "InChIKey",
                                 "formula", "neutralMass", "molNeutralized", "adduct", "fragments_mz",
                                 "fragments_formula")
suspAnnCols <- function() c("formRank", "compRank", "annSimForm", "annSimComp", "annSimBoth", "maxFrags",
                            "maxFragMatches", "maxFragMatchesRel", "estIDLevel")

# UNDONE: remove this function? If not and it changes, then update assignSetsIDLs()
getAllSuspCols <- function(targetCols, allCols, mConsNames) getMergedConsCols(targetCols, allCols, mConsNames)

doFGroupScreeningClearMobilities <- function(obj)
{
    obj <- callNextMethod()
    cols <- getMergedConsCols(c("mobility", "CCS", "d_mob", "d_mob_rel", "d_CCS", "d_CCS_rel"), names(obj@screenInfo),
                              mergedConsensusNames(obj))
    if (length(cols) > 0)
    {
        scr <- copy(screenInfo(obj))
        scr[, (cols) := NULL]
        obj@screenInfo <- scr[]
    }
    return(obj)
}

doScreeningShow <- function(obj)
{
    printf("Suspects: %s (%d hits total)\n", getStrListWithMax(unique(screenInfo(obj)$name), 6, ", "),
           nrow(screenInfo(obj)))
    printf("Suspects annotated: %s\n", if (isSuspAnnotated(obj)) "yes" else "no")
}

calcSuspMZs <- function(suspects, adduct)
{
    suspects <- copy(suspects)
    
    if (is.null(suspects[["mz"]]))
        suspects[, mz := NA_real_] # make it present to simplify code below
    
    if (!is.null(adduct))
        suspects[is.na(mz), mz := calculateMasses(neutralMass, ..adduct, type = "mz")]
    else if (!is.null(suspects[["adduct"]]))
    {
        unAdducts <- sapply(unique(suspects[is.na(mz) & !is.na(adduct) & nzchar(adduct)]$adduct), as.adduct)
        suspects[is.na(mz) & !is.na(adduct) & nzchar(adduct), mz := calculateMasses(neutralMass, unAdducts[adduct], type = "mz")]
    }
    
    return(suspects[])
}

prepareSuspectList <- function(suspects, adduct, skipInvalid, checkDesc, prefCalcChemProps, neutralChemProps,
                               calcMZs = TRUE)
{
    hash <- makeHash(suspects, adduct, skipInvalid, checkDesc, prefCalcChemProps, neutralChemProps, calcMZs)
    cd <- loadCacheData("screenSuspectsPrepList", hash)
    if (!is.null(cd))
        suspects <- cd
    else
    {
        if (is.data.table(suspects))
            suspects <- copy(suspects)
        else
            suspects <- as.data.table(suspects)
        
        # convert to character in case factors are given...
        for (col in c("name", "formula", "SMILES", "InChI", "adduct"))
        {
            if (!is.null(suspects[[col]]))
                suspects[, (col) := as.character(get(col))]
        }
        
        # convert to numerics, may be logical if all are NA...
        for (col in c("mz", "rt"))
        {
            if (!is.null(suspects[[col]]))
                suspects[, (col) := as.numeric(get(col))]
        }
        
        # make name column file safe and unique
        sanNames <- fs::path_sanitize(suspects$name, replacement = "_")
        sanNames <- strtrim(sanNames, 150) # UNDONE: make max length configurable?
        sanNames <- make.unique(sanNames, sep = "-")
        changedNames <- which(sanNames != suspects$name)
        if (length(changedNames) > 0)
        {
            many <- length(changedNames) > 100
            if (many)
                changedNames <- changedNames[seq_len(100)]
            msg <- paste0("The following suspect names were changed to make them file compatible and/or unique:\n",
                          paste0(suspects$name[changedNames], " --> ", sanNames[changedNames], collapse = "\n"))
            if (many)
                msg <- paste0(msg, "\n", "(only the first 100 occurences are shown).")
            warning(msg, call. = FALSE)
            suspects[, name_orig := name]
            suspects[, name := sanNames]
        }
        
        if (checkDesc)
            suspects <- prepareChemTable(suspects, prefCalcChemProps = prefCalcChemProps,
                                         neutralChemProps = neutralChemProps)
        
        # calculate ionic masses if possible (not possible if no adducts are given and fGroups are annotated)
        if (calcMZs && (is.null(suspects[["mz"]]) || any(is.na(suspects[["mz"]]))) &&
            (!is.null(adduct) || !is.null(suspects[["adduct"]])))
            suspects <- calcSuspMZs(suspects, adduct)
        else if (is.null(suspects[["mz"]]))
            suspects[, mz := NA_real_]

        saveCacheData("screenSuspectsPrepList", suspects, hash)
    }        
    
    if (calcMZs)
    {
        # check for any suspects without proper mass info
        isNA <- is.na(suspects$neutralMass) & is.na(suspects$mz)
        if (any(isNA))
        {
            wrong <- paste0(sprintf("%s (line %d)", suspects$name[isNA], which(isNA)), collapse = "\n")
            if (skipInvalid)
            {
                warning(paste("Ignored following suspects for which no mass could be calculated:",
                              wrong), call. = FALSE)
                suspects <- suspects[!isNA]
            }
            else
                stop(paste("Could not calculate ion masses for the following suspects: ", wrong), call. = FALSE)
        }
    }
    
    return(suspects)
}

expandSuspMobilities <- function(suspects)
{
    hasMob <- !is.null(suspects[["mobility_susp"]]); hasCCS <- !is.null(suspects[["CCS_susp"]])
    if (!hasMob && !hasCCS)
        return(copy(suspects))

    verifyCol <- function(col)
    {
        # NOTE: also allow all NA columns, which are mostly logical by default
        if (!is.numeric(suspects[[col]]) && !is.character(suspects[[col]]) && !all(is.na(suspects[[col]])))
            stop(sprintf("%s column must be numeric or character (now %s)", col, class(suspects[[col]])), call. = FALSE)
    }
    
    if (hasMob)
        verifyCol("mobility_susp")
    if (hasCCS)
        verifyCol("CCS_susp")
    
    doSplit <- function(x)
    {
        if (is.na(x))
            return(NA_real_)
        if (is.character(x))
            return(if (nzchar(x)) as.numeric(unlist(strsplit(x, ";"))) else NA_real_)
        return(x)
    }
    
    suspNoMobs <- removeDTColumnsIfPresent(suspects, c("mobility", "CCS")) # for expanding below
    
    return(rbindlist(lapply(seq_len(nrow(suspects)), function(i)
    {
        mobs <- if (hasMob) doSplit(suspects$mobility_susp[i]) else NA_real_
        CCSs <- if (hasCCS) doSplit(suspects$CCS_susp[i]) else NA_real_
        if (length(mobs) != length(CCSs) && !checkmate::testScalarNA(mobs) && !checkmate::testScalarNA(CCSs))
            stop(sprintf("The length of mobility and CCS values for suspect '%s' (row %d) differs: %d/%d",
                         suspects$name[i], i, length(mobs), length(CCSs)), call. = FALSE)
        data.table(suspNoMobs[i], mobility = mobs, CCS = CCSs)
    })))
}

assignFeatureMobilitiesSuspects <- function(features, scr, IMSWindow, selectFunc = NULL)
{
    printf("Finding mobilities for all features from suspects...\n")
    oldCount <- countMobilityFeatures(features)
    
    if (is.null(scr[["mobility"]]))
    {
        warning("Cannot load mobilities from suspects: No suspect mobility data", call. = FALSE)
        return(features)
    }
    
    assignedMobilities <- copy(scr)
    assignedMobilities[, ngroup := .N, by = "group"]
    assignedMobilities <- assignedMobilities[ngroup == 1][, -"ngroup"]
    # UNDONE: message which were omitted
    assignedMobilities <- expandSuspMobilities(assignedMobilities)
    assignedMobilities <- assignedMobilities[!is.na(mobility)]
    assignedMobilities <- assignedMobilities[, .(group, mobility, mobmin = mobility - IMSWindow, mobmax = mobility + IMSWindow)]
    
    features@features <- Map(features@features, names(features@features), f = function(fTable, ana)
    {
        mobTable <- copy(assignedMobilities)
        mobTable <- mobTable[group %chin% fTable$group]
        if (!is.null(selectFunc))
            mobTable <- selectFunc(mobTable, ana)
        if (nrow(mobTable) == 0)
            mobTable[, ims_parent_ID := character()]
        else
            mobTable[, ims_parent_ID := fTable[match(mobTable$group, group, nomatch = 0)]$ID][, group := NULL]
        mobTable[, c("mob_area", "mob_intensity", "mob_assign_method") := .(NA_real_, NA_real_, "suspect")]
        return(doAssignFeatureMobilities(fTable, mobTable))
    })
    
    printf("Assigned %d mobility features.\n", countMobilityFeatures(features) - oldCount)
    
    return(features)
}

expandAndUpdateScreenInfoForIMS <- function(scr, gInfo)
{
    scr <- expandTableForIMSFGroups(scr, gInfo)
    scr[, d_rt := gInfo$ret[match(group, gInfo$group)] - rt] # update
    return(scr[])
}

finalizeScreenInfoForIMS <- function(scr, gInfo, IMSMatchParams)
{
    # shared code for assignMobilities() and screenSuspects()

    scrOrigExp <- expandSuspMobilities(scr)
    
    # pick mobility that is closest to that of a feature group
    scr[, mob_group := gInfo$mobility[match(group, gInfo$group)]]
    scr[, ims_parent_group := gInfo$ims_parent_group[match(group, gInfo$group)]]
    
    # assign closest mobility/CCS from suspect list
    scr[, c("mobility", "CCS") := NA_real_]
    scr[!is.na(mob_group), c("mobility", "CCS") := {
        g <- group; n <- name
        soe <- scrOrigExp[group == g & name == n]
        wh <- which.min(abs(soe$mobility - mob_group))
        if (length(wh) == 0)
            list(NA_real_, NA_real_)
        else
            list(soe$mobility[wh], soe$CCS[wh])
    }, by = seq_len(nrow(scr[!is.na(mob_group)]))]

    scr <- assignTabIMSDeviations(scr, gInfo)

    if (!is.null(IMSMatchParams))
    {
        checkCol <- if (IMSMatchParams$param == "mobility") "d_mob" else "d_CCS"
        if (IMSMatchParams$relative)
            checkCol <- paste0(checkCol, "_rel")
        
        scr <- scr[is.na(get(checkCol)) | abs(get(checkCol)) <= IMSMatchParams$window]
        
        if (IMSMatchParams$minMatches > 0)
        {
            scr[, keep := (.N - 1) >= min(IMSMatchParams$minMatches, scrOrigExp[ims_parent_group == group, .N]),
                by = c("ims_parent_group", "name")]
            scr <- scr[keep == TRUE, -"keep"]
        }
    }
    
    return(removeDTColumnsIfPresent(scr, c("mob_group", "ims_parent_group")))
}

selectFromSuspAdductCol <- function(tab, col, fgAnn, adductChrDef)
{
    isUsable <- function(v) !is.na(v) & (!is.character(v) | nzchar(v))
    
    # preference order: regular non-adduct col > col from (adduct arg > adduct column > feat annotations)
    
    # special case: we can use all data from regular non-adduct column
    if (!is.null(tab[[col]]) && all(isUsable(tab[[col]])))
        return(tab[[col]])
    
    tryAddCol <- function(add, row)
    {
        cl <- if (!is.null(add)) paste0(col, "_", add) else col
        if (!is.null(tab[[cl]]) && isUsable(tab[[cl]][row]))
            return(tab[[cl]][row])
        return(NA_real_)
    }
    
    # UNDONE: would be nice to vectorize this somehow...
    return(sapply(seq_len(nrow(tab)), function(i)
    {
        v <- NA_real_
        
        if (!is.null(adductChrDef))
            v <- tryAddCol(adductChrDef, i)
        if (is.na(v) && !is.null(tab[["adduct"]]) && !is.na(tab$adduct[i]) && nzchar(tab$adduct[i]))
            v <- tryAddCol(tab$adduct[i], i)
        if (is.na(v) && !is.null(tab[["group"]]) && nrow(fgAnn) > 0)
            v <- tryAddCol(fgAnn[group == tab$group[i]]$adduct, i)
        
        return(v)
    }))
}

doScreenSuspects <- function(fGroups, suspects, rtWindow, mzWindow, IMSMatchParams, adduct, skipInvalid)
{
    gInfo <- groupInfo(fGroups)
    annTbl <- annotations(fGroups)
    
    # NOTE: rt is always included
    metaDataCols <- union("rt", intersect(suspMetaDataCols(), names(suspects)))
    
    # HACK: the mobility columns are handled differently
    metaDataCols <- setdiff(metaDataCols, c("mobility", "mobility_susp", "CCS", "CCS_mobility"))
    
    emptyResult <- data.table()
    for (col in c(metaDataCols, "group", "d_rt", "d_mz"))
    {
        if (col %in% c("rt", "mz", "neutralMass", "d_rt", "d_mz"))
            emptyResult[, (col) := numeric()]
        else
            emptyResult[, (col) := character()]
    }
    
    setMetaData <- function(t, suspRow)
    {
        for (col in metaDataCols)
        {
            if (col == "rt" && is.null(suspects[["rt"]]))
                set(t, 1L, col, NA_real_) # exception: always want this column
            else
                set(t, 1L, col, suspRow[[col]])
        }
        return(t)
    }        
    
    if (length(fGroups) > 0)
    {
        adductTxt <- if (!is.null(adduct)) as.character(adduct)
            
        prog <- openProgBar(0, nrow(suspects))
        
        retlist <- lapply(seq_len(nrow(suspects)), function(ti)
        {
            hasRT <- !is.null(suspects[["rt"]]) && !is.na(suspects$rt[ti])
            
            gi <- gInfo
            
            # only consider nearby eluting fGroups if RTs are available
            if (hasRT)
                gi <- gInfo[numLTE(abs(ret - suspects$rt[ti]), rtWindow)]
            
            # match by mz, this is either done by...
            #   - fGroup ionized mass and 'mz' column from suspect data if the latter is available
            #   - fGroup ionized mass and calculated suspect ionized mass, only if adducts were specified
            #     (mandatory if no adduct annotations available). Note that ionized masses are already calculated by
            #     prepareSuspectList() and stored in the mz column.
            #   - Neutralized fGroup/suspect mass (only if adduct annotations are available, this case mz column is NA)
            
            if (is.na(suspects$mz[ti])) # no ionized suspect available, must use annotation data to compare neutral masses
            {
                at <- annTbl[group %chin% gi$group & numLTE(abs(neutralMass - suspects[ti]$neutralMass), mzWindow)]
                gi <- gi[group %chin% at$group]
            }
            else
                gi <- gi[numLTE(abs(gi$mz - suspects$mz[ti]), mzWindow)]
            
            if (nrow(gi) == 0)
                hits <- copy(emptyResult) # no hits
            else
            {
                hits <- rbindlist(Map(gi$group, gi$ret, gi$mz, f = function(g, gret, gmz)
                {
                    ret <- data.table()
                    setMetaData(ret, suspects[ti])

                    # copy the right mobility and CCS columns from the suspect list
                    ret[, mobility_susp := selectFromSuspAdductCol(suspects[ti], "mobility", annTbl, adductTxt)]
                    ret[, CCS_susp := selectFromSuspAdductCol(suspects[ti], "CCS", annTbl, adductTxt)]
                    
                    ret[, c("group", "d_rt", "d_mz") := .(g, d_rt = if (hasRT) gret - rt else NA_real_,
                                                          ifelse(is.na(mz), annTbl[group == g]$neutralMass - neutralMass,
                                                                 gmz - mz))]
                    
                    return(ret)
                }), fill = TRUE)
            }
            
            setTxtProgressBar(prog, ti)
            return(hits)
        })
        ret <- rbindlist(retlist, fill = TRUE)
        setcolorder(ret, "name")
        
        if (hasMobilities(fGroups))
            ret <- finalizeScreenInfoForIMS(ret, gInfo, IMSMatchParams)
        
        setTxtProgressBar(prog, nrow(suspects))
        close(prog)
    }
    else
        ret <- copy(emptyResult)
    
    suspectsn <- nrow(suspects)
    foundn <- uniqueN(ret$name)
    printf("Found %d/%d suspects (%.2f%%)\n", foundn, suspectsn, foundn * 100 / suspectsn)
    
    return(ret[])
}

# method definition for screening methods of screenSuspects()
doScreenSuspectsAmend <- function(fGroups, suspects, rtWindow, mzWindow, IMSMatchParams, adduct, skipInvalid,
                                  prefCalcChemProps, neutralChemProps, onlyHits, amend = FALSE)
{
    aapply(checkmate::assertFlag, . ~ onlyHits + amend)
    
    fGroupsScreened <- callNextMethod(fGroups, suspects, rtWindow, mzWindow, IMSMatchParams, adduct, skipInvalid,
                                      prefCalcChemProps, neutralChemProps, onlyHits)
    if (!amend)
        return(fGroupsScreened)
    
    # amend screening results
    
    fGroups@screenInfo <- rbind(fGroups@screenInfo, fGroupsScreened@screenInfo, fill = TRUE)
    fGroups@screenInfo <- unique(fGroups@screenInfo, by = c("name", "group"))
    
    if (onlyHits)
        fGroups <- fGroups[, fGroups@screenInfo$group]
    
    return(fGroups)
}

delScreening <- function(fGroups, j, k)
{
    scr <- copy(screenInfo(fGroups))
    
    ac <- checkmate::makeAssertCollection()
    if (!is.function(j))
        j <- assertDeleteArgAndToChr(j, names(fGroups), add = ac)
    checkmate::assert(
        checkmate::checkFunction(k),
        checkmate::checkChoice(k, unique(scr$name)),
        checkmate::checkScalarNA(k),
        .var.name = "k", add = ac
    )
    checkmate::reportAssertions(ac)
    
    if (is.null(j))
        j <- names(fGroups)
    
    scr[, keep := TRUE]
    
    if (is.atomic(k) && is.na(k))
        scr[group %chin% j, keep := FALSE]
    else if (is.function(k))
        scr[group %chin% j, keep := !k(.SD)]
    else
        scr[group %chin% j, keep := !name %chin% k]
    
    scr <- scr[keep == TRUE][, keep := NULL][]
    fGroups@screenInfo <- scr
    return(fGroups)
}

doSFGroupsScreeningDelete <- function(obj, i = NULL, j = NULL, k = NULL, ...)
{
    if (!is.null(k))
    {
        if (!is.null(i))
            stop("Cannot specify i and k arguments simultaneously.", call. = FALSE)
        return(delScreening(obj, j, k))
    }
    
    # NOTE: this function is a method definition, so we can call callNextMethod()
    obj <- callNextMethod()
    obj@screenInfo <- obj@screenInfo[group %chin% names(obj)]
    return(obj)
}

# used by groupFeatures methods for featuresSuspects
doGroupSuspects <- function(feat, groupFunc, ..., verbose = TRUE)
{
    anaInfo <- analysisInfo(feat)
    fTable <- featureTable(feat)
    detectedSusps <- unique(unlist(lapply(fTable, "[[", "suspect")))
    hasMob <- hasMobilities(feat)
    fgSusps <- sapply(detectedSusps, function(susp)
    {
        featSub <- delete(feat, j = function(ft, ...) ft$suspect != susp)
        ret <- groupFunc(featSub, ..., verbose = FALSE)
        if (hasMob)
        {
            ft <- as.data.table(getFeatures(ret))
            ret@groupInfo[, mobility := mean(ft$mobility[ft$group == group]), by = "group"]
            ret@groupInfo[, ims_parent_group := NA_character_]
        }
        return(ret)
    }, simplify = FALSE)
    
    fgInfoAll <- rbindlist(lapply(fgSusps, function(fg) copy(groupInfo(fg))), idcol = "suspect")
    fgInfoAll[, group_susp := {
        if (.N > 1)
            paste0(suspect, "-", seq_len(.N))
        else
            suspect
    }, by = "suspect"]

    gTable <- data.table()
    gTable[, (fgInfoAll$group_susp) := Map(fgInfoAll$suspect, fgInfoAll$group, f = function(s, g)
    {
        ints <- numeric(length(fTable))
        fg <- fgSusps[[s]][, g]
        ints[match(analyses(fg), anaInfo$analysis)] <- fg[[1]]
        return(ints)
    })]
    
    ftind <- data.table()
    ftind[, (fgInfoAll$group_susp) := Map(fgInfoAll$suspect, fgInfoAll$group, f = function(s, g)
    {
        inds <- integer(length(fTable))
        fg <- removeEmptyAnalyses(fgSusps[[s]][, g])
        anaInds <- match(analyses(fg), anaInfo$analysis)
        inds[anaInds] <- mapply(featureTable(fg), anaInds, FUN = function(ft, ai) match(ft$ID, fTable[[ai]]$ID))
        return(inds)
    })]
    
    gInfo <- subsetDTColumnsIfPresent(copy(fgInfoAll), c("group_susp", "ret", "mz", "mobility", "ims_parent_group"))
    setnames(gInfo, "group_susp", "group")

    # set screenInfo
    sInfo <- fgInfoAll[, c("suspect", "group_susp"), with = FALSE]
    setnames(sInfo, c("suspect", "group_susp"),  c("name", "group"))
    
    metaDataCols <- union("rt", intersect(suspMetaDataCols(), names(feat@suspects)))
    susp <- copy(feat@suspects)
    if (is.null(susp[["rt"]]))
        susp[, rt := NA_real_]
    sInfo <- merge(sInfo, susp[, metaDataCols, with = FALSE], by = "name")
    sInfo[, d_rt := gInfo[match(sInfo$group, group)]$ret - rt]
    sInfo[, d_mz := gInfo[match(sInfo$group, group)]$mz - mz]
    if (hasMob)
        sInfo[, d_mob := gInfo[match(sInfo$group, group)]$mobility - mobility]
    return(featureGroupsScreening(screenInfo = sInfo, groups = gTable, groupInfo = gInfo, features = feat,
                                  ftindex = ftind))
}

doSuspectFilter <- function(obj, onlyHits, selectHitsBy, selectBestFGroups, maxLevel, maxFormRank, maxCompRank,
                            minAnnSimForm, minAnnSimComp, minAnnSimBoth, absMinFragMatches, relMinFragMatches, minRF,
                            maxLC50, negate, applyIMS)
{
    if (nrow(screenInfo(obj)) > 0)
    {
        filteredSI <- copy(screenInfo(obj))
        filteredSI[, rowID := seq_len(nrow(filteredSI))] # so we can see what was removed (for applyIMS)
        
        colFilter <- function(si, pred, what, col, dataWhich, funcToRun, ac = TRUE)
        {
            val <- get(what)
            if (!is.null(val))
            {
                allCols <- if (ac)
                    getAllSuspCols(col, names(screenInfo(obj)), mergedConsensusNames(obj))
                else
                    intersect(col, names(screenInfo(obj)))
                
                if (length(allCols) == 0)
                    warning(sprintf("Cannot apply %s filter: no %s data available (did you run %s()?).", what,
                                    dataWhich, funcToRun), call. = FALSE)
                else
                {
                    if (negate)
                        doPred <- function(x, v) is.na(x) | !nzchar(x) | !pred(x, v)
                    else
                        doPred <- function(x, v) !is.na(x) & nzchar(x) & pred(x, v)
                    
                    # ensure at least one column follows predicate
                    si <- si[rowSums(sapply(mget(allCols), doPred, val)) >= 1]
                }
            }
            return(si)
        }
        colFilterAnn <- function(...) colFilter(dataWhich = "annotation", funcToRun = "annotateSuspects", ...)
        minPred <- function(x, v) x >= v
        maxPred <- function(x, v) x <= v
        levPred <- function(x, v) maxPred(numericIDLevel(x), v)
        
        filteredSI <- colFilterAnn(filteredSI, levPred, "maxLevel", "estIDLevel", ac = FALSE)
        filteredSI <- colFilterAnn(filteredSI, maxPred, "maxFormRank", "formRank", ac = FALSE)
        filteredSI <- colFilterAnn(filteredSI, maxPred, "maxCompRank", "compRank", ac = FALSE)
        filteredSI <- colFilterAnn(filteredSI, minPred, "minAnnSimForm", "annSimForm")
        filteredSI <- colFilterAnn(filteredSI, minPred, "minAnnSimComp", "annSimComp")
        filteredSI <- colFilterAnn(filteredSI, minPred, "minAnnSimBoth", "annSimBoth")
        filteredSI <- colFilterAnn(filteredSI, minPred, "absMinFragMatches", "maxFragMatches")
        filteredSI <- colFilterAnn(filteredSI, minPred, "relMinFragMatches", "maxFragMatchesRel")
        
        filteredSI <- colFilter(filteredSI, minPred, "minRF", "RF_SMILES", dataWhich = "response factor", funcToRun = "predictRespFactors")
        filteredSI <- colFilter(filteredSI, maxPred, "maxLC50", "LC50_SMILES", dataWhich = "LC50", funcToRun = "predictTox")
        
        # do here so that only duplicates not yet filtered out in previous steps are considered
        # NOTE for sets: for ID levels only the regular (non-set) estIDLevel column is used
        if (!is.null(selectHitsBy) || selectBestFGroups)
        {
            doKeep <- function(v, d) is.na(v) | length(v) == 1 | seq_along(v) == order(v, decreasing = d)[1]
            doSelectFilter <- function(si, by, byCol)
            {
                if (by == "level" && is.null(si[["estIDLevel"]]))
                    warning("Cannot select by identification level: no annotation data available (did you run annotateSuspects()?).")
                else
                {
                    gTab <- as.data.table(obj, collapseSuspects = NULL, onlyHits = TRUE)
                    # equalize names with screenInfo
                    if (!is.null(gTab[["adduct"]]))
                    {
                        # may be there if adduct annotations are available, remove to not interfere with susp_adduct
                        gTab[, adduct := NULL]
                    }
                    suspnames <- grep("^susp_", names(gTab), value = TRUE)
                    setnames(gTab, suspnames, sub("^susp_", "", suspnames))
                    
                    if (by == "intensity")
                    {
                        gTab[, avgInts := rowMeans(.SD), .SDcol = getADTIntCols(analyses(obj))]
                        gTab <- gTab[, keep := doKeep(avgInts, !negate), by = byCol]
                    }
                    else # select by best hit
                        gTab <- gTab[, keep := doKeep(estIDLevel, negate), by = byCol]
                    
                    if (any(!gTab$keep))
                    {
                        # merge-in keep column so we can subset screenInfo
                        si <- copy(si)
                        si[gTab, keep := i.keep, on = c("group", "name")]
                        setorderv(si, "name")
                        si <- si[keep == TRUE, -"keep"]
                    }
                }
                return(si)
            }
            
            if (!is.null(selectHitsBy))
                filteredSI <- doSelectFilter(filteredSI, selectHitsBy, "name")
            if (selectBestFGroups)
                filteredSI <- doSelectFilter(filteredSI, "level", "group")
        }
        
        if (applyIMS != "both" && hasMobilities(obj))
        {
            origSI <- screenInfo(obj)
            rowsRM <- setdiff(seq_len(nrow(origSI)), filteredSI$rowID)
            groupsKeep <- if (applyIMS) groupInfo(obj)[is.na(mobility)]$group else groupInfo(obj)[!is.na(mobility)]$group
            rowsRM <- rowsRM[!origSI[rowsRM]$group %in% groupsKeep]
            if (length(rowsRM) > 0)
                obj@screenInfo <- obj@screenInfo[-rowsRM]
        }
        else
            obj@screenInfo <- filteredSI[, -"rowID"]
    }
    
    # NOTE: do last in case previous steps removed hits 
    if (!is.null(onlyHits))
    {
        obj <- doFGroupsFilter(obj, "suspects", list(), function(fGroups)
        {
            sGroups <- unique(screenInfo(fGroups)$group)
            if (negate && onlyHits)
                fGroups <- fGroups[, setdiff(names(fGroups), sGroups)]
            else
                fGroups <- fGroups[, sGroups]
            return(fGroups)
        }, applyIMS = applyIMS)
    }
    
    return(obj)
}

annotatedMSMSSimilarity <- function(annPL, specSimParams)
{
    if (is.null(annPL)) # mainly to handle formula candidates w/out MS/MS
        return(0)
    
    annPL <- prepSpecSimilarityPL(annPL, specSimParams$removePrecursor, specSimParams$relMinIntensity,
                                  specSimParams$minPeaks)
    
    if (nrow(annPL) == 0 || !any(annPL$annotated))
        return(0)
    
    annPLAnn <- annPL[annotated == TRUE]
    
    return(drop(specDistRect(list(annPLAnn), list(annPL), specSimParams$method, "none", 0, 0,
                             specSimParams$mzWeight, specSimParams$intWeight, specSimParams$absMzDev)))
}

# used as predictTox method definitions
doPredictToxSuspects <-  function(obj, LC50Mode = "static", concUnit = "ugL")
{
    checkPackage("MS2Tox", "kruvelab/MS2Tox")
    
    checkmate::assertChoice(LC50Mode, c("static", "flow"))
    assertConcUnit(concUnit)
    
    if (length(obj) == 0)
        return(obj)
    
    scr <- screenInfo(obj)
    if (is.null(scr[["SMILES"]]) || all(is.na(scr$SMILES)))
        stop("Suspects lack necessary SMILES information to perform calculations, aborting...", call. = FALSE)
    if (any(is.na(scr$SMILES)))
        warning("Some suspect SMILES are NA and will be ignored", call. = FALSE)
    
    # avoid duplicate calculations if there happen to be suspects with the same SMILES
    inpSMILES <- unique(screenInfo(obj)$SMILES)
    inpSMILES <- inpSMILES[!is.na(inpSMILES)]
    
    printf("Predicting LC50 values from SMILES with MS2Tox for %d suspects...\n", length(inpSMILES))
    pr <- predictLC50SMILES(inpSMILES, LC50Mode, concUnit)
    
    if (!is.null(scr[["LC50_SMILES"]]))
        scr[, LC50_SMILES := NULL] # clearout for merge below
    scr <- merge(scr, pr, by = "SMILES", sort = FALSE, all.x = TRUE)
    
    obj@screenInfo <- scr
    
    return(obj)
}
