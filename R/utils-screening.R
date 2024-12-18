# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include utils.R
#' @include utils-compounds.R
NULL

isScreening <- function(fGroups) inherits(fGroups, c("featureGroupsScreening", "featureGroupsScreeningSet"))
isSuspAnnotated <- function(fGroups) isScreening(fGroups) && !is.null(screenInfo(fGroups)[["estIDLevel"]])

suspMetaDataCols <- function() c("name", "rt", "name_orig", "mz", "SMILES", "InChI", "InChIKey", "formula",
                                 "neutralMass", "molNeutralized", "adduct", "fragments_mz", "fragments_formula")
suspAnnCols <- function() c("formRank", "compRank", "annSimForm", "annSimComp", "annSimBoth", "maxFrags",
                            "maxFragMatches", "maxFragMatchesRel", "estIDLevel")

getAllSuspCols <- function(targetCols, allCols, mConsNames)
{
    targetCols <- c(targetCols, sapply(targetCols, function(cl) paste0(cl, "-", mConsNames),
                                       USE.NAMES = FALSE))
    return(intersect(targetCols, allCols))
}

doScreeningShow <- function(obj)
{
    printf("Suspects: %s (%d hits total)\n", getStrListWithMax(unique(screenInfo(obj)$name), 6, ", "),
           nrow(screenInfo(obj)))
    printf("Suspects annotated: %s\n", if (isSuspAnnotated(obj)) "yes" else "no")
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
        {
            if (is.null(suspects[["mz"]]))
                suspects[, mz := NA_real_] # make it present to simplify code below
            
            if (!is.null(adduct))
                suspects[is.na(mz), mz := calculateMasses(neutralMass, ..adduct, type = "mz")]
            else
            {
                unAdducts <- sapply(unique(suspects[is.na(mz)]$adduct), as.adduct)
                suspects[is.na(mz) & !is.na(adduct), mz := calculateMasses(neutralMass, unAdducts[adduct], type = "mz")]
            }
        }
        else if (is.null(suspects[["mz"]]))
        {
            # NOTE: if mz column is already available it either contains user values or already NAs
            suspects[, mz := NA_real_]
        }
        
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
                              wrong))
                suspects <- suspects[!isNA]
            }
            else
                stop(paste("Could not calculate ion masses for the following suspects: "), wrong)
        }
    }
    
    return(suspects)
}

doScreenSuspects <- function(fGroups, suspects, rtWindow, mzWindow, skipInvalid)
{
    gInfo <- groupInfo(fGroups)
    annTbl <- annotations(fGroups)
    
    # NOTE: rt is always included
    metaDataCols <- union("rt", intersect(suspMetaDataCols(), names(suspects)))
    
    emptyResult <- function()
    {
        ret <- data.table()
        for (col in c(metaDataCols, "group", "d_rt", "d_mz"))
        {
            if (col %in% c("rt", "mz", "neutralMass", "d_rt", "d_mz"))
                ret[, (col) := numeric()]
            else
                ret[, (col) := character()]
        }
        return(ret)
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
        prog <- openProgBar(0, nrow(suspects))
        
        retlist <- lapply(seq_len(nrow(suspects)), function(ti)
        {
            hasRT <- !is.null(suspects[["rt"]]) && !is.na(suspects$rt[ti])
            
            gi <- gInfo
            
            # only consider nearby eluting fGroups if RTs are available
            if (hasRT)
                gi <- gInfo[numLTE(abs(gInfo$rts - suspects$rt[ti]), rtWindow), ]
            
            # match by mz, this is either done by...
            #   - fGroup ionized mass and 'mz' column from suspect data if the latter is available
            #   - fGroup ionized mass and calculated suspect ionized mass, only if adducts were specified
            #     (mandatory if no adduct annotations available). Note that ionized masses are already calculated by
            #     prepareSuspectList() and stored in the mz column.
            #   - Neutralized fGroup/suspect mass (only if adduct annotations are available, this case mz column is NA)
            
            if (is.na(suspects$mz[ti])) # no ionized suspect available, must use annotation data to compare neutral masses
            {
                at <- annTbl[group %in% rownames(gi) & numLTE(abs(neutralMass - suspects[ti]$neutralMass), mzWindow)]
                gi <- gi[at$group, ]
            }
            else
                gi <- gi[numLTE(abs(gi$mzs - suspects$mz[ti]), mzWindow), ]
            
            if (nrow(gi) == 0)
                hits <- emptyResult() # no hits
            else
            {
                hits <- rbindlist(lapply(rownames(gi), function(g)
                {
                    ret <- data.table()
                    setMetaData(ret, suspects[ti])
                    ret[, c("group", "d_rt", "d_mz") := .(g, d_rt = if (hasRT) gi[g, "rts"] - rt else NA_real_,
                                                          ifelse(is.na(mz), annTbl[group == g]$neutralMass - neutralMass,
                                                                 gi[g, "mzs"] - mz))]
                    return(ret)
                }), fill = TRUE)
            }
            
            setTxtProgressBar(prog, ti)
            return(hits)
        })
        ret <- rbindlist(retlist, fill = TRUE)
        setcolorder(ret, "name")
        
        setTxtProgressBar(prog, nrow(suspects))
        close(prog)
    }
    else
        ret <- emptyResult()
    
    suspectsn <- nrow(suspects)
    foundn <- uniqueN(ret$name)
    printf("Found %d/%d suspects (%.2f%%)\n", foundn, suspectsn, foundn * 100 / suspectsn)
    
    return(ret[])
}

doSuspectFilter <- function(obj, onlyHits, selectHitsBy, selectBestFGroups, maxLevel, maxFormRank, maxCompRank,
                            minAnnSimForm, minAnnSimComp, minAnnSimBoth, absMinFragMatches, relMinFragMatches, minRF,
                            maxLC50, negate)
{
    if (nrow(screenInfo(obj)) > 0)
    {
        colFilter <- function(pred, what, col, dataWhich, funcToRun, ac = TRUE)
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
                    obj@screenInfo <- screenInfo(obj)[rowSums(sapply(mget(allCols), doPred, val)) >= 1]
                }
            }
            return(obj)
        }
        colFilterAnn <- function(...) colFilter(dataWhich = "annotation", funcToRun = "annotateSuspects", ...)
        minPred <- function(x, v) x >= v
        maxPred <- function(x, v) x <= v
        levPred <- function(x, v) maxPred(numericIDLevel(x), v)
        
        obj <- colFilterAnn(levPred, "maxLevel", "estIDLevel", ac = FALSE)
        obj <- colFilterAnn(maxPred, "maxFormRank", "formRank", ac = FALSE)
        obj <- colFilterAnn(maxPred, "maxCompRank", "compRank", ac = FALSE)
        obj <- colFilterAnn(minPred, "minAnnSimForm", "annSimForm")
        obj <- colFilterAnn(minPred, "minAnnSimComp", "annSimComp")
        obj <- colFilterAnn(minPred, "minAnnSimBoth", "annSimBoth")
        obj <- colFilterAnn(minPred, "absMinFragMatches", "maxFragMatches")
        obj <- colFilterAnn(minPred, "relMinFragMatches", "maxFragMatchesRel")
        
        obj <- colFilter(minPred, "minRF", "RF_SMILES", dataWhich = "response factor", funcToRun = "predictRespFactors")
        obj <- colFilter(maxPred, "maxLC50", "LC50_SMILES", dataWhich = "LC50", funcToRun = "predictTox")
        
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
                        gTab[, avgInts := rowMeans(.SD), .SDcol = analyses(obj)]
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
                        obj@screenInfo <- si[keep == TRUE, -"keep"]
                    }
                }
                return(obj@screenInfo)
            }
            
            if (!is.null(selectHitsBy))
                obj@screenInfo <- doSelectFilter(obj@screenInfo, selectHitsBy, "name")
            if (selectBestFGroups)
                obj@screenInfo <- doSelectFilter(obj@screenInfo, "level", "group")
        }
    }
    
    # NOTE: do last in case previous steps removed hits 
    if (!is.null(onlyHits))
    {
        sGroups <- unique(screenInfo(obj)$group)
        if (negate && onlyHits)
            obj <- obj[, setdiff(names(obj), sGroups)]
        else
            obj <- obj[, sGroups]
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

estimateIdentificationLevel <- function(suspectName, suspectFGroup, suspectRTDev, suspectInChIKey1, suspectFormula,
                                        suspectAnnSimForm, suspectAnnSimComp, suspectAnnSimBoth,
                                        maxSuspFrags, maxFragMatches, formTable, formRank, mFormNames, formScoreRanges,
                                        formulasNormalizeScores, compTable, compRank, mCompNames, compScoreRanges,
                                        compoundsNormalizeScores, absMzDev, IDLevelRules, logPath)
{
    formScores <- formScoreNames(FALSE); formNormScores <- formScoreNames(TRUE)
    compScores <- compScoreNames(FALSE); compNormScores <- compScoreNames(TRUE)
    fRow <- cRow <- NULL
    
    if (!is.null(formTable) && !is.null(suspectFormula))
    {
        formTableNorm <- normalizeAnnScores(formTable, formNormScores, formScoreRanges, mFormNames,
                                            formulasNormalizeScores == "minmax")
        if (!is.na(formRank))
        {
            fRow <- formTable[formRank]
            fRowNorm <- formTableNorm[formRank]
        }
    }
    
    if (!is.null(compTable) && !is.null(suspectInChIKey1))
    {
        compTableNorm <- normalizeAnnScores(compTable, compNormScores, compScoreRanges, mCompNames, compoundsNormalizeScores == "minmax")
        if (!is.na(compRank))
        {
            cRow <- compTable[compRank]
            cRowNorm <- compTableNorm[compRank]
        }
    }
    
    getValType <- function(val, IDType)
    {
        if (!is.list(val) || is.null(val[["type"]]) || !val[["type"]] %in% c("formula", "compound"))
            stop(sprintf("Need to specify the type (formula/compound) for %s!", IDType))
        return(val[["type"]])
    }
    
    getVal <- function(val, IDType, valType)
    {
        if (is.list(val))
        {
            if (is.null(val[[valType]]))
                stop(sprintf("Need to specify %s for %s!", valType, IDType))
            return(val[[valType]])
        }
        return(val)
    }
    
    getOptVal <- function(val, valType, default)
    {
        if (is.list(val) && !is.null(val[[valType]]))
            return(val[[valType]])
        return(default)
    }
    
    checkAnnotationScore <- function(val, scType, rank, annRow, annTable, annRowNorm, annTableNorm, scCols)
    {
        scCols <- scCols[!is.na(unlist(annRow[, scCols, with = FALSE]))]
        if (length(scCols) == 0)
            return("score not available")
        
        minValue <- getVal(val, scType, "min")
        
        if (getOptVal(val, "relative", FALSE))
        {
            annRow <- annRowNorm
            annTable <- annTableNorm
        }
        
        scoreVal <- rowMeans(annRow[, scCols, with = FALSE], na.rm = TRUE)
        if (is.infinite(scoreVal)) # only NA values
            return("no value(s) available for this score")
        if (scoreVal < minValue)
            return(sprintf("(average) score too low: %f/%f", scoreVal, minValue))
        
        htn <- getOptVal(val, "higherThanNext", 0)
        if (htn > 0 && nrow(annTable) > 1)
        {
            otherVals <- rowMeans(annTable[-rank, scCols, with = FALSE], na.rm = TRUE)
            otherVals <- otherVals[is.finite(otherVals)] # remove results for NA rows
            if (length(otherVals) > 0)
            {
                otherHighest <- max(otherVals)
                if (is.infinite(htn)) # special case: should be highest
                {
                    if (otherHighest > 0)
                        return("not the highest score")
                }
                else if ((scoreVal - otherHighest) < htn)
                    return(sprintf("difference with highest score from other candidates too low: %f/%f", scoreVal - otherHighest, htn))
            }
        }
        
        return(TRUE)            
    }

    if (!is.null(logPath))
    {
        logDir <- R.utils::getAbsolutePath(logPath) # take absolute path for length calculation below
        logFile <- paste0(suspectName, "-", suspectFGroup, ".txt")
        
        # check if path would be too long for e.g Windows systems, which may happen with very long suspect names
        pathLen <- nchar(logDir) + nchar(logFile) + 1 # +1: path separator
        if (pathLen > 255)
        {
            # truncate end part of suspect name
            logFile <- paste0(substr(suspectName, 1, nchar(suspectName) - (pathLen - 255)), "-", suspectFGroup, ".txt")
        }
        logOut <- file.path(logDir, logFile)
        
        mkdirp(dirname(logOut))
        logFile <- withr::local_connection(file(logOut, "w"))
        doLog <- function(indent, s, ...) fprintf(logFile, paste0(strrep(" ", indent * 4), s), ...)
    }
    else
        doLog <- function(...) NULL
    
    formCompScores <- intersect(formScores, compScores)
    allScores <- union(formScores, compScores)
    
    checkLevelOK <- function(IDL, indent = 0)
    {
        indent <- indent + 1
        for (type in names(IDL))
        {
            levelOK <- NULL
            levelFailed <- NULL
            val <- IDL[[type]]
            
            if (type %in% c("rank", "annMSMSSim", formCompScores))
                doLog(indent, "Checking ID level type '%s' (for %s)\n", type,
                      getValType(val, type))
            else
                doLog(indent, "Checking ID level type '%s'\n", type)
            
            if (type %in% c("or", "and"))
            {
                if (!is.list(val) || checkmate::testNamed(val))
                    stop(sprintf("Specify a list with '%s'", type))
                checkf <- if (type == "or") any else all
                levelOK <- checkf(mapply(val, seq_along(val), FUN = function(IDL, i)
                {
                    doLog(indent + 1, "check %s condition %d/%d\n", toupper(type), i, length(val))
                    return(checkLevelOK(IDL, indent + 2))
                }))
            }
            else if (type == "all" && val == TRUE)
                levelOK <- TRUE # special case: this level is always valid
            else if (type == "suspectFragments")
            {
                minFrags <- min(val, maxSuspFrags, na.rm = TRUE)
                levelOK <- !is.na(maxFragMatches) && maxFragMatches >= minFrags
                if (!levelOK)
                {
                    if (is.na(maxFragMatches))
                        levelFailed <- "no fragments to match"
                    else
                        levelFailed <- sprintf("not enough fragments: %d/%d", maxFragMatches, minFrags)
                }
            }
            else if (type == "retention")
            {
                rtm <- getVal(val, type, "max")
                levelOK <- !is.null(suspectRTDev) && !is.na(suspectRTDev) && numLTE(abs(suspectRTDev), rtm)
                if (!levelOK)
                {
                    if (is.null(suspectRTDev) && is.na(suspectRTDev))
                        levelFailed <- "no retention time information available"
                    else
                        levelFailed <- sprintf("too high retention time deviation: %f/%f",
                                               abs(suspectRTDev), rtm)
                }
            }
            else if (type == "rank")
            {
                r <- if (getValType(val, type) == "formula") formRank else compRank
                maxR <- getVal(val, type, "max")
                levelOK <- !is.na(r) && r <= maxR
                if (!levelOK)
                    levelFailed <- if (is.na(r)) "candidate not ranked" else sprintf("ranked too low: %d/%d", r, maxR)
            }
            else if (type == "annMSMSSim")
            {
                sim <- if (getValType(val, type) == "formula") suspectAnnSimForm else suspectAnnSimComp
                minSim <- getVal(val, type, "min")
                levelOK <- !is.na(sim) && numGTE(sim, minSim)
                if (!levelOK)
                    levelFailed <- if (is.na(sim)) "no calculated similarity" else sprintf("similarity too low: %f/%f", sim, minSim)
            }
            else if (type == "annMSMSSimBoth")
            {
                minSim <- getVal(val, type, "min")
                levelOK <- !is.na(suspectAnnSimBoth) && numGTE(suspectAnnSimBoth, minSim)
                if (!levelOK)
                    levelFailed <- if (is.na(suspectAnnSimBoth)) "no calculated similarity" else
                        sprintf("similarity too low: %f/%f", suspectAnnSimBoth, minSim)
            }
            else if (type %in% allScores)
            {
                if (type %in% formCompScores)
                    isForm <- getValType(val, type) == "formula"
                else
                    isForm <- type %in% formScores
                
                if (isForm)
                    levelOK <- checkAnnotationScore(val, type, formRank, fRow, formTable, fRowNorm, formTableNorm,
                                                    getAllMergedConsCols(type, names(formTable), mFormNames))
                else
                    levelOK <- checkAnnotationScore(val, type, compRank, cRow, compTable, cRowNorm, compTableNorm,
                                                    getAllMergedConsCols(type, names(compTable), mCompNames))
                
                if (!isTRUE(levelOK))
                {
                    levelFailed <- levelOK
                    levelOK <- FALSE
                }
            }
            else
                stop(paste("Unknown ID level type:", type))
            
            if (!levelOK)
            {
                doLog(indent, "ID level failed: %s\n", levelFailed)
                return(FALSE)
            }
            doLog(indent, "ID level type passed!\n")
        }
        
        return(TRUE)
    }
    
    doLog(0, "Estimating identification level for '%s' to feature group '%s'\n---\n", suspectName, suspectFGroup)
    
    if (is.null(fRow))
    {
        if (is.null(suspectFormula))
            doLog(0, "NOTE: there is no formula data for available this suspect.\n")
        else
            doLog(0, "NOTE: the suspect formula could not be matched with the formula annotations.\n")
    }
    if (is.null(cRow))
    {
        if (is.null(suspectInChIKey1))
            doLog(0, "NOTE: there is no compound data (eg InChIKey) available for this suspect.\n")
        else
            doLog(0, "NOTE: the suspect compound could not be matched with the compound annotations.\n")
    }
    
    for (lvl in names(IDLevelRules))
    {
        doLog(0, "Checking level '%s'\n", lvl)
        if (checkLevelOK(IDLevelRules[[lvl]]))
        {
            doLog(0, "assigned level '%s'!\n", lvl)
            return(lvl)
        }
    }
    
    return(NA_character_)
}

# method definition for as.data.table, both non-sets and sets
doFGScrAsDataTable <- function(x, ..., collapseSuspects = ",", onlyHits = FALSE)
{
    return(doFGAsDataTable(x, ..., collapseSuspects = collapseSuspects, onlyHits = onlyHits))
}

mergeScreenInfoWithDT <- function(tab, scrInfo, collapseSuspects, onlyHits)
{
    scrInfo <- copy(scrInfo)
    setnames(scrInfo, names(scrInfo), paste0("susp_", names(scrInfo))) # add susp_ column prefixes
    
    if (!is.null(collapseSuspects))
    {
        scrInfo[, susp_name := paste0(susp_name, collapse = collapseSuspects), by = "susp_group"]
        # only keep unique and remove suspect specific columns
        scrInfo <- unique(scrInfo[, c("susp_group", "susp_name"), with = FALSE], by = "susp_group")
    }
    
    ret <- merge(tab, scrInfo, by.x = "group", by.y = "susp_group", all.x = !onlyHits, sort = FALSE,
                 allow.cartesian = is.null(collapseSuspects))
    return(ret)
}
