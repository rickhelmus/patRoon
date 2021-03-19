#' @include utils.R
#' @include utils-compounds.R
NULL

isScreening <- function(fGroups) inherits(fGroups, c("featureGroupsScreening", "featureGroupsScreeningSet"))

getAllSuspSetCols <- function(targetCols, allCols, sets)
{
    targetCols <- c(targetCols, sapply(targetCols, function(cl) paste0(cl, "-", sets),
                                       USE.NAMES = FALSE))
    return(intersect(targetCols, allCols))
}

convertSuspDataIfNeeded <- function(scr, destFormat, destCol, fromFormats, fromCols)
{
    hasData <- function(x) !is.na(x) & nzchar(x)
    missingInScr <- function(x) if (is.null(scr[[x]])) rep(TRUE, nrow(scr)) else !hasData(scr[[x]])

    countEntries <- function() if (is.null(scr[[destCol]])) 0 else sum(hasData(scr[[destCol]]))
    curEntryCount <- countEntries()
    if (curEntryCount < nrow(scr))
    {
        printf("Trying to calculate missing %s data in suspect list... ", destCol)
        
        if (destFormat == "formula")
            doConv <- function(inp, f) convertToFormulaBabel(inp, f, mustWork = FALSE)
        else
            doConv <- function(inp, f) babelConvert(inp, f, destFormat, mustWork = FALSE)
        
        for (i in seq_along(fromFormats))
        {
            if (!is.null(scr[[fromCols[i]]]))
                scr[missingInScr(destCol) & !missingInScr(fromCols[i]), (destCol) := doConv(get(fromCols[i]), fromFormats[i])]
        }
     
        newEntryCount <- countEntries() - curEntryCount
        printf("Done! Filled in %d (%.1f%%) entries.\n", newEntryCount,
               if (newEntryCount > 0) newEntryCount * 100 / nrow(scr) else 0)
    }
    return(scr)
}

prepareSuspectList <- function(suspects, adduct, skipInvalid, calcMZs = TRUE)
{
    hash <- makeHash(suspects, adduct, skipInvalid, calcMZs)
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

        # make name column file safe and unique
        sanNames <- fs::path_sanitize(suspects$name, replacement = "_")
        sanNames <- make.unique(sanNames, sep = "-")
        changedNames <- which(sanNames != suspects$name)
        if (length(changedNames) > 0)
        {
            warning(paste0("The following suspect names were changed to make them file compatible and/or unique:\n",
                           paste0(suspects$name[changedNames], " --> ", sanNames[changedNames], collapse = "\n")))
            suspects[, name_orig := name]
            suspects[, name := sanNames]
        }
        
        # get missing identifiers & formulae if necessary and possible
        suspects <- convertSuspDataIfNeeded(suspects, destFormat = "smi", destCol = "SMILES",
                                            fromFormats = "inchi", fromCols = "InChI")
        suspects <- convertSuspDataIfNeeded(suspects, destFormat = "inchi", destCol = "InChI",
                                            fromFormats = "smi", fromCols = "SMILES")
        suspects <- convertSuspDataIfNeeded(suspects, destFormat = "inchikey", destCol = "InChIKey",
                                            fromFormats = c("smi", "inchi"), fromCols = c("SMILES", "InChI"))
        suspects <- convertSuspDataIfNeeded(suspects, destFormat = "formula", destCol = "formula",
                                            fromFormats = c("smi", "inchi"), fromCols = c("SMILES", "InChI"))
        
        # neutral masses given for all?
        if (!is.null(suspects[["neutralMass"]]) && !any(is.na(suspects[["neutralMass"]])))
            neutralMasses <- suspects[["neutralMass"]]
        else
        {
            printf("Calculating neutral masses for each suspect...\n")
            prog <- openProgBar(0, nrow(suspects))
            
            canUse <- function(v) !is.null(v) && !is.na(v) && (!is.character(v) || nzchar(v))
            neutralMasses <- sapply(seq_len(nrow(suspects)), function(i)
            {
                if (canUse(suspects[["neutralMass"]][i]))
                    ret <- suspects$neutralMass[i]
                else if (canUse(suspects[["formula"]][i]))
                    ret <- getFormulaMass(suspects$formula[i])
                else if (canUse(suspects[["SMILES"]][i]))
                    ret <- getNeutralMassFromSMILES(suspects$SMILES[i], mustWork = FALSE)[[1]]
                else
                    ret <- NA
                
                setTxtProgressBar(prog, i)
                return(ret)
            })
            
            close(prog)
        }

        # calculate ionic masses if possible (not possible if no adducts are given and fGroups are annotated)
        if (calcMZs && (is.null(suspects[["mz"]]) || any(is.na(suspects[["mz"]]))) &&
            (!is.null(adduct) || !is.null(suspects[["adduct"]])))
        {
            if (!is.null(adduct))
                addMZs <- adductMZDelta(adduct)
            else
                addMZs <- sapply(suspects[["adduct"]], function(a) adductMZDelta(as.adduct(a)))
            
            if (!is.null(suspects[["mz"]]))
                suspects[, mz := ifelse(!is.na(suspects$mz), suspects$mz, neutralMasses + addMZs)]
            else
                suspects[, mz := neutralMasses + addMZs]
        }
        else if (is.null(suspects[["mz"]]))
        {
            # NOTE: if mz column is already available it either contains user values or already NAs
            suspects[, mz := NA_real_]
        }
        
        suspects[, neutralMass := neutralMasses]
        
        saveCacheData("screenSuspectsPrepList", suspects, hash)
    }        
    
    if (!calcMZs)
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
    
    setMetaData <- function(t, suspRow)
    {
        for (col in c("name", "name_orig", "rt", "mz", "SMILES", "InChI", "InChIKey", "formula", "neutralMass", "adduct",
                      "fragments_mz", "fragments_formula"))
        {
            if (!is.null(suspects[[col]]))
                set(t, 1L, col, suspRow[[col]])
            else if (col == "rt")
                set(t, 1L, col, NA_real_) # exception: always want this column
        }
        return(t)
    }        
    
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
        
        if (nrow(gi) == 0) # no results? --> add NA result
        {
            ret <- data.table()
            setMetaData(ret, suspects[ti])
            ret[, c("group", "d_rt", "d_mz") := list(NA_character_, NA_real_, NA_real_)]
            return(ret)
        }
        
        hits <- rbindlist(lapply(rownames(gi), function(g)
        {
            ret <- data.table()
            setMetaData(ret, suspects[ti])
            ret[, c("group", "d_rt", "d_mz") := .(g, d_rt = if (hasRT) gi[g, "rts"] - rt else NA_real_,
                                                  ifelse(is.na(mz), annTbl[group == g]$neutralMass - neutralMass,
                                                         gi[g, "mzs"] - mz))]
            return(ret)
        }), fill = TRUE)
        
        setTxtProgressBar(prog, ti)
        return(hits)
    })
    ret <- rbindlist(retlist, fill = TRUE)
    
    setcolorder(ret, "name")
    
    setTxtProgressBar(prog, nrow(suspects))
    close(prog)
    
    suspectsn <- nrow(suspects)
    foundn <- suspectsn - sum(is.na(ret$group))
    printf("Found %d/%d suspects (%.2f%%)\n", foundn, suspectsn, foundn * 100 / suspectsn)
    
    ret <- ret[!is.na(group), ]
    
    return(ret[])
}

annotatedMSMSSimilarity <- function(annPL, absMzDev, relMinIntensity, method)
{
    if (is.null(annPL[["formula"]]))
        return(0)

    # remove precursor, as eg MetFrag doesn't include this and it's not so interesting anyway
    annPL <- annPL[precursor == FALSE]
    
    if (nrow(annPL) > 0)
    {
        maxInt <- max(annPL$intensity)
        annPL <- annPL[(intensity / maxInt) >= relMinIntensity]
    }
    
    if (nrow(annPL) == 0 || sum(!is.na(annPL$formula)) == 0)
        return(0)
    
    if (method == "cosine")
    {
        annPL[, intensityAnn := ifelse(is.na(formula), 0, intensity)]
        return(as.vector((annPL$intensityAnn %*% annPL$intensity) /
                             (sqrt(sum(annPL$intensityAnn^2)) * sqrt(sum(annPL$intensity^2)))))
    }
    else # jaccard
        return(sum(!is.na(annPL$formula)) / nrow(annPL))
}

estimateIdentificationLevel <- function(suspectName, suspectFGroup, suspectRTDev, suspectInChIKey1, suspectFormula,
                                        suspectAnnSimForm, suspectAnnSimComp, suspectAnnSimBoth,
                                        maxSuspFrags, maxFragMatches, formTable, formRank, mFormNames,
                                        formScoreRanges, formulasNormalizeScores, compTable,
                                        compRank, mCompNames, compScoreRanges, compoundsNormalizeScores,
                                        absMzDev, IDLevelRules)
{
    fRow <- cRow <- NULL
    if (!is.null(formTable) && !is.null(suspectFormula))
    {
        formTableNorm <- normalizeFormScores(formTable, formScoreRanges, mFormNames,
                                             formulasNormalizeScores == "minmax")
        unFTable <- unique(formTable, by = "formula"); unFTableNorm <- unique(formTableNorm, by = "formula")
        if (!is.na(formRank))
        {
            fRow <- unFTable[formRank]
            fRowNorm <- unFTableNorm[formRank]
        }
    }
    
    if (!is.null(compTable) && !is.null(suspectInChIKey1))
    {
        compTableNorm <- normalizeCompScores(compTable, compScoreRanges, mCompNames, compoundsNormalizeScores == "minmax")
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
        
        scoreVal <- rowMeans(annRow[, scCols, with = FALSE])
        if (scoreVal < minValue)
            return(sprintf("(average) score too low: %f/%f", scoreVal, minValue))
        
        htn <- getOptVal(val, "higherThanNext", 0)
        if (htn > 0 && nrow(annTable) > 1)
        {
            otherHighest <- max(rowMeans(annTable[-rank, scCols, with = FALSE]))
            if (is.infinite(htn)) # special case: should be highest
            {
                if (otherHighest > 0)
                    return("not the highest score")
            }
            else if ((scoreVal - otherHighest) < htn)
                return(sprintf("difference with highest score from other candidates too low: %f/%f", scoreVal - otherHighest, htn))
        }
        
        return(TRUE)            
    }

    logDir <- R.utils::getAbsolutePath(file.path("log", "ident")) # take absolute path for length calculation below
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
    
    formScores <- formulaScorings()$name
    compScores <- compoundScorings()$name
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
            
            if (type == "or")
            {
                if (!is.list(val) || checkmate::testNamed(val))
                    stop("Specify a list with 'or'")
                levelOK <- any(mapply(val, seq_along(val), FUN = function(IDL, i)
                {
                    doLog(indent + 1, "check OR condition %d/%d\n", i, length(val))
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
                    levelOK <- checkAnnotationScore(val, type, formRank, fRow, unFTable, fRowNorm, unFTableNorm,
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

mergeScreenInfoWithDT <- function(tab, scrInfo, collapseSuspects, onlyHits)
{
    scrInfo <- copy(scrInfo)
    setnames(scrInfo, names(scrInfo), paste0("susp_", names(scrInfo))) # add susp_ column prefixes
    
    if (!is.null(collapseSuspects))
    {
        scrInfo[, susp_name := paste0(susp_name, collapse = collapseSuspects), by = "susp_group"]
        # only keep unique and remove suspect specific columns
        # UNDONE: keep specific columns if only one suspect?
        scrInfo <- unique(scrInfo[, c("susp_group", "susp_name"), with = FALSE], by = "susp_group")
    }
    
    return(merge(tab, scrInfo, by.x = "group", by.y = "susp_group", all.x = !onlyHits, sort = FALSE))
}
