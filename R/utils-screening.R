#' @include utils.R
#' @include utils-compounds.R
NULL

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

prepareSuspectList <- function(suspects, adduct, skipInvalid)
{
    hash <- makeHash(suspects, adduct, skipInvalid)
    cd <- loadCacheData("screenSuspectsPrepList", hash)
    if (!is.null(cd))
        suspects <- cd
    else
    {
        # UNDONE: check if/make name column is file safe/unique
        
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
        
        # get missing identifiers & formulae if necessary and possible
        suspects <- convertSuspDataIfNeeded(suspects, destFormat = "smi", destCol = "SMILES",
                                            fromFormats = "inchi", fromCols = "InChI")
        suspects <- convertSuspDataIfNeeded(suspects, destFormat = "inchi", destCol = "InChI",
                                            fromFormats = "smi", fromCols = "SMILES")
        suspects <- convertSuspDataIfNeeded(suspects, destFormat = "inchikey", destCol = "InChIKey",
                                            fromFormats = c("smi", "inchi"), fromCols = c("SMILES", "InChI"))
        suspects <- convertSuspDataIfNeeded(suspects, destFormat = "formula", destCol = "formula",
                                            fromFormats = c("smi", "inchi"), fromCols = c("SMILES", "InChI"))
        
        if (!is.null(suspects[["mz"]]) && !any(is.na(suspects[["mz"]])))
        {
            saveCacheData("screenSuspectsPrepList", suspects, hash)
            return(suspects) # no further need for calculation of ion masses
        }
        
        # neutral masses given for all?
        if (!is.null(suspects[["neutralMass"]]) && !any(is.na(suspects[["neutralMass"]])))
            neutralMasses <- suspects[["neutralMass"]]
        else
        {
            printf("Calculating ion masses for each suspect...\n")
            prog <- openProgBar(0, nrow(suspects))
            
            canUse <- function(v) !is.null(v) && !is.na(v) && (!is.character(v) || nzchar(v))
            neutralMasses <- sapply(seq_len(nrow(suspects)), function(i)
            {
                if (canUse(suspects[["neutralMass"]][i]))
                    ret <- suspects$neutralMass[i]
                else if (canUse(suspects[["formula"]][i]))
                    ret <- rcdk::get.formula(suspects$formula[i])@mass
                else if (canUse(suspects[["SMILES"]][i]))
                    ret <- getNeutralMassFromSMILES(suspects$SMILES[i], mustWork = FALSE)[[1]]
                else
                    ret <- NA
                
                setTxtProgressBar(prog, i)
                return(ret)
            })
            
            close(prog)
        }
        
        if (!is.null(adduct))
            addMZs <- adductMZDelta(adduct)
        else
            addMZs <- sapply(suspects[["adduct"]], function(a) adductMZDelta(as.adduct(a)))
        
        if (!is.null(suspects[["mz"]]))
            suspects[, mz := ifelse(!is.na(suspects$mz), suspects$mz, neutralMasses + addMZs)]
        else
            suspects[, mz := neutralMasses + addMZs]
        
        suspects[, neutralMass := neutralMasses]
        
        saveCacheData("screenSuspectsPrepList", suspects, hash)
    }        
    
    # check for any suspects without proper mass info
    isNA <- is.na(suspects$mz)
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
    
    return(suspects)
}

doScreenSuspects <- function(fGroups, suspects, rtWindow, mzWindow, adduct, skipInvalid)
{
    gInfo <- groupInfo(fGroups)
    
    setMetaData <- function(t, suspRow)
    {
        for (col in c("name", "rt", "mz", "SMILES", "InChI", "InChIKey", "formula", "neutralMass", "adduct",
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
        
        # find related feature group(s)
        gi <- gInfo
        if (hasRT)
            gi <- gInfo[numLTE(abs(gInfo$rts - suspects$rt[ti]), rtWindow) & numLTE(abs(gInfo$mzs - suspects$mz[ti]), mzWindow), ]
        else
            gi <- gInfo[numLTE(abs(gInfo$mzs - suspects$mz[ti]), mzWindow), ]
        
        if (nrow(gi) == 0) # no results? --> add NA result
        {
            ret <- data.table()
            setMetaData(ret, suspects[ti])
            ret[, c("group", "d_rt", "d_mz") := NA]
            return(ret)
        }
        
        hits <- rbindlist(lapply(rownames(gi), function(g)
        {
            ret <- data.table()
            setMetaData(ret, suspects[ti])
            ret[, c("group", "d_rt", "d_mz") := .(g, d_rt = if (hasRT) gi[g, "rts"] - rt else NA, gi[g, "mzs"] - mz)]
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

# UNDONE: export?
numericIDLevel <- function(l) as.integer(gsub("[[:alpha:]]*", "", l))

genIDLevelRulesFile <- function(out, inLevels = NULL, exLevels = NULL)
{
    aapply(checkmate::assertCharacter, . ~ inLevels + exLevels, null.ok = TRUE)
    checkmate::assertPathForOutput(basename(out), overwrite = TRUE)
    
    defFile <- system.file("inst", "misc", "IDLevelRules.yml", package = "patRoon")
    
    if (is.null(inLevels) && is.null(exLevels))
        file.copy(defFile, out, overwrite = TRUE)
    else
    {
        rules <- yaml::yaml.load_file(defFile)
        if (!is.null(inLevels))
            rules <- rules[grepl(inLevels, names(rules))]
        if (!is.null(exLevels))
            rules <- rules[!grepl(exLevels, names(rules))]
        # UNDONE: this quotes ID levels without sub-level, fix?
        yaml::write_yaml(rules, out, indent = 4)
    }
    invisible(NULL)
}

# UNDONE/NOTE: mustExist/relative fields only used for scorings of compound/formulas
estimateIdentificationLevel <- function(suspectName, suspectFGroup, suspectRTDev, suspectInChIKey1, suspectFormula,
                                        suspectAnnSimForm, suspectAnnSimComp, suspectAnnSimBoth,
                                        maxSuspFrags, maxFragMatches, formTable, formRank,
                                        formScoreRanges, formulasNormalizeScores, compTable,
                                        compRank, mCompNames, compScoreRanges, compoundsNormalizeScores,
                                        absMzDev, IDLevelRules)
{
    fRow <- cRow <- NULL
    if (!is.null(formTable) && !is.null(suspectFormula))
    {
        formTableNorm <- normalizeFormScores(formTable, formScoreRanges, formulasNormalizeScores == "minmax")
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

    logOut <- file.path("log", "ident", paste0(suspectName, "-", suspectFGroup, ".txt"))
    mkdirp(dirname(logOut))
    unlink(logOut)
    doLog <- function(indent, s, ...) fprintf(logOut, paste0(strrep(" ", indent * 4), s), ..., append = TRUE)
    
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
                # UNDONE: if suspect fragments are less than the rule value then the
                # former is used as minimum, make this configurable?
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
                                                    getAllFormulasCols(type, names(formTable)))
                else
                    levelOK <- checkAnnotationScore(val, type, compRank, cRow, compTable, cRowNorm, compTableNorm,
                                                    getAllCompCols(type, names(compTable), mCompNames))
                
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
