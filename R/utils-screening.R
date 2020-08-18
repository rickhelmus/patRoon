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
        printf("Done! Filled in %d (%.1f) entries.\n", newEntryCount,
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
    
    ret[, name_unique := make.unique(name)] # UNDONE: do we still need this?
    setcolorder(ret, c("name", "name_unique"))
    
    setTxtProgressBar(prog, nrow(suspects))
    close(prog)
    
    suspectsn <- nrow(suspects)
    foundn <- suspectsn - sum(is.na(ret$group))
    printf("Found %d/%d suspects (%.2f%%)\n", foundn, suspectsn, foundn * 100 / suspectsn)
    
    ret <- ret[!is.na(group), ]
    
    return(ret[])
}

annotatedMSMSSimilarity <- function(fragInfo, MSMSList, absMzDev, relMinIntensity)
{
    if (nrow(MSMSList) == 0 || nrow(fragInfo) == 0)
        return(0)
    
    MSMSList <- MSMSList[, c("mz", "intensity")]
    annMSMSList <- MSMSList[fragInfo$PLIndex]
    return(OrgMassSpecR::SpectrumSimilarity(annMSMSList, MSMSList, t = absMzDev,
                                            b = relMinIntensity, print.graphic = FALSE))
}

defaultIDLevelRules <- function(inLevels = NULL, exLevels = NULL)
{
    aapply(checkmate::assertCharacter, . ~ inLevels + exLevels, null.ok = TRUE)
    
    ret <- defIDLevelRules # stored inside R/sysdata.rda
    
    pred <- function(p, l) grepl(p, l)
    if (!is.null(inLevels))
        ret <- ret[grepl(inLevels, paste0(ret$level, ret$subLevel)), ]
    if (!is.null(exLevels))
        ret <- ret[!grepl(exLevels, paste0(ret$level, ret$subLevel)), ]
    return(ret)
}

# UNDONE/NOTE: mustExist/relative fields only used for scorings of compound/formulas
estimateIdentificationLevel <- function(suspectRTDev, suspectInChIKey1, suspectFormula, suspectAnnSim,
                                        suspectFragmentsMZ, suspectFragmentsForms,
                                        checkFragments, MSMSList,
                                        formTable, formScoreRanges, formulasNormalizeScores,
                                        compTable, mCompNames, compScoreRanges, compoundsNormalizeScores,
                                        absMzDev, IDLevelRules)
{
    if (!is.null(suspectFragmentsMZ))
        suspectFragmentsMZ <- as.numeric(unlist(strsplit(suspectFragmentsMZ, ";")))
    if (!is.null(suspectFragmentsForms))
        suspectFragmentsForms <- unlist(strsplit(suspectFragmentsForms, ";"))
    
    fRow <- cRow <- NULL
    if (!is.null(formTable) && !is.null(suspectFormula))
    {
        formTableNorm <- normalizeFormScores(formTable, formScoreRanges, formulasNormalizeScores == "minmax")
        unFTable <- unique(formTable, by = "formula"); unFTableNorm <- unique(formTableNorm, by = "formula")
        formRank <- which(suspectFormula == unFTable$neutral_formula)
        if (length(formRank) != 0)
        {
            formRank <- formRank[1]
            fRow <- unFTable[formRank]; fRowNorm <- unFTableNorm[formRank]
        }
    }
    
    if (!is.null(compTable) && !is.null(suspectInChIKey1))
    {
        compTableNorm <- normalizeCompScores(compTable, compScoreRanges, mCompNames, compoundsNormalizeScores == "minmax")
        compRank <- which(suspectInChIKey1 == compTable$InChIKey1)
        if (length(compRank) != 0)
        {
            compRank <- compRank[1]
            cRow <- compTable[compRank]; cRowNorm <- compTableNorm[compRank]
        }
    }
    
    if (!is.null(MSMSList))
        MSMSList <- MSMSList[precursor == FALSE]

    IDLevelRules <- if (is.data.table(IDLevelRules)) copy(IDLevelRules) else as.data.table(IDLevelRules)
    setorderv(IDLevelRules, c("level", "subLevel"))
    IDLevelList <- split(IDLevelRules, by = c("level", "subLevel"))

    mzWithin <- function(mz1, mz2) abs(mz1 - mz2) <= absMzDev
    
    checkAnnotationScore <- function(ID, rank, annRow, annTable, annRowNorm, annTableNorm, scCols)
    {
        # special case: rank
        if (ID$score == "rank")
            return(rank >= ID$value)
        
        scCols <- scCols[!is.na(unlist(annRow[, scCols, with = FALSE]))]
        if (length(scCols) == 0)
            return(!ID$mustExist)
        
        if (ID$relative)
        {
            annRow <- annRowNorm
            annTable <- annTableNorm
        }

        scoreVal <- rowMeans(annRow[, scCols, with = FALSE])
        if (scoreVal < ID$value)
            return(FALSE)
        
        if (!is.na(ID$higherThanNext) && ID$higherThanNext > 0 && nrow(annTable) > 1)
        {
            otherHighest <- max(rowMeans(annTable[-rank, scCols, with = FALSE]))
            if (is.infinite(ID$higherThanNext)) # special case: should be highest
            {
                if (otherHighest > 0)
                    return(FALSE)
            }
            else if ((scoreVal - otherHighest) < ID$higherThanNext)
                return(FALSE)
        }
        
        return(TRUE)            
    }
    checkScore <- function(ID)
    {
        if (ID$type == "retention" && ID$score == "maxDeviation")
            return(!is.na(suspectRTDev) && numLTE(abs(suspectRTDev), ID$value))
        if (ID$type == "formula")
            return(checkAnnotationScore(ID, formRank, fRow, unFTable, fRowNorm, unFTableNorm,
                                        getAllFormulasCols(ID$score, names(formTable))))
        if (ID$type == "compound")
            return(checkAnnotationScore(ID, compRank, cRow, compTable, cRowNorm, compTableNorm,
                                        getAllCompCols(ID$score, names(compTable), mCompNames)))
        if (ID$type == "suspectFragments")
        {
            suspMSMSMatchesMZ <- suspMSMSMatchesFormF <- suspMSMSMatchesFormC <- 0
            if (!is.null(suspectFragmentsMZ))
                suspMSMSMatchesMZ <- sum(sapply(MSMSList$mz, function(mz1) any(sapply(suspectFragmentsMZ, mzWithin, mz1 = mz1))))
            if (!is.null(suspectFragmentsForms))
            {
                if (!is.null(fRow) && "formula" %in% checkFragments)
                {
                    frTable <- formTable[byMSMS == TRUE & suspectFormula == neutral_formula]
                    if (nrow(frTable) > 0)
                    {
                        fi <- getFragmentInfoFromForms(MSMSList, frTable)
                        suspMSMSMatchesFormF <- sum(suspectFragmentsForms %in% fi$formula)
                    }
                }
                if (!is.null(cRow) && "compound" %in% checkFragments && !is.null(cRow[["fragInfo"]]))
                    suspMSMSMatchesFormC <- sum(suspectFragmentsForms %in% cRow$fragInfo$formula)
            }
            # UNDONE: make min(length...) configurable?
            return(max(suspMSMSMatchesMZ, suspMSMSMatchesFormF, suspMSMSMatchesFormC) >= min(ID$value, length(suspectFragmentsMZ)))
        }
        if (ID$type == "annotatedMSMSSimilarity")
            return(suspectAnnSim >= ID$value)
        stop(paste("Unknown ID level type:", ID$type))
    }

    for (IDL in IDLevelList)
    {
        if ("none" %in% IDL$type) # special case: always valid
            levelOK <- TRUE
        else
        {
            if ("suspectFragments" %in% IDL$type)
            {
                if (is.null(MSMSList) || nrow(MSMSList) == 0)
                    next
                hasFRMZ <- "mz" %in% checkFragments && !is.null(suspectFragmentsMZ) &&
                    !is.na(suspectFragmentsMZ) && length(suspectFragmentsMZ) > 0
                hasFRForms <- !is.null(suspectFragmentsForms) && !is.na(suspectFragmentsForms) &&
                    length(suspectFragmentsForms) > 0
                if (!hasFRMZ && !hasFRForms)
                    next
                if (!hasFRMZ &&
                    (is.null(formTable) || !"formula" %in% checkFragments) &&
                    (is.null(compTable) || !"compound" %in% checkFragments))
                    next
            }
            if ("retention" %in% IDL$type && is.null(suspectRTDev))
                next
            if ("formula" %in% IDL$type && (is.null(formTable) || nrow(formTable) == 0 ||
                                            is.null(suspectFormula) || is.null(fRow) ||
                                            nrow(fRow) == 0))
                next
            if (any(c("compound", "annotatedMSMSSimilarity") %in% IDL$type) &&
                (is.null(compTable) || nrow(compTable) == 0 || is.null(cRow) || nrow(cRow) == 0))
                next
            
            levelOK <- all(sapply(split(IDL, seq_len(nrow(IDL))), checkScore))
        }
        if (levelOK)
            return(paste0(IDL$level[1], IDL$subLevel[1]))
    }
    
    return(NA_character_)
}

