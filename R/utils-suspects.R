#' @include utils.R
#' @include utils-compounds.R
NULL

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
            InChISMILES <- NULL
            
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
                else if (canUse(suspects[["InChI"]][i]))
                {
                    # it's more efficient to calculate all at once with obabel --> cache result
                    if (is.null(InChISMILES))
                    {
                        doInChI <- suspects$InChI[!is.na(suspects$InChI) & nzchar(suspects$InChI)]
                        InChISMILES <<- babelConvert(doInChI, "inchi", "smi", mustWork = FALSE)
                        names(InChISMILES) <<- doInChI
                    }
                    ret <- getNeutralMassFromSMILES(InChISMILES[[suspects$InChI[i]]], mustWork = FALSE)[[1]]
                }
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

annotatedMSMSSimilarity <- function(fragInfo, MSMSList, absMzDev, relMinIntensity)
{
    if (nrow(MSMSList) == 0 || nrow(fragInfo) == 0)
        return(0)
    
    MSMSList <- MSMSList[, c("mz", "intensity")]
    annMSMSList <- MSMSList[fragInfo$PLIndex]
    return(OrgMassSpecR::SpectrumSimilarity(annMSMSList, MSMSList, t = absMzDev,
                                            b = relMinIntensity, print.graphic = FALSE))
}

estimateIdentificationLevel <- function(suspectInChIKey1, suspectFormula, suspectAnnSim,
                                        suspectFragments, MSMSList, formTable, formScoreRanges,
                                        minFormScores, minFormScoresToNext, compTable, absMzDev)
{
    # Level 2a: suspect is the only candidate with >= 0.9 MoNa score
    # Level 3a: suspect is with (multiple candidates of) MoNa >= 0.4
    # Level 3b: suspect matches >= 3 fragments (or all available) in suspect list
    # Level 3c: suspect has sufficient in-silico score --> MFSim >= 0.7
    #   - MFSim: similarity between reconstructed spectrum from annotated peaks and actual MS/MS spectrum
    # Level 4: Clear single formula result
    #   - score >= minScore && minScoreNext higher than next candidate (if multiple)
    # Level 5: Anything else
    
    if (!is.null(suspectFragments))
        suspectFragments <- as.numeric(unlist(strsplit(suspectFragments, ";")))
    
    # compound info for suspect
    cRow <- if (!is.null(compTable) && !is.null(suspectInChIKey1)) compTable[suspectInChIKey1 == InChIKey1] else
        data.table::data.table()
    
    if (nrow(cRow) != 0 && !is.null(cRow[["individualMoNAScore"]]))
    {
        goodMoNas <- compTable[individualMoNAScore >= 0.9]
        if (nrow(goodMoNas) == 1 && suspectInChIKey1 %in% goodMoNas$InChIKey1)
            return("2a")
        
        if (cRow$individualMoNAScore >= 0.4)
            return("3a")
    }
    
    if (!is.null(MSMSList))
        MSMSList <- MSMSList[precursor == FALSE]
    if (!is.null(MSMSList) && nrow(MSMSList) > 0 && length(suspectFragments) > 0)
    {
        suspMSMSMatches <- sapply(MSMSList$mz, function(mz1) any(sapply(suspectFragments, mzWithin, mz1 = mz1,
                                                                        absMzDev = absMzDev)))
        if (length(suspectFragments) > 0 && sum(suspMSMSMatches) >= min(3, length(suspectFragments)))
            return("3b")
    }
    
    if (suspectAnnSim >= 0.7)
        return("3c")
    
    # no (good) hit --> level 4/5
    if (!is.null(formTable) && !is.null(suspectFormula))
    {
        formTable <- normalizeFormScores(formTable, formScoreRanges, FALSE)
        unFTable <- unique(formTable, by = "formula")
        fRow <- unFTable[suspectFormula == neutral_formula]
        
        # suspect formula present and top ranked?
        if (suspectFormula %in% fRow$neutral_formula && suspectFormula == unFTable$neutral_formula[1])
        {
            if (nrow(unFTable) > 1) # multiple hits
            {
                avgScoreVals <- function(sc, fr)
                {
                    scCols <- getAllFormulasCols(sc, names(formTable))
                    if (length(scCols) == 0)
                        return(NA)
                    return(rowMeans(fr[, scCols, with = FALSE]))
                }
                    
                scoreVals <- sapply(names(minFormScores), avgScoreVals, fRow)
                nextScoreVals <- sapply(names(minFormScores), avgScoreVals, unFTable[2])
                
                if (all(is.na(scoreVals) | (scoreVals > minFormScores & (scoreVals - nextScoreVals) > minFormScoresToNext)))
                    return(4)
            }
            else
                return("4") # only hit
        }
    }
    
    return("5") # no (good) formula match
}

# NOTE: no-result rows are removed
# UNDONE: proper defaults for minFormScores/minFormScoresToNext
# UNDONE: allow compound scoring thresholds for ID level estimation
annotateSuspectList <- function(scr, MSPeakLists = NULL, formulas = NULL, compounds = NULL,
                                absMzDev = 0.005, relMinMSMSIntensity = 0.05,
                                minFormScores, minFormScoresToNext)
{
    ac <- checkmate::makeAssertCollection()
    assertScreeningResults(scr, fromFGroups = TRUE, add = ac)
    aapply(checkmate::assertClass, . ~ MSPeakLists + formulas + compounds,
           c("MSPeakLists", "formulas", "compounds"), null.ok = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ absMzDev + relMinMSMSIntensity, lower = 0,
           finite = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertNumeric, . ~ minFormScores + minFormScoresToNext, finite = TRUE,
           any.missing = FALSE, names = "unique", fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    hash <- makeHash(scr, MSPeakLists, formulas, compounds, absMzDev, relMinMSMSIntensity, minFormScores,
                     minFormScoresToNext)
    cd <- loadCacheData("annotateSuspects", hash)
    if (!is.null(cd))
        return(cd)
    
    scr <- copy(scr)
    scr <- scr[!is.na(group)]
    
    # get InChIKeys/Formulas if necessary and possible

    hasData <- function(x) !is.na(x) & nzchar(x)
    missingInScr <- function(what) if (is.null(scr[[what]])) rep(TRUE, nrow(scr)) else !hasData(scr[[what]])
    
    if (is.null(scr[["InChIKey"]]) || any(!hasData(scr$InChIKey)))
    {
        printf("Trying to calculate missing InChIKeys...\n")
        
        if (!is.null(scr[["SMILES"]]))
            scr[missingInScr("InChIKey") & !missingInScr("SMILES"), InChIKey := babelConvert(SMILES, "smi", "inchikey", mustWork = FALSE)]
        
        # re-try from InChI for the results not yet available
        if (!is.null(scr[["InChI"]]))
            scr[missingInScr("InChIKey") & !missingInScr("InChI"), InChIKey := babelConvert(InChI, "inchi", "inchikey", mustWork = FALSE)]
    }
    if (is.null(scr[["formula"]]) || any(!hasData(scr$formula)))
    {
        printf("Trying to calculate missing formulas...\n")
        missingSMILES <- missingInScr("SMILES")
        missingInChIs <- missingInScr("InChI")
        
        if (!is.null(scr[["SMILES"]]))
            scr[!missingSMILES & missingInScr("formula"), formula := convertToFormulaBabel(SMILES, "smi", mustWork = FALSE)]
        if (!is.null(scr[["InChI"]]))
            scr[!missingInChIs & missingInScr("formula"), formula := convertToFormulaBabel(SMILES, "inchi", mustWork = FALSE)]
    }
    
    for (i in seq_len(nrow(scr)))
    {
        gName <- scr$name[i]
        MSMSList <- if (!is.null(MSPeakLists)) MSPeakLists[[gName]][["MSMS"]] else NULL
        fTable <- if (!is.null(formulas)) formulas[[gName]] else NULL
        fScRanges <- if (!is.null(formulas)) formulas@scoreRanges[[gName]] else NULL
        cTable <- if (!is.null(compounds)) compounds[[gName]] else NULL
        
        suspIK1 <- if (!is.null(scr[["InChIKey"]])) getIKBlock1(scr$InChIKey[i]) else NULL
        annSim <- 0; suspRank <- NA
        if (!is.null(MSMSList) && !is.null(cTable) && !is.null(suspIK1))
        {
            suspRank <- which(suspIK1 == cTable$InChIKey1)
            suspRank <- if (length(suspRank) > 0) suspRank[1] else NA
            
            if (!is.na(suspRank) && !is.null(cTable[["fragInfo"]][[suspRank]]))
                annSim <- annotatedMSMSSimilarity(cTable[["fragInfo"]][[suspRank]],
                                                  MSPeakLists[[gName]][["MSMS"]],
                                                  absMzDev, relMinMSMSIntensity)
        }
        
        set(scr, i, c("suspCompAnnRank", "annotatedMSMSSimilarity"), list(suspRank, annSim))
        set(scr, i, "estIDLevel", estimateIdentificationLevel(suspIK1, scr$formula[i], annSim,
                                                              if (!is.null(scr[["fragments"]])) scr$fragments[i] else NULL,
                                                              MSMSList, fTable, fScRanges, minFormScores,
                                                              minFormScoresToNext, cTable, absMzDev))
    }
    
    saveCacheData("annotateSuspects", scr, hash)
    
    return(scr[])
}
