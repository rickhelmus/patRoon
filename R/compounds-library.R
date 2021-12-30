#' @include main.R
#' @include compounds.R
#' @include feature_groups-set.R
#' @include MSLibrary.R
NULL

unifyLibNames <- function(cTab)
{
    unNames <- c(Name = "compoundName",
                 Synon = "compoundName2",
                 DB_ID = "identifier",
                 SMILES = "SMILES",
                 InChI = "InChI",
                 InChIKey = "InChIKey",
                 Formula = "neutral_formula",
                 ExactMass = "neutralMass",
                 Precursor_type = "precursorType",
                 Spectrum_Type = "spectrumType",
                 PrecursorMZ = "ion_formula_mz",
                 Instrument_Type = "instrumentType"
    )
    
    unNames <- unNames[names(unNames) %in% names(cTab)] # filter out missing
    setnames(cTab, names(unNames), unNames)
    
    return(cTab[, unNames, with = FALSE]) # filter out any other columns
}


#' @export
setMethod("generateCompoundsLibrary", "featureGroups", function(fGroups, MSPeakLists, MSLibrary, minSim = 0.75,
                                                                minAnnSim = minSim, absMzDev = 0.002, adduct = NULL,
                                                                checkIons = "adduct",
                                                                specSimParams = getDefSpecSimParams())
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertClass(MSLibrary, "MSLibrary", add = ac)
    aapply(checkmate::assertNumber, . ~ minSim + minAnnSim + absMzDev, lower = 0, finite = TRUE, fixed = list(add = ac))
    checkmate::assertChoice(checkIons, c("adduct", "polarity", "none"), add = ac)
    assertSpecSimParams(specSimParams, add = ac)
    checkmate::reportAssertions(ac)
    
    if (specSimParams$shift != "none")
        stop("Spectral shifting not supported", call. = FALSE)
    
    if (length(fGroups) == 0)
        return(compounds(algorithm = "library"))
    
    adduct <- checkAndToAdduct(adduct, fGroups)
    
    gCount <- length(fGroups)
    gInfo <- groupInfo(fGroups)
    annTbl <- annotations(fGroups)
    libRecs <- records(MSLibrary)
    libSpecs <- spectra(MSLibrary)
    libAnn <- annotations(MSLibrary)
    
    libRecs <- libRecs[!is.na(PrecursorMZ) & !is.na(SMILES) & !is.na(InChI) & !is.na(InChIKey) & !is.na(Formula)]
    if (checkIons == "adduct")
        libRecs <- libRecs[!is.na(Precursor_type)]
    if (checkIons != "none")
        libRecs <- libRecs[!is.na(Ion_mode)]
    
    getRecsForAdduct <- function(recs, add, addChr)
    {
        if (checkIons == "none")
            return(recs)
        pos <- add@charge > 0
        recs <- recs[((pos & Ion_mode == "POSITIVE") | (!pos & Ion_mode == "NEGATIVE"))]
        if (checkIons == "adduct")
            recs <- recs[Precursor_type == addChr]
        return(recs)
    }
    
    if (!is.null(adduct))
        libRecs <- getRecsForAdduct(libRecs, adduct, as.character(adduct))
    else
        allAdducts <- sapply(unique(annTbl$adduct), as.adduct)
    
    printf("Processing %d feature groups with a library of %d records...\n", gCount, nrow(libRecs))
    
    compList <- withProg(length(fGroups), FALSE, sapply(names(fGroups), function(grp)
    {
        # HACK: call here since there are quite a few early returns below
        doProgress()
        
        if (is.null(MSPeakLists[[grp]]) || is.null(MSPeakLists[[grp]][["MS"]]))
            return(NULL)
        spec <- MSPeakLists[[grp]][["MSMS"]]
        if (is.null(spec))
            return(NULL)
        
        precMZ <- MSPeakLists[[grp]]$MS[precursor == TRUE]$mz
        
        spec <- prepSpecSimilarityPL(spec, removePrecursor = specSimParams$removePrecursor,
                                     relMinIntensity = specSimParams$relMinIntensity, minPeaks = specSimParams$minPeaks)
        if (nrow(spec) == 0)
            return(NULL)
        
        cTab <- libRecs[numLTE(abs(precMZ - PrecursorMZ), absMzDev)]
        
        if (is.null(adduct))
        {
            addChr <- annTbl[group == grp]$adduct
            libRecs <- getRecsForAdduct(libRecs, allAdducts[addChr], addChr)
        }
        
        if (nrow(cTab) == 0)
            return(NULL)

        cTab <- unifyLibNames(cTab)
        cTab[, InChIKey1 := getIKBlock1(InChIKey)]
        lspecs <- Map(libSpecs[cTab$identifier], cTab$ion_formula_mz, cTab$identifier, f = function(sp, pmz, lid)
        {
            # convert to MSPeakLists format
            ret <- as.data.table(sp)
            ret[, ID := seq_len(.N)]
            ret <- assignPrecursorToMSPeakList(ret, pmz)
            ret <- prepSpecSimilarityPL(ret, removePrecursor = specSimParams$removePrecursor,
                                        relMinIntensity = specSimParams$relMinIntensity, minPeaks = specSimParams$minPeaks)
            
            # add annotations, if any
            if (!is.null(libAnn[[lid]]) && length(libAnn[[lid]]) > 0)
                ret[, annotation := libAnn[[lid]][ID]]
            
            return(ret)
        })
        lspecs <- pruneList(lspecs, checkZeroRows = TRUE)
        cTab <- cTab[identifier %in% names(lspecs)]
        if (nrow(cTab) == 0)
            return(NULL)
        
        sims <- specDistRect(list(spec), lspecs, specSimParams$method, specSimParams$shift, 0,
                             0, specSimParams$mzWeight, specSimParams$intWeight, specSimParams$absMzDev)
        
        cTab[, c("score", "libMatch") := sims[1, ]]
        
        if (length(libAnn) > 0)
        {
            cTabAnn <- cTab[numGTE(score, minAnnSim)] # store before filtering/de-duplicating
            # needed below for neutral loss calculation
            thisAdduct <- if (!is.null(adduct)) adduct else as.adduct(annTbl[group == grp]$adduct)
        }
        
        cTab <- cTab[numGTE(score, minSim)]
        
        setorderv(cTab, "score", -1)
        
        cTab <- unique(cTab, by = "InChIKey") # NOTE: prior sorting ensure top ranked stays

        # fill in fragInfos
        cTab[, fragInfo := list(Map(lspecs[identifier], InChIKey, neutral_formula, f = function(ls, ik, form)
        {
            bsp <- as.data.table(binSpectra(spec, ls, "none", 0, specSimParams$absMzDev))
            bsp <- bsp[intensity_1 != 0 & intensity_2 != 0] # overlap
            
            # NOTE: the mz values from the binned spectra could be slightly different --> take the original values
            fi <- data.table(mz = spec[match(bsp$ID_1, ID)]$mz, PLID = bsp$ID_1)
            
            fi[, c("ion_formula", "neutral_loss") := NA_character_]
            setorderv(fi, "PLID")
            
            if (length(libAnn) > 0)
            {
                cta <- cTabAnn[InChIKey == ik]
                specMatched <- spec[ID %in% fi$PLID] # get peak list with only matched peaks, we need it for binning
                ann <- rbindlist(lapply(lspecs[cta$identifier], function(lsp)
                {
                    if (is.null(lsp[["annotation"]]))
                        return(NULL)
                    
                    # fixup annotations
                    lsp <- copy(lsp)
                    lsp[, annotation := trimws(annotation)]
                    lsp[, annotation := sub("\\+|\\-$", "", annotation)] # remove any trailing charge
                    lsp <- lsp[verifyFormulas(annotation)]
                    
                    if (nrow(lsp) == 0)
                        return(NULL)
                    
                    # find overlapping peaks
                    bsp2 <- as.data.table(binSpectra(specMatched, lsp, "none", 0, specSimParams$absMzDev))
                    bsp2 <- bsp2[intensity_1 != 0 & intensity_2 != 0] # overlap
                    
                    if (nrow(bsp2) == 0)
                        return(NULL)

                    return(data.table(ID = bsp2$ID_1, annotation = lsp[match(bsp2$ID_2, ID)]$annotation))
                }))
                
                if (nrow(ann) > 0)
                {
                    # resolve conflicts in annotation: keep most abundant
                    ann[, count := .N, by = c("ID", "annotation")]
                    setorderv(ann, c("ID", "count"), order = c(1, -1))
                    ann <- unique(ann, by = "ID")
                    
                    fi[match(ann$ID, PLID), ion_formula := ann$annotation]
                    fi[!is.na(ion_formula), ion_formula := sub("\\+|\\-$", "", ion_formula)] # remove any trailing charge
                    ionform <- calculateIonFormula(form, thisAdduct)
                    fi[!is.na(ion_formula), neutral_loss := sapply(ion_formula, subtractFormula, formula1 = ionform)]
                }
            }
            
            return(fi)
        }))]
        
        cTab[, explainedPeaks := sapply(fragInfo, nrow)]
        cTab[, libPeaksCompared := sapply(lspecs[identifier], nrow)]
        cTab[, libPeaksTotal := sapply(libSpecs[identifier], nrow)]
        cTab[, database := "library"]
        
        return(cTab)
    }, simplify = FALSE))
    
    compList <- pruneList(compList, checkZeroRows = TRUE)
    
    ngrp <- length(compList)
    printf("Loaded %d compounds from %d features (%.2f%%).\n", sum(unlist(lapply(compList, nrow))),
           ngrp, if (gCount == 0) 0 else ngrp * 100 / gCount)
    
    return(compounds(groupAnnotations = compList, scoreTypes = c("score", "libMatch"),
                     scoreRanges = sapply(compList, function(ct) list(score = range(ct$score), libMatch = range(ct$libMatch)), simplify = FALSE),
                     algorithm = "library"))
})

#' @rdname compound-generation
#' @export
setMethod("generateCompoundsLibrary", "featureGroupsSet", function(fGroups, MSPeakLists, MSLibrary, minSim = 0.75,
                                                                   minAnnSim = minSim, absMzDev = 0.002, adduct = NULL,
                                                                  ..., setThreshold = 0, setThresholdAnn = 0)
{
    generateCompoundsSet(fGroups, MSPeakLists, adduct, generateCompoundsLibrary, MSLibrary = MSLibrary, minSim = minSim,
                         minAnnSim = minAnnSim, absMzDev = absMzDev, ..., setThreshold = setThreshold,
                         setThresholdAnn = setThresholdAnn)
})
