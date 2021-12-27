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
                                                                absMzDev = 0.002, adduct = NULL, checkIons = "adduct",
                                                                specSimParams = getDefSpecSimParams())
{
    # UNDONE: cache
    # UNDONE: show mirror spectrum in report? Would need library data somehow
    # UNDONE: don't normalize scores (or already not done?)
    # UNDONE: separate specSimParams for lib? E.g. to assume that lib spectra are cleaner and don't need intensity cleaning

    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertClass(MSLibrary, "MSLibrary", add = ac)
    aapply(checkmate::assertNumber, . ~ minSim + absMzDev, lower = 0, finite = TRUE, fixed = list(add = ac))
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
    
    # UNDONE: ms level? or just a filter?
    # UNDONE: support entries without SMILES/InChI(keys)/Formulas?
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
    
    printf("Processing %d feature groups...\n", gCount)
    
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
        lspecs <- Map(libSpecs[cTab$identifier], cTab$ion_formula_mz, f = function(sp, pmz)
        {
            # convert to MSPeakLists format
            ret <- as.data.table(sp)
            ret[, ID := seq_len(.N)]
            ret <- assignPrecursorToMSPeakList(ret, pmz)
            ret <- prepSpecSimilarityPL(ret, removePrecursor = specSimParams$removePrecursor,
                                        relMinIntensity = specSimParams$relMinIntensity, minPeaks = specSimParams$minPeaks)
            return(ret)
        })
        lspecs <- pruneList(lspecs, checkZeroRows = TRUE)
        cTab <- cTab[identifier %in% names(lspecs)]
        if (nrow(cTab) == 0) # UNDONE: allow results without annotations (i.e. like MF)? Would interfere with min sim score though
            return(NULL)
        
        sims <- specDistRect(list(spec), lspecs, specSimParams$method, specSimParams$shift, 0,
                             0, specSimParams$mzWeight, specSimParams$intWeight, specSimParams$absMzDev)
        
        cTab[, score := sims[1, ]]
        
        cTab <- cTab[numGTE(score, minSim)]
        
        setorderv(cTab, "score", -1)
        
        # UNDONE: make optional and specify on which column to collapse
        cTab <- unique(cTab, by = "InChIKey") # NOTE: prior sorting ensure top ranked stays

        # fill in fragInfos
        # UNDONE: support libraries with annotations
        cTab[, fragInfo := list(lapply(lspecs[identifier], function(ls)
        {
            bsp <- as.data.table(binSpectra(spec, ls, "none", 0, specSimParams$absMzDev))
            bsp <- bsp[intensity_1 != 0 & intensity_2 != 0] # overlap
            
            # NOTE: the mz values from the binned spectra could be slightly different --> take the original values
            fi <- data.table(mz = spec[match(bsp$ID_1, ID)]$mz, PLID = bsp$ID_1)
            
            fi[, c("ion_formula", "neutral_loss") := NA_character_]
            setorderv(fi, "PLID")
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
    
    return(compounds(groupAnnotations = compList, scoreTypes = "score",
                     scoreRanges = sapply(compList, function(ct) list(score = range(ct$score)), simplify = FALSE),
                     algorithm = "library"))
})

if (F) {
#' @rdname compound-generation
#' @export
setMethod("generateCompoundsSIRIUS", "featureGroupsSet", function(fGroups, MSPeakLists, relMzDev = 5, adduct = NULL,
                                                                  ..., setThreshold = 0, setThresholdAnn = 0)
{
    generateCompoundsSet(fGroups, MSPeakLists, adduct, generateCompoundsSIRIUS, relMzDev = relMzDev, ...,
                         setThreshold = setThreshold, setThresholdAnn = setThresholdAnn)
})
}