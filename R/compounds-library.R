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
                 Formula = "formula",
                 ExactMass = "neutralMass",
                 Precursor_Type = "precursorType",
                 Spectrum_Type = "spectrumType",
                 PrecursorMZ = "precursorMZ",
                 Instrument_Type = "instrumentType"
    )
    
    unNames <- unNames[names(unNames) %in% names(cTab)] # filter out missing
    setnames(cTab, names(unNames), unNames)
    
    return(cTab[, unNames, with = FALSE]) # filter out any other columns
}


#' @export
setMethod("generateCompoundsLibrary", "featureGroups", function(fGroups, MSPeakLists, MSLibrary, absMzDev = 0.002,
                                                                adduct = NULL, specSimParams = getDefSpecSimParams())
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertClass(MSLibrary, "MSLibrary", add = ac)
    checkmate::assertNumber(absMzDev, lower = 0, finite = TRUE, add = ac)
    assertSpecSimParams(specSimParams, add = ac)
    checkmate::reportAssertions(ac)
    
    if (length(fGroups) == 0)
        return(compounds(algorithm = "library"))
    
    adduct <- checkAndToAdduct(adduct, fGroups)
    
    gCount <- length(fGroups)
    gInfo <- groupInfo(fGroups)
    annTbl <- annotations(fGroups)
    recTab <- records(MSLibrary)
    recTab <- recTab[!is.na(PrecursorMZ) & !is.na(Precursor_Type) & !is.na(Ion_mode)]
    
    getRecsForAdduct <- function(recs, add, addChr)
    {
        pos <- add@charge > 0
        recs[Precursor_Type == addChr & ((pos & Ion_mode == "POSITIVE") | (!pos & Ion_mode == "NEGATIVE"))]
    }
    
    if (!is.null(adduct))
        recTab <- getRecsForAdduct(recTab, adduct, as.character(adduct))
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
        
        cTab <- recTab[numLTE(abs(precMZ - PrecursorMZ), absMzDev)]
        
        if (is.null(adduct))
        {
            addChr <- annTbl[group == grp]$adduct
            recTab <- getRecsForAdduct(recTab, allAdducts[addChr], addChr)
        }
        
        if (nrow(cTab) == 0)
            return(NULL)

        cTab <- unifyLibNames(cTab)
        cTab[, InChIKey1 := getIKBlock1(InChIKey)]
        
        libSpecs <- lapply(spectra(MSLibrary)[cTab$identifier], function(sp)
        {
            # convert to MSPeakLists format
            ret <- as.data.table(sp)
            ret[, ID := seq_len(.N)]
            return(ret)
        })
        
        # UNDONE: ensure that no shift is applied in specSimParams
        sims <- specDistRect(list(spec), libSpecs, specSimParams$method, specSimParams$shift, 0,
                             0, specSimParams$mzWeight, specSimParams$intWeight, specSimParams$absMzDev)
        
        cTab[, score := sims[1, ]]
        
        cTab <- cTab[score >= 0.75] # UNDONE
        
        setorderv(cTab, "score", -1)
        
        # UNDONE: make optional and specify on which column to collapse
        cTab <- unique(cTab, by = "InChIKey") # NOTE: prior sorting ensure top ranked stays
        
        return(cTab)
    }, simplify = FALSE))
    
    compList <- pruneList(compList, checkZeroRows = TRUE)
    
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