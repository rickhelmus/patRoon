#' @include main.R
#' @include compounds.R
#' @include feature_groups-set.R
#' @include mslibrary.R
NULL

unifyLibNames <- function(cTab)
{
    unNames <- c(Name = "compoundName",
                 Synon = "compoundName2",
                 DB_ID = "identifier",
                 SMILES = "SMILES",
                 InChI = "InChI",
                 InChIKey = "InChIKey",
                 formula = "neutral_formula",
                 neutralMass = "neutralMass",
                 molNeutralized = "molNeutralized",
                 Precursor_type = "precursorType",
                 Spectrum_type = "spectrumType",
                 PrecursorMZ = "ion_formula_mz",
                 Instrument_Type = "instrumentType"
    )
    
    unNames <- unNames[names(unNames) %in% names(cTab)] # filter out missing
    setnames(cTab, names(unNames), unNames)
    
    return(cTab[, unNames, with = FALSE]) # filter out any other columns
}


#' Compound annotation with an MS library
#'
#' Uses a MS library loaded by \code{\link{loadMSLibrary}} for compound annotation.
#'
#' @templateVar algo MS library spectra
#' @templateVar do generate compound candidates
#' @templateVar generic generateCompounds
#' @templateVar algoParam library
#' @template algo_generator
#'
#' @details This method matches measured MS/MS data (peak lists) with those from an MS library to find candidate
#'   structures. Hence, only feature groups with MS/MS peak list data are annotated.
#'
#'   The library is searched for candidates with the following criteria: \enumerate{
#'
#'   \item Only records with ion \emph{m/z} (\code{PrecursorMZ}), \acronym{SMILES}, \acronym{InChI}, \acronym{InChIKey}
#'   and \code{formula} data are considered.
#'
#'   \item Depending on the value of the \code{checkIons} argument, records with different adduct
#'   (\code{Precursor_type}) or polarity (\code{Ion_mode}) may be ignored.
#'
#'   \item The \emph{m/z} values of the candidate and feature group should match (tolerance set by \code{absMzDev}
#'   argument).
#'
#'   \item The spectral similarity should not be lower than the value defined for the \code{minSim} argument.
#'
#'   \item If multiple candidates with the same first-block \acronym{InChIKey} are found then only the candidate with
#'   the best spectral match is kept.
#'
#'   }
#'
#'   If the library contains annotations these will be added to the matched MS/MS peaks. However, since the candidate
#'   selected from criterion #5 above may not contain all the annotation data available from the MS library, annotations
#'   from other records are also considered (controlled by the \code{minAnnSim} argument). If this leads to different
#'   annotations for the same mass peak then only the most abundant annotation is kept.
#'
#' @param MSLibrary The \code{\link{MSLibrary}} object that should be used to find candidates.
#' @param minSim The minimum spectral similarity for candidate records.
#' @param minAnnSim The minimum spectral similarity of a record for it to be used to find annotations (see the
#'   \verb{Details} section).
#' @param absMzDev The maximum absolute \emph{m/z} deviation between the feature group and library record \emph{m/z}
#'   values for candidate selection.
#' @param checkIons A \code{character} that excludes library records with different adduct (\code{checkIons="adduct"})
#'   or MS ionization polarity (\code{checkIons="polarity"}). If \code{checkIons="none"} then these filters are not
#'   applied.
#' @param specSimParamsLib Like \code{specSimParams}, but these parameters \emph{only} influence pre-treatment of
#'   library spectra (only the \code{removePrecursor}, \code{relMinIntensity} and \code{minPeaks} parameters are used).
#'
#' @template adduct-arg
#' @template specSimParams-arg
#' @template spectrumType-arg
#'
#' @inheritParams generateCompounds
#'
#' @seealso \code{\link{loadMSLibrary}} to obtain MS library data and the methods for \code{\link{MSLibrary}} to treat
#'   the data before using it for annotation.
#'
#' @templateVar what generateCompoundsLibrary
#' @template main-rd-method
#' @export
setMethod("generateCompoundsLibrary", "featureGroups", function(fGroups, MSPeakLists,
                                                                specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
                                                                MSLibrary, minSim = 0.75,
                                                                minAnnSim = minSim, absMzDev = 0.002, adduct = NULL,
                                                                checkIons = "adduct", spectrumType = "MS2",
                                                                specSimParamsMatch = getDefSpecSimParams(),
                                                                specSimParamsLib = getDefSpecSimParams())
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    aapply(assertSpecSimParams, . ~ specSimParams + specSimParamsMatch + specSimParamsLib, fixed = list(add = ac))
    checkmate::assertClass(MSLibrary, "MSLibrary", add = ac)
    aapply(checkmate::assertNumber, . ~ minSim + minAnnSim + absMzDev, lower = 0, finite = TRUE, fixed = list(add = ac))
    checkmate::assertChoice(checkIons, c("adduct", "polarity", "none"), add = ac)
    checkmate::assertCharacter(spectrumType, min.len = 1, min.chars = 1, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (specSimParamsMatch$shift != "none")
        stop("Spectral shifting not supported", call. = FALSE)
    
    if (length(fGroups) == 0)
        return(compounds(algorithm = "library"))
    
    adduct <- checkAndToAdduct(adduct, fGroups)
    
    gCount <- length(fGroups)
    gInfo <- groupInfo(fGroups)
    annTbl <- annotations(fGroups)
    libRecs <- records(MSLibrary)
    libSpecs <- spectra(MSLibrary)

    libRecs <- libRecs[!is.na(PrecursorMZ) & !is.na(SMILES) & !is.na(InChI) & !is.na(InChIKey) & !is.na(formula)]
    if (checkIons == "adduct")
        libRecs <- libRecs[!is.na(Precursor_type)]
    else if (checkIons == "polarity")
        libRecs <- libRecs[!is.na(Ion_mode)]
    
    getRecsForAdduct <- function(recs, add, addChr)
    {
        if (checkIons == "none")
            return(recs)
        else if (checkIons == "polarity")
        {
            pos <- add@charge > 0
            recs <- recs[((pos & Ion_mode == "POSITIVE") | (!pos & Ion_mode == "NEGATIVE"))]
        }
        else # if (checkIons == "adduct")
            recs <- recs[Precursor_type == addChr]
        return(recs)
    }
    
    if (!is.null(adduct))
        libRecs <- getRecsForAdduct(libRecs, adduct, as.character(adduct))
    else
        allAdducts <- sapply(unique(annTbl$adduct), as.adduct)

    if (!is.null(spectrumType))
        libRecs <- libRecs[Spectrum_type %chin% spectrumType]
    
    cacheDB <- openCacheDBScope()
    baseHash <- makeHash(minSim, minAnnSim, absMzDev, adduct, checkIons, specSimParamsMatch, specSimParamsLib)
    setHash <- makeHash(fGroups, MSPeakLists, MSLibrary, baseHash)
    cachedSet <- loadCacheSet("compoundsLibrary", setHash, cacheDB)
    resultHashes <- vector("character", gCount)
    resultHashCount <- 0
    
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
        
        spec <- prepSpecSimilarityPL(spec, removePrecursor = specSimParamsMatch$removePrecursor,
                                     relMinIntensity = specSimParamsMatch$relMinIntensity,
                                     minPeaks = specSimParamsMatch$minPeaks)
        if (nrow(spec) == 0)
            return(NULL)
        
        cTab <- libRecs[numLTE(abs(precMZ - PrecursorMZ), absMzDev)]
        if (is.null(adduct))
        {
            addChr <- annTbl[group == grp]$adduct
            cTab <- getRecsForAdduct(cTab, allAdducts[[addChr]], addChr)
        }
        
        if (nrow(cTab) == 0)
            return(NULL)

        hash <- makeHash(baseHash, spec, cTab, libSpecs[cTab$identifier])
        resultHashCount <<- resultHashCount + 1
        resultHashes[resultHashCount] <<- hash
        
        cached <- NULL
        if (!is.null(cachedSet))
            cached <- cachedSet[[hash]]
        if (is.null(cached))
            cached <- loadCacheData("compoundsLibrary", hash, cacheDB)
        if (!is.null(cached))
            return(cached)
        
        cTab <- unifyLibNames(cTab)
        cTab[, InChIKey1 := getIKBlock1(InChIKey)]
        lspecs <- Map(libSpecs[cTab$identifier], cTab$ion_formula_mz, cTab$identifier, f = function(sp, pmz, lid)
        {
            # convert to MSPeakLists format
            ret <- copy(sp)
            ret[, ID := seq_len(.N)]
            ret <- assignPrecursorToMSPeakList(ret, pmz)
            ret <- prepSpecSimilarityPL(ret, removePrecursor = specSimParamsLib$removePrecursor,
                                        relMinIntensity = specSimParamsLib$relMinIntensity,
                                        minPeaks = specSimParamsLib$minPeaks)
            return(ret)
        })
        lspecs <- pruneList(lspecs, checkZeroRows = TRUE)
        cTab <- cTab[identifier %in% names(lspecs)]
        if (nrow(cTab) == 0)
            return(NULL)
        
        sims <- specDistRect(list(spec), lspecs, specSimParamsMatch$method, specSimParamsMatch$shift, 0,
                             0, specSimParamsMatch$mzWeight, specSimParamsMatch$intWeight, specSimParamsMatch$absMzDev)
        
        cTab[, c("score", "libMatch") := sims[1, ]]
        
        hasAnnons <- FALSE
        for (sp in lspecs)
        {
            if (!is.null(sp[["annotation"]]))
            {
                hasAnnons <- TRUE
                break
            }
        }
        
        if (hasAnnons)
        {
            cTabAnn <- cTab[numGTE(score, minAnnSim)] # store before filtering/de-duplicating
            # needed below for neutral loss calculation
            thisAdduct <- if (!is.null(adduct)) adduct else as.adduct(annTbl[group == grp]$adduct)
        }
        
        cTab <- cTab[numGTE(score, minSim)]
        
        setorderv(cTab, "score", -1)
        
        cTab <- unique(cTab, by = "InChIKey1") # NOTE: prior sorting ensure top ranked stays

        # fill in fragInfos
        cTab[, fragInfo := list(Map(lspecs[identifier], InChIKey1, neutral_formula, f = function(ls, ik1, form)
        {
            bsp <- as.data.table(binSpectra(spec, ls, "none", 0, specSimParamsMatch$absMzDev))
            bsp <- bsp[intensity_1 != 0 & intensity_2 != 0] # overlap
            
            # NOTE: the mz values from the binned spectra could be slightly different --> take the original values
            fi <- data.table(mz = spec[match(bsp$ID_1, ID)]$mz, PLID = bsp$ID_1)
            
            fi[, c("ion_formula", "neutral_loss") := NA_character_]
            setorderv(fi, "PLID")
            
            if (hasAnnons)
            {
                cta <- cTabAnn[InChIKey1 == ik1]
                specMatched <- spec[ID %in% fi$PLID] # get peak list with only matched peaks, we need it for binning
                ann <- rbindlist(lapply(lspecs[cta$identifier], function(lsp)
                {
                    if (is.null(lsp[["annotation"]]))
                        return(NULL)
                    
                    # verify annotations
                    lsp <- lsp[verifyFormulas(annotation)]
                    
                    if (nrow(lsp) == 0)
                        return(NULL)
                    
                    # find overlapping peaks
                    bsp2 <- as.data.table(binSpectra(specMatched, lsp, "none", 0, specSimParamsMatch$absMzDev))
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
        
        saveCacheData("compoundsLibrary", cTab, hash, cacheDB)
        
        return(cTab)
    }, simplify = FALSE))
    
    if (is.null(cachedSet))
        saveCacheSet("compoundsLibrary", resultHashes[seq_len(resultHashCount)], setHash, cacheDB)
    
    compList <- pruneList(compList, checkZeroRows = TRUE)
    
    ngrp <- length(compList)
    printf("Loaded %d compounds from %d features (%.2f%%).\n", sum(unlist(lapply(compList, nrow))),
           ngrp, if (gCount == 0) 0 else ngrp * 100 / gCount)
    
    return(compounds(groupAnnotations = compList, scoreTypes = c("score", "libMatch"),
                     scoreRanges = sapply(compList, function(ct) list(score = range(ct$score),
                                                                      libMatch = range(ct$libMatch)), simplify = FALSE),
                     algorithm = "library", MSPeakLists = MSPeakLists, specSimParams = specSimParams))
})

#' @template featAnnSets-gen_args
#' @param \dots \setsWF Further arguments passed to the non-sets workflow method.
#' @rdname generateCompoundsLibrary
#' @export
setMethod("generateCompoundsLibrary", "featureGroupsSet", function(fGroups, MSPeakLists,
                                                                   specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
                                                                   MSLibrary, minSim = 0.75,
                                                                   minAnnSim = minSim, absMzDev = 0.002, adduct = NULL,
                                                                   ..., setThreshold = 0, setThresholdAnn = 0,
                                                                   setAvgSpecificScores = FALSE)
{
    generateCompoundsSet(fGroups, MSPeakLists, specSimParams, adduct, generateCompoundsLibrary, MSLibrary = MSLibrary,
                         minSim = minSim, minAnnSim = minAnnSim, absMzDev = absMzDev, ..., setThreshold = setThreshold,
                         setThresholdAnn = setThresholdAnn, setAvgSpecificScores = setAvgSpecificScores)
})
