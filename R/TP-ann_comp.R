#' @include main.R
#' @include TP-structure.R
NULL

#' @rdname transformationProductsStructure-class
transformationProductsAnnComp <- setClass("transformationProductsAnnComp", contains = "transformationProductsStructure")

setMethod("initialize", "transformationProductsAnnComp",
          function(.Object, ...) callNextMethod(.Object, algorithm = "ann_comp", ...))

# NOTE: this function is called by a withProg() block, so handles progression updates
getTPsCompounds <- function(annTable, parName, parFormula, parSMILES, parLogP, extraOptsFMCSR, suspTPsSMILES, fpType,
                            fpSimMethod, minFitFormula, minFitCompound, minSimSusp, parallel)
{
    tab <- copy(annTable)
    setnames(tab, "neutral_formula", "formula")
    
    tab[, c("name", "ID", "parent_ID", "chem_ID", "generation") := .(paste0(parName, "-TP", .I), .I, NA_integer_, .I, 1)]
    
    # UNDONE: apply thresholds
    # UNDONE: unset annotation objects?
    
    # UNDONE: also support XLogP from feat annotations?
    # UNDONE: logP tolerance
    if (!is.null(annTable[["logP"]]))
    {
        tab[, retDir := fcase(logP < parLogP, -1,
                              logP > parLogP, 1,
                              default = 0)] # NOTE: default should also handle NAs
    }
    
    tab[, fitFormula := sapply(formula, calcFormulaFit, parFormula)]
    tab <- tab[numGTE(fitFormula, minFitFormula)]
    
    if (nrow(tab) > 0 && !is.null(parSMILES))
    {
        compFits <- do.call(calcStructFitFMCS, c(list(parSMILES, tab$SMILES, parallel), extraOptsFMCSR))
        tab[SMILES %chin% names(compFits), fitCompound := compFits[SMILES]]
        
        if (!is.null(suspTPsSMILES) && length(suspTPsSMILES) > 0)
        {
            tab[, c("suspSim", "suspSimSMILES") := rbindlist(doApply("lapply", parallel, SMILES, function(SMI)
            {
                dists <- sapply(suspTPsSMILES, patRoon:::distSMILES, SMI1 = SMI, fpType = fpType,
                                fpSimMethod = fpSimMethod)
                wh <- which.max(dists)
                return(list(dists[wh], suspTPsSMILES[wh]))
            }, prog = FALSE))]
            tab <- tab[numGTE(fitCompound, minFitCompound) | numGTE(suspSim, minSimSusp)]
        }
    }
    
    # UNDONE: also do TP score w/out suspSim (ie when TPs=NULL)?
    if (!is.null(tab[["fitCompound"]]) && !is.null(tab[["suspSim"]]))
        tab[, TP_score := pmax(NAToZero(fitCompound), NAToZero(suspSim)) + NAToZero(annSim)]
    else
        tab[, TP_score := NAToZero(fitFormula) + NAToZero(annSim)]
    
    tab <- subsetDTColumnsIfPresent(tab, c("group", "name", "ID", "parent_ID", "chem_ID", "generation", "compoundName",
                                           "SMILES", "InChI", "InChIKey", "formula", "logP", "retDir", "annSim",
                                           "fitFormula", "fitCompound", "suspSim", "suspSimSMILES", "TP_score"))
    
    doProgress()
    
    return(tab)
}


#' @export
generateTPsAnnComp <- function(parents, compounds, TPsRef = NULL, minFitFormula = 0, minFitCompound = 0, minSimSusp = 0,
                               extraOptsFMCSR = NULL, skipInvalid = TRUE, prefCalcChemProps = TRUE,
                               neutralChemProps = FALSE, calcLogP = "rcdk", calcSims = FALSE, fpType = "extended",
                               fpSimMethod = "tanimoto", parallel = TRUE)
{
    # UNDONE: support >1 generations? Probably not really worthwhile...
    
    checkmate::assert(
        checkmate::checkClass(parents, "data.frame"),
        checkmate::checkClass(parents, "compounds"), # UNDONE: doesn't make much sense here, but keep for consistency?
        checkmate::checkClass(parents, "featureGroupsScreening"),
        checkmate::checkClass(parents, "featureGroupsScreeningSet"),
        .var.name = "parents"
    )
    
    ac <- checkmate::makeAssertCollection()
    if (is.data.frame(parents))
        assertSuspectList(parents, needsAdduct = FALSE, skipInvalid = TRUE, add = ac)
    checkmate::assertClass(compounds, "compounds", add = ac)
    checkmate::assertClass(TPsRef, "transformationProductsStructure", null.ok = TRUE, add = ac)
    aapply(checkmate::assertNumber, . ~ minFitFormula + minFitCompound + minSimSusp, fixed = list(add = ac))
    checkmate::assertList(extraOptsFMCSR, null.ok = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ skipInvalid + prefCalcChemProps + neutralChemProps + calcSims + parallel,
           fixed = list(add = ac))
    assertXLogPMethod(calcLogP, add = ac)
    aapply(checkmate::assertString, . ~ fpType + fpSimMethod, min.chars = 1, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    parents <- getTPParents(parents, skipInvalid, prefCalcChemProps, neutralChemProps)
    
    if (nrow(parents) == 0)
        results <- list()
    else
    {
        cacheDB <- openCacheDBScope()
        annTable <- as.data.table(compounds)
        
        # UNDONE: move to general code used also by other algos
        if (calcLogP != "none")
        {
            allSMILESLogP <- union(parents$SMILES, annTable$SMILES)
            
            ph <- makeHash(allSMILESLogP, calcLogP)
            allLogP <- loadCacheData("TPsAnnCompLogP", ph, cacheDB)
            if (is.null(allLogP))
            {
                printf("Calculating log P values... ")
                allLogP <- calculateLogP(allSMILESLogP, method = calcLogP, mustWork = FALSE)
                names(allLogP) <- allSMILESLogP
                saveCacheData("TPsAnnCompLogP", allLogP, ph, cacheDB)
                printf("Done!\n")
            }
            annTable[, logP := allLogP[SMILES]]
            parents[, logP := allLogP[SMILES]]
        }
        
        parsSplit <- split(parents, seq_len(nrow(parents)))
        names(parsSplit) <- parents$name
        
        baseHash <- makeHash(compounds, TPsRef, extraOptsFMCSR, skipInvalid, prefCalcChemProps, neutralChemProps,
                             calcLogP, calcSims, fpType, fpSimMethod)
        setHash <- makeHash(parents, baseHash)
        cachedSet <- loadCacheSet("TPsAnnComp", setHash, cacheDB)
        hashes <- sapply(parsSplit, function(par) makeHash(baseHash, par[, c("name", "SMILES", "formula")],
                                                           with = FALSE))
        cachedResults <- pruneList(sapply(hashes, function(h)
        {
            result <- NULL
            if (!is.null(cachedSet))
                result <- cachedSet[[h]]
            if (is.null(result))
                result <- loadCacheData("TPsAnnComp", h, cacheDB)
            return(result)
        }, simplify = FALSE))
        
        parsTBD <- setdiff(parents$name, names(cachedResults))
        newResults <- list()
        if (length(parsTBD) > 0)
        {
            newResults <- withProg(length(parsTBD), FALSE, sapply(parsSplit[parsTBD], function(par)
            {
                nr <- getTPsCompounds(annTable, par$name, par$formula, par$SMILES, par[["logP"]], extraOptsFMCSR,
                                      if (!is.null(TPsRef)) TPsRef[[par$name]]$SMILES else NULL,
                                      fpType, fpSimMethod, minFitFormula, minFitCompound, minSimSusp, parallel)
                saveCacheData("TPsAnnComp", nr, hashes[[par$name]], cacheDB)
                return(nr)
            }, simplify = FALSE))
            newResults <- pruneList(newResults, checkZeroRows = TRUE)
        }
        
        if (is.null(cachedSet))
            saveCacheSet("TPsAnnComp", hashes, setHash, cacheDB)
        
        results <- c(cachedResults, newResults)
        results <- results[intersect(parents$name, names(results))]
        parents <- parents[name %in% names(results)]
    }
    
    return(transformationProductsAnnComp(calcSims = calcSims, fpType = fpType, fpSimMethod = fpSimMethod,
                                         parents = parents, products = results))
}
