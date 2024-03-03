#' @include main.R
#' @include TP-structure.R
NULL

#' @rdname transformationProductsStructure-class
transformationProductsAnnComp <- setClass("transformationProductsAnnComp", contains = "transformationProductsStructure")

setMethod("initialize", "transformationProductsAnnComp",
          function(.Object, ...) callNextMethod(.Object, algorithm = "ann_comp", ...))

# NOTE: this function is called by a withProg() block, so handles progression updates
getTPsCompounds <- function(annTable, parentRow, TPStructParams, extraOptsFMCSR, simSuspSMILES, minRTDiff,
                            minFitFormula, minFitCompound, minSimSusp, minFitCompOrSimSusp, parallel)
{
    tab <- copy(annTable)
    setnames(tab, "neutral_formula", "formula")
    
    tab[, c("name", "ID", "parent_ID", "chem_ID", "generation") := .(paste0(parentRow$name, "-TP", .I), .I, NA_integer_, .I, 1)]
    
    # UNDONE: apply thresholds
    # UNDONE: unset annotation objects?

    tab <- getTPsRetDir(tab, parentRow, TPStructParams$calcLogP, TPStructParams$forceCalcLogP,
                        TPStructParams$minLogPDiff)

    if (!is.null(minRTDiff) && !is.null(parentRow[["rt"]]) && !is.na(parentRow$rt) && !is.null(tab[["ret"]]))
    {
        tab[, retDiff := ret - parentRow$rt]
        tab[, featRetDir := fcase((retDiff + minRTDiff) < 0, -1,
                                  (retDiff - minRTDiff) > 0, 1,
                                  default = 0)]
        tab <- tab[retDir == 0 | featRetDir == 0 | retDir == featRetDir]
    }
    
    tab[, fitFormula := sapply(formula, calcFormulaFit, parentRow$formula)]
    tab <- tab[numGTE(fitFormula, minFitFormula)]
    
    if (nrow(tab) > 0)
    {
        compFits <- do.call(calcStructFitFMCS, c(list(parentRow$SMILES, tab$SMILES, parallel), extraOptsFMCSR))
        tab[SMILES %chin% names(compFits), fitCompound := compFits[SMILES]]
        tab <- tab[numGTE(fitCompound, minFitCompound)]
        
        if (nrow(tab) >0 && !is.null(simSuspSMILES) && length(simSuspSMILES) > 0)
        {
            tab[, c("simSusp", "simSuspSMILES") := rbindlist(doApply("lapply", parallel, SMILES, function(SMI)
            {
                dists <- sapply(simSuspSMILES, patRoon:::distSMILES, SMI1 = SMI, fpType = TPStructParams$fpType,
                                fpSimMethod = TPStructParams$fpSimMethod)
                wh <- which.max(dists)
                return(list(dists[wh], simSuspSMILES[wh]))
            }, prog = FALSE))]
            tab <- tab[numGTE(NAToZero(simSusp), minSimSusp)]
        }
        else
            tab[, c("simSusp", "simSuspSMILES") := .(NA_real_, NA_character_)]

        tab <- tab[numGTE(fitCompound, minFitCompOrSimSusp[1]) | numGTE(NAToZero(simSusp), minFitCompOrSimSusp[2])]
        tab[, TPScore := pmax(fitCompound, NAToZero(simSusp)) + NAToZero(annSim)]
    }
    else
        tab[, c("fitCompound", "simSusp", "TPScore") := numeric()]
    
        
    tab <- subsetDTColumnsIfPresent(tab, c("group", "name", "ID", "parent_ID", "chem_ID", "generation", "compoundName",
                                           "SMILES", "InChI", "InChIKey", "formula", "logP", "retDir", "annSim",
                                           "fitFormula", "fitCompound", "simSusp", "simSuspSMILES", "TPScore"))
    
    doProgress()
    
    return(tab)
}


#' @export
generateTPsAnnComp <- function(parents, compounds, TPsRef = NULL, fGroupsComps = NULL, minRTDiff = 20, minFitFormula = 0,
                               minFitCompound = 0, minSimSusp = 0, minFitCompOrSimSusp = c(0, 0),
                               extraOptsFMCSR = NULL, skipInvalid = TRUE, prefCalcChemProps = TRUE,
                               neutralChemProps = FALSE, TPStructParams = getDefTPStructParams(), parallel = TRUE)
{
    # UNDONE: support >1 generations? Probably not really worthwhile...
    # UNDONE: doc that FP params are also used for suspSim
    
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
    checkmate::assertClass(fGroupsComps, "featureGroups", null.ok = TRUE, add = ac)
    checkmate::assertNumber(minRTDiff, lower = 0, finite = TRUE, null.ok = TRUE, add = add)
    aapply(checkmate::assertNumber, . ~ minFitFormula + minFitCompound + minSimSusp, lower = 0, finite = TRUE,
           fixed = list(add = ac))
    checkmate::assertNumeric(minFitCompOrSimSusp, lower = 0, finite = TRUE, len = 2, add = ac)
    checkmate::assertList(extraOptsFMCSR, null.ok = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ skipInvalid + prefCalcChemProps + neutralChemProps + parallel,
           fixed = list(add = ac))
    assertTPStructParams(TPStructParams, add = ac)
    checkmate::reportAssertions(ac)
    
    parentsTab <- getTPParents(parents, skipInvalid, prefCalcChemProps, neutralChemProps)
    
    if (nrow(parentsTab) == 0)
        results <- list()
    else
    {
        cacheDB <- openCacheDBScope()
        annTable <- as.data.table(compounds, fGroups = fGroupsComps)
        
        if (TPStructParams$calcLogP != "none")
        {
            if (!is.null(annTable[["XlogP"]]) && !is.null(annTable[["AlogP"]]))
                setnames(annTable, "XlogP", "logP") # prefer XlogP if both are available
            else
                setnames(annTable, c("XlogP", "AlogP"), c("logP", "logP"), skip_absent = TRUE)
            
            if (!is.null(fGroupsComps) && isScreening(parents))
            {
                # figure out parent RTs & groups from screening data
                gi <- groupInfo(parents)
                hasRT <- !is.null(parentsTab[["rt"]])
                parentsTab[, rt := {
                    parn <- name
                    scr <- screenInfo(parents)[name == parn]
                    if (nrow(scr) != 1)
                        if (hasRT) hasRT else NA_real_ # unclear which feature group to consider
                    else
                        gi[scr$group, "rts"]
                }, by = seq_len(nrow(parentsTab))]
            }
            
            # try to figure out log Ps from annotation data
            if (!is.null(parentsTab[["InChIKey"]]) && !is.null(annTable[["logP"]]))
            {
                parentsTab <- merge(parentsTab, annTable[, c("InChIKey", "logP"), with = FALSE], by = "InChIKey",
                                    all.x = TRUE, sort = FALSE)
            }
            
            parentsTab <- maybeCalcTPLogPs(parentsTab, TPStructParams$calcLogP, TPStructParams$forceCalcLogP)
            annTable <- maybeCalcTPLogPs(annTable, TPStructParams$calcLogP, TPStructParams$forceCalcLogP)
        }
        
        parsSplit <- split(parentsTab, seq_len(nrow(parentsTab)))
        names(parsSplit) <- parentsTab$name
        
        baseHash <- makeHash(compounds, TPsRef, fGroupsComps, minRTDiff, minFitFormula, minFitCompound, minSimSusp,
                             minFitCompOrSimSusp, extraOptsFMCSR, skipInvalid, prefCalcChemProps, neutralChemProps,
                             TPStructParams)
        setHash <- makeHash(parentsTab, baseHash)
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
        
        parsTBD <- setdiff(parentsTab$name, names(cachedResults))
        newResults <- list()
        if (length(parsTBD) > 0)
        {
            newResults <- withProg(length(parsTBD), FALSE, sapply(parsSplit[parsTBD], function(par)
            {
                sss <- if (!is.null(TPsRef) && !is.null(TPsRef[[par$name]])) TPsRef[[par$name]]$SMILES else NULL
                nr <- getTPsCompounds(annTable, par, TPStructParams, extraOptsFMCSR, sss, minRTDiff, minFitFormula,
                                      minFitCompound, minSimSusp, minFitCompOrSimSusp, parallel)
                saveCacheData("TPsAnnComp", nr, hashes[[par$name]], cacheDB)
                return(nr)
            }, simplify = FALSE))
            newResults <- pruneList(newResults, checkZeroRows = TRUE)
        }
        
        if (is.null(cachedSet))
            saveCacheSet("TPsAnnComp", hashes, setHash, cacheDB)
        
        results <- c(cachedResults, newResults)
        results <- results[intersect(parentsTab$name, names(results))]
        parentsTab <- parentsTab[name %in% names(results)]
    }
    
    # NOTE: retDirs were already calculated when gathering TPs so these could be used as a pre-filter. Hence, skip
    # retDir calculation in the constructor.
    return(transformationProductsAnnComp(doRetDirs = FALSE, TPStructParams = TPStructParams, parents = parentsTab,
                                         products = results))
}
