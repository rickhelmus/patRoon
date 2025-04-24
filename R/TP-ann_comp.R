#' @include main.R
#' @include TP-structure.R
NULL

#' Transformation products obtained from compound annotations
#'
#' Class to store results transformation products (TPs) obtained from compound annotations.
#'
#' This class is derived from the \code{\link{transformationProductsStructure}} base class, please see its documentation
#' for more details. Objects from this class are returned by \code{\link{generateTPsAnnComp}}.
#'
#' @seealso The base class \code{\link{transformationProductsStructure}} for more relevant methods and
#'   \code{\link{generateTPsAnnComp}}
#'
#' @templateVar class transformationProductsAnnComp
#' @template class-hierarchy
#'
#' @export
transformationProductsAnnComp <- setClass("transformationProductsAnnComp", contains = "transformationProductsStructure")

setMethod("initialize", "transformationProductsAnnComp",
          function(.Object, ...) callNextMethod(.Object, algorithm = "ann_comp", ...))

#' @describeIn transformationProductsAnnComp Performs rule-based filtering. Useful to simplify and clean-up the data.
#'
#' @param obj The \code{transformationProductsAnnComp} object that should be filtered.
#' @param \dots Further arguments passed to the {\link[=filter,transformationProductsStructure-method]{parent filter
#'   method}}.
#' @param minFitFormula,minFitCompound,minSimSusp,minFitCompOrSimSusp,minTPScore Thresholds related to TP scoring. See
#'   \code{\link{generateTPsAnnComp}} for more details.
#' @param topMost Only keep this number of top-most TPs (based on \code{TPScore}) for each parent/feature group
#'   combination. Set to \code{NULL} to skip this step.
#' 
#' @inheritParams filter,transformationProducts-method
#'
#' @export
setMethod("filter", "transformationProductsAnnComp", function(obj, ..., minFitFormula = 0, minFitCompound = 0,
                                                              minSimSusp = 0, minFitCompOrSimSusp = c(0, 0),
                                                              minTPScore = 0, topMost = NULL, verbose = TRUE,
                                                              negate = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertNumber, . ~ minFitFormula + minFitCompound + minSimSusp, lower = 0, finite = TRUE,
           fixed = list(add = ac))
    checkmate::assertNumeric(minFitCompOrSimSusp, lower = 0, finite = TRUE, len = 2, add = ac)
    checkmate::assertCount(topMost, positive = TRUE, null.ok = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ verbose + negate, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    if (length(obj) == 0)
        return(obj)
    
    oldn <- length(obj)
    
    hash <- makeHash(obj, minFitFormula, minFitCompound, minSimSusp, minFitCompOrSimSusp, minTPScore, topMost, negate)
    cache <- loadCacheData("filterTPs", hash)
    if (!is.null(cache))
        obj <- cache
    else
    {
        aboveThr <- function(x, t) t == 0 | is.na(x) | (negate & (x < t)) | (!negate & numGTE(x, t))

        obj <- delete(obj, j = function(tab, ...)
        {
            return(!(aboveThr(tab$fitFormula, minFitFormula) &
                         aboveThr(tab$fitCompound, minFitCompound) &
                         aboveThr(tab$simSusp, minSimSusp) &
                         (aboveThr(tab$fitCompound, minFitCompOrSimSusp[1]) | aboveThr(tab$simSusp, minFitCompOrSimSusp[2])) &
                         aboveThr(tab$TPScore, minTPScore)))
        })
        
        if (!is.null(topMost))
        {
            obj <- delete(obj, j = function(tab, ...)
            {
                if (nrow(tab) <= topMost)
                    return(FALSE)
                tab <- copy(tab)
                tab[, ord := data.table::frank(-TPScore, ties.method = "first"), by = "group"]
                return(if (negate) tab[ord <= topMost, which = TRUE] else tab[ord > topMost, which = TRUE])
            })
        }

        saveCacheData("filterTPs", obj, hash)
    }
    
    if (...length() > 0)
        obj <- callNextMethod(obj, ..., verbose = FALSE, negate = negate)
    
    if (verbose)
    {
        newn <- length(obj)
        printf("Done! Filtered %d (%.2f%%) TPs. Remaining: %d\n", oldn - newn, if (oldn == 0) 0 else (1-(newn/oldn))*100, newn)
    }
    
    return(obj)
})

setMethod("TPsFromScreening", "transformationProductsAnnComp", function(obj) FALSE)

setMethod("linkTPsToFGroups", "transformationProductsAnnComp", function(TPs, fGroups)
{
    ret <- as.data.table(TPs)[group %chin% names(fGroups), c("group", "name"), with = FALSE]
    setnames(ret, "name", "TP_name")
    return(ret)
})


# NOTE: this function is called by a withProg() block, so handles progression updates
getTPsCompounds <- function(annTable, parentRow, TPStructParams, extraOptsFMCSR, simSusps, minRTDiff,
                            minFitFormula, minFitCompound, minSimSusp, minFitCompOrSimSusp, neutralizeTPs, parallel)
{
    if (nrow(annTable) == 0)
        return(data.table(group = character(), name = character(), ID = integer(), parent_ID = integer(),
                          chem_ID = integer(), generation = integer(), compoundName = character(), SMILES = character(),
                          InChI = character(), InChIKey = character(), formula = character(), annSim = numeric(),
                          fitFormula = numeric(), fitCompound = numeric(), simSusp = numeric(),
                          simSuspSMILES = numeric(), "simSuspInChI", simSuspInChIKey = numeric(), TPScore = numeric()))
    
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
        
        if (nrow(tab) > 0 && !is.null(simSusps) && nrow(simSusps) > 0)
        {
            simSusps <- copy(simSusps)
            setnames(simSusps, paste0("simSusp", names(simSusps)))
            tab[, c("simSusp", names(simSusps)) := rbindlist(doApply("lapply", parallel, SMILES, function(SMI)
            {
                dists <- sapply(simSusps$simSuspSMILES, patRoon:::distSMILES, SMI1 = SMI, fpType = TPStructParams$fpType,
                                fpSimMethod = TPStructParams$fpSimMethod)
                wh <- which.max(dists)
                return(c(list(dists[wh]), as.list(simSusps[wh])))
            }, prog = FALSE))]
            tab <- tab[numGTE(NAToZero(simSusp), minSimSusp)]
        }
        else
            tab[, simSusp := NA_real_]

        tab <- tab[numGTE(fitCompound, minFitCompOrSimSusp[1]) | numGTE(NAToZero(simSusp), minFitCompOrSimSusp[2])]
        tab[, TPScore := pmax(fitCompound, NAToZero(simSusp)) + NAToZero(annSim)]
    }
    else
        tab[, c("fitCompound", "simSusp", "TPScore") := numeric()]
    
    # NOTE: we also need to do this for eg SIRIUS which doesn't report InChIKeys
    tab <- prepareChemTable(tab, prefCalcChemProps = FALSE, neutralChemProps = neutralizeTPs, verbose = FALSE)
    
    tab <- subsetDTColumnsIfPresent(tab, c("group", "name", "ID", "parent_ID", "chem_ID", "generation", "compoundName",
                                           "SMILES", "InChI", "InChIKey", "formula", "logP", "retDir", "annSim",
                                           "fitFormula", "fitCompound", "simSusp", "simSuspSMILES", "simSuspInChI",
                                           "simSuspInChIKey", "TPScore"))
    
    doProgress()
    
    return(tab)
}


#' Obtain transformation products (TPs) from compound annotation candidates
#'
#' Transforms and prioritizes \link[=compounds]{compound annotation candidates} to obtain TPs.
#'
#' @templateVar algo compound annotations
#' @templateVar do obtain transformation products
#' @templateVar generic generateTPs
#' @templateVar algoParam ann_comp
#' @template algo_generator
#'
#' @details The \code{generateTPsAnnComp} function implements the unknown TP screening from compound candidates approach
#'   as described in \insertCite{Helmus2025}{patRoon}. This algorithm does not rely on any known or predicted TPs and is
#'   therefore suitable for 'full non-target' workflows. \emph{All} compound candidates are considered as potential TPs
#'   and are ranked by the \code{TP score}:
#'
#'   \deqn{TP score = max(fitCompound,simSusp) + annSim}
#'
#' With: \itemize{
#'
#'   \item \code{annSim}: the \link[=id-conf]{annotation similarity}
#'
#'   \item \code{fitCompound}: the structural fit of the compound candidate into the parent (or vice versa, maximum is
#'   taken). Calculated as the \code{"Overlap coefficient"} with \code{\link[fmcsR:fmcs]{fmcsR::fmcs}}.
#'
#'   \item \code{simSusp}: the maximum structural similarity with TP suspect candidates for this parent, \emph{i.e.}
#'   obtained from other algorithms of \code{\link{generateTPs}}). The calculation is configured by the
#'   \link[=getDefTPStructParams]{TPStructParams}.
#'
#' }
#'
#'   To speed up the calculation process, several thresholds are applied to rule out unlikely candidates. These
#'   thresholds are defaulted to those derived in \insertCite{Helmus2025}{patRoon}. Nevertheless, calculations can take
#'   a very long time (multiple hours), especially when processing large numbers of candidates from \emph{e.g.} PubChem.
#'
#'   Unlike most other TP generation algorithms, no additional suspect screening step is required.
#'
#' @param TPsRef A \code{\link{transformationProductsStructure}} object containing suspect TP candidates obtained for
#'   the same \code{parents} from another TP generation algorithm (\emph{e.g.} \code{\link{generateTPsBioTransformer}}).
#'   This is used for the calculation of \code{simSusps} (see Details). Set to \code{NULL} to skip its calculation.
#' @param compounds The \code{\link{compounds}} object containing the compound candidates.
#' @param fGroupsComps The \code{\link{featureGroups}} object for which the \code{compounds} were generated. This is
#'   used to obtain retention times for the calculation for \link[=retDir]{retention order directions}. Set to
#'   \code{NULL} to skip its calculation.
#' @param minRTDiff Minimum retention time (in seconds) difference between the parent and a TP to calculate the
#'   \link[=retDir]{retention order direction}. Candidates with unexpected retention orders are filtered out.
#' @param minFitFormula,minFitCompound,minSimSusp Thresholds to filter out unlikely candidates. For \code{fitFormula}:
#'   see \code{\link{generateTPsAnnForm}}, for the others see the Details section.
#' @param minFitCompOrSimSusp A two-sized numeric vector specifying the thresholds for \code{fitCompound} \emph{or}
#'   \code{simSusp}, respectively.
#' @param extraOptsFMCSR A \code{list} with additional options passed to the \code{\link[fmcsR:fmcs]{fmcsR::fmcs}}
#'   function. The following defaults are set: \code{au=1, bu=4, matching.mode="aromatic"}.
#'
#' @templateVar req \acronym{SMILES} or \acronym{InChI}
#' @templateVar neutTPs TRUE
#' @template tp_gen-scr
#' 
#' @template parallel-arg
#' @template tp_gen-struct_params
#'
#' @return \code{generateTPsAnnComp} returns an object of the class \code{\link{transformationProductsAnnComp}}. Please
#'   see its documentation for \emph{e.g.} filtering steps that can be performed on this object.
#'
#' @note Setting \code{parallel=TRUE} can speed up calculations considerably on multi-core systems. but will also add to RAM
#'   usage. Furthermore, parallelization is only favorable for long calculations due to the overhead of setting up
#'   multiple \R processes.
#'
#'   It is possible that candidates are equal to their parent. To remove these the \code{removeParentIsomers}
#'   \code{\link[=filter,transformationProductsStructure-method]{filter}} can be used afterwards.
#'
#' @references \insertAllCited{} \cr\cr \insertRef{Wang2013}{patRoon}
#'
#' @export
generateTPsAnnComp <- function(parents, compounds, TPsRef = NULL, fGroupsComps = NULL, minRTDiff = 20, minFitFormula = 0,
                               minFitCompound = 0, minSimSusp = 0, minFitCompOrSimSusp = c(0, 0),
                               extraOptsFMCSR = NULL, skipInvalid = TRUE, prefCalcChemProps = TRUE,
                               neutralChemProps = FALSE, neutralizeTPs = TRUE, TPStructParams = getDefTPStructParams(),
                               parallel = TRUE)
{
    # UNDONE: support >1 generations? Probably not really worthwhile...
    
    checkmate::assert(
        checkmate::checkClass(parents, "data.frame"),
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
    aapply(checkmate::assertFlag, . ~ skipInvalid + prefCalcChemProps + neutralChemProps + neutralizeTPs + parallel,
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
                        gi[group == scr$group]$ret
                }, by = seq_len(nrow(parentsTab))]
            }
            
            # try to figure out log Ps from annotation data
            if (!is.null(parentsTab[["InChIKey"]]) && !is.null(annTable[["logP"]]))
            {
                parentsTab[, IK1 := getIKBlock1(parentsTab$InChIKey)]
                parentsTab <- merge(parentsTab, annTable[, c("UID", "logP"), with = FALSE], by.x = "IK1", by.y = "UID",
                                    all.x = TRUE, sort = FALSE)
                parentsTab[, IK1 := NULL]
            }
            
            parentsTab <- maybeCalcTPLogPs(parentsTab, TPStructParams$calcLogP, TPStructParams$forceCalcLogP)
            annTable <- maybeCalcTPLogPs(annTable, TPStructParams$calcLogP, TPStructParams$forceCalcLogP)
        }
        
        parsSplit <- split(parentsTab, seq_len(nrow(parentsTab)))
        names(parsSplit) <- parentsTab$name
        
        baseHash <- makeHash(compounds, TPsRef, fGroupsComps, minRTDiff, minFitFormula, minFitCompound, minSimSusp,
                             minFitCompOrSimSusp, extraOptsFMCSR, neutralizeTPs, TPStructParams)
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
                simSusps <- if (!is.null(TPsRef) && !is.null(TPsRef[[par$name]]))
                    subsetDTColumnsIfPresent(TPsRef[[par$name]], c("SMILES", "InChI", "InChIKey"))
                else
                    NULL
                nr <- getTPsCompounds(annTable, par, TPStructParams, extraOptsFMCSR, simSusps, minRTDiff, minFitFormula,
                                      minFitCompound, minSimSusp, minFitCompOrSimSusp, neutralizeTPs, parallel)
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
