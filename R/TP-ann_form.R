#' @include main.R
#' @include TP-formula.R
NULL

#' @rdname transformationProductsFormula-class
transformationProductsAnnForm <- setClass("transformationProductsAnnForm", contains = "transformationProductsFormula")

setMethod("initialize", "transformationProductsAnnForm",
          function(.Object, ...) callNextMethod(.Object, algorithm = "ann_form", ...))

# NOTE: this function is called by a withProg() block, so handles progression updates
# UNDONE: does parallel make sense here?
getTPsFormulas <- function(annTable, parName, parFormula, parallel)
{
    tab <- copy(annTable)
    setnames(tab, "neutral_formula", "formula")
    
    tab[, c("name", "ID", "parent_ID", "chem_ID", "generation") := .(paste0(parName, "-TP", .I), .I, NA_integer_, .I, 1)]
    
    # UNDONE: apply thresholds
    # UNDONE: unset annotation objects?
    
    tab[, formulaDiff := sapply(formula, getFormulaDiffText, form2 = parFormula)]
    tab[, fitFormula := sapply(formula, calcFormulaFit, parFormula)]

    # UNDONE: don't do here? Maybe just in components and for all TPs types    
    # if (!is.null(parentAnn) && nrow(parentAnn) == 1)
    # {
    #     tab[, c("fragMatches", "NLMatches") := {
    #         fi <- fragInfo[[1]]
    #         if (nrow(fi) > 0)
    #             list(sum(fi$ion_formula %chin% parentAnn$fragInfo[[1]]$ion_formula),
    #                  sum(fi$neutral_loss %chin% parentAnn$fragInfo[[1]]$neutral_loss))
    #         else
    #             list(NA_integer_, NA_integer_)
    #     }, by = seq_len(nrow(tab))]
    # }
    
    tab[, TP_score := NAToZero(fitFormula) + NAToZero(annSim)]
    
    tab <- subsetDTColumnsIfPresent(tab, c("group", "ID", "parent_ID", "chem_ID", "generation", "formula", "annSim",
                                           "fitFormula", "TP_score"))
    
    doProgress()
    
    return(tab)
}


#' @export
generateTPsAnnForm <- function(parents, formulas, skipInvalid = TRUE, prefCalcChemProps = TRUE,
                               neutralChemProps = FALSE, parallel = TRUE)
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
    checkmate::assertClass(formulas, "formulas", add = ac)
    aapply(checkmate::assertFlag, . ~ skipInvalid + prefCalcChemProps + neutralChemProps + parallel,
           fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    parents <- getTPParents(parents, skipInvalid, prefCalcChemProps, neutralChemProps)
    
    if (nrow(parents) == 0)
        results <- list()
    else
    {
        cacheDB <- openCacheDBScope()
        annTable <- as.data.table(formulas)
        
        parsSplit <- split(parents, seq_len(nrow(parents)))
        names(parsSplit) <- parents$name
        
        baseHash <- makeHash(formulas, skipInvalid, prefCalcChemProps, neutralChemProps)
        setHash <- makeHash(parents, baseHash)
        cachedSet <- loadCacheSet("TPsAnnForm", setHash, cacheDB)
        hashes <- sapply(parsSplit, function(par) makeHash(baseHash, par[, c("name", "SMILES", "formula")],
                                                           with = FALSE))
        cachedResults <- pruneList(sapply(hashes, function(h)
        {
            result <- NULL
            if (!is.null(cachedSet))
                result <- cachedSet[[h]]
            if (is.null(result))
                result <- loadCacheData("TPsAnnForm", h, cacheDB)
            return(result)
        }, simplify = FALSE))
        
        parsTBD <- setdiff(parents$name, names(cachedResults))
        newResults <- list()
        if (length(parsTBD) > 0)
        {
            newResults <- withProg(length(parsTBD), FALSE, sapply(parsSplit[parsTBD], function(par)
            {
                nr <- getTPsFormulas(annTable, par$name, par$formula, parallel)
                saveCacheData("TPsAnnForm", nr, hashes[[par$name]], cacheDB)
                return(nr)
            }, simplify = FALSE))
            newResults <- pruneList(newResults, checkZeroRows = TRUE)
        }
        
        if (is.null(cachedSet))
            saveCacheSet("TPsAnnForm", hashes, setHash, cacheDB)
        
        results <- c(cachedResults, newResults)
        results <- results[intersect(parents$name, names(results))]
        parents <- parents[name %in% names(results)]
    }
    
    return(transformationProductsAnnForm(parents = parents, products = results))
}
