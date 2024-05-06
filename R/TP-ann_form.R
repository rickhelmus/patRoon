#' @include main.R
#' @include TP-formula.R
NULL

#' @rdname transformationProductsFormula-class
transformationProductsAnnForm <- setClass("transformationProductsAnnForm", contains = "transformationProductsFormula")

setMethod("initialize", "transformationProductsAnnForm",
          function(.Object, ...) callNextMethod(.Object, algorithm = "ann_form", ...))

#' @export
setMethod("filter", "transformationProductsAnnForm", function(obj, ..., minFitFormula = 0, minTPScore = 0,
                                                              topMost = NULL, verbose = TRUE, negate = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertNumber(minFitFormula, lower = 0, finite = TRUE, add = ac)
    checkmate::assertCount(topMost, positive = TRUE, null.ok = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ verbose + negate, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    if (length(obj) == 0)
        return(obj)
    
    oldn <- length(obj)
    
    hash <- makeHash(obj, minFitFormula, minTPScore, topMost, negate)
    cache <- loadCacheData("filterTPs", hash)
    if (!is.null(cache))
        obj <- cache
    else
    {
        aboveThr <- function(x, t) t == 0 | is.na(x) | (negate & (x < t)) | (!negate & numGTE(x, t))
        
        obj <- delete(obj, j = function(tab, ...)
        {
            return(!(aboveThr(tab$fitFormula, minFitFormula) & aboveThr(tab$TPScore, minTPScore)))
        })
        
        if (!is.null(topMost))
        {
            obj <- delete(obj, j = function(tab, ...)
            {
                if (nrow(tab) <= topMost)
                    return(FALSE)
                ord <- order(tab$TPScore, decreasing = !negate)
                return(ord[-seq_len(topMost)])
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

setMethod("TPsFromScreening", "transformationProductsAnnForm", function(obj) FALSE)

setMethod("linkTPsToFGroups", "transformationProductsAnnForm", function(TPs, fGroups)
{
    ret <- as.data.table(TPs)[group %chin% names(fGroups), c("group", "name"), with = FALSE]
    setnames(ret, "name", "TP_name")
    return(ret)
})


# NOTE: this function is called by a withProg() block, so handles progression updates
getTPsFormulas <- function(annTable, parName, parFormula, minFitFormula)
{
    tab <- copy(annTable)
    setnames(tab, "neutral_formula", "formula")
    
    tab[, c("name", "ID", "parent_ID", "chem_ID", "generation") := .(paste0(parName, "-TP", .I), .I, NA_integer_, .I, 1)]
    
    # UNDONE: apply thresholds
    # UNDONE: unset annotation objects?
    
    tab[, fitFormula := sapply(formula, calcFormulaFit, parFormula)]
    tab[, TPScore := fitFormula + NAToZero(annSim)]
    
    tab <- subsetDTColumnsIfPresent(tab, c("group", "name", "ID", "parent_ID", "chem_ID", "generation", "formula",
                                           "annSim", "fitFormula", "TPScore"))

    tab <- tab[numGTE(fitFormula, minFitFormula)]
    
    doProgress()
    
    return(tab)
}


#' @export
generateTPsAnnForm <- function(parents, formulas, minFitFormula = 0, skipInvalid = TRUE, prefCalcChemProps = TRUE,
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
    checkmate::assertNumber(minFitFormula, add = ac)
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
            newResults <- doApply("sapply", parallel, parsSplit[parsTBD], function(par)
            {
                nr <- getTPsFormulas(annTable, par$name, par$formula, minFitFormula)
                saveCacheData("TPsAnnForm", nr, hashes[[par$name]], cacheDB)
                return(nr)
            }, simplify = FALSE)
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
