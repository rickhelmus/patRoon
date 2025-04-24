#' @include main.R
#' @include TP-formula.R
NULL

#' Transformation products obtained from formula annotations
#'
#' Class to store results transformation products (TPs) obtained from formula annotations.
#'
#' This class is derived from the \code{\link{transformationProductsFormula}} base class, please see its documentation
#' for more details. Objects from this class are returned by \code{\link{generateTPsAnnForm}}.
#'
#' @seealso The base class \code{\link{transformationProductsFormula}} for more relevant methods and
#'   \code{\link{generateTPsAnnForm}}
#'
#' @templateVar class transformationProductsAnnForm
#' @template class-hierarchy
#'
#' @export
transformationProductsAnnForm <- setClass("transformationProductsAnnForm", contains = "transformationProductsFormula")

setMethod("initialize", "transformationProductsAnnForm",
          function(.Object, ...) callNextMethod(.Object, algorithm = "ann_form", ...))

#' @describeIn transformationProductsAnnForm Performs rule-based filtering. Useful to simplify and clean-up the data.
#'
#' @param obj The \code{transformationProductsAnnForm} object that should be filtered.
#' @param \dots Further arguments passed to the {\link[=filter,transformationProductsFormula-method]{parent filter
#'   method}}.
#' @param minFitFormula,minTPScore Thresholds related to TP scoring. See \code{\link{generateTPsAnnForm}} for more
#'   details.
#' @param topMost Only keep this number of top-most TPs (based on \code{TPScore}) for each parent/feature group
#'   combination. Set to \code{NULL} to skip this step.
#'
#' @inheritParams filter,transformationProducts-method
#'
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

setMethod("TPsFromScreening", "transformationProductsAnnForm", function(obj) FALSE)

setMethod("linkTPsToFGroups", "transformationProductsAnnForm", function(TPs, fGroups)
{
    ret <- as.data.table(TPs)[group %chin% names(fGroups), c("group", "name"), with = FALSE]
    setnames(ret, "name", "TP_name")
    return(ret)
})


# NOTE: this function is called by a withProg() block, so handles progression updates
# NOTE: parGroup may be NULL, eg when parents are from a suspect list
getTPsFormulas <- function(annTable, parName, parFormula, minFitFormula)
{
    if (nrow(annTable) == 0)
        return(data.table(group = character(), name = character(), ID = integer(), parent_ID = integer(),
                           chem_ID = integer(), generation = integer(), formula = character(), annSim = numeric(),
                           fitFormula = numeric(), TPScore = numeric()))
    
    tab <- copy(annTable)
    setnames(tab, "neutral_formula", "formula")
    
    tab[, c("name", "ID", "parent_ID", "chem_ID", "generation") := .(paste0(parName, "-TP", .I), .I, NA_integer_, .I, 1)]
    
    # UNDONE: unset annotation objects?
    
    tab[, fitFormula := sapply(formula, calcFormulaFit, parFormula)]
    tab[, TPScore := fitFormula + NAToZero(annSim)]
    
    tab <- subsetDTColumnsIfPresent(tab, c("group", "name", "ID", "parent_ID", "chem_ID", "generation", "formula",
                                           "annSim", "fitFormula", "TPScore"))

    tab <- tab[numGTE(fitFormula, minFitFormula)]
    
    doProgress()
    
    return(tab)
}


#' Obtain transformation products (TPs) from formula annotation candidates
#'
#' Transforms and prioritizes \link[=formulas]{formula annotation candidates} to obtain TPs.
#'
#' @templateVar algo formula annotations
#' @templateVar do obtain transformation products
#' @templateVar generic generateTPs
#' @templateVar algoParam ann_form
#' @template algo_generator
#'
#' @details The \code{generateTPsAnnForm} function implements the unknown TP screening from formula candidates approach
#'   as described in \insertCite{Helmus2025}{patRoon}. This algorithm does not rely on any known or predicted TPs and is
#'   therefore suitable for 'full non-target' workflows. \emph{All} formula candidates are considered as potential TPs
#'   and are ranked by the \code{TP score}: \deqn{TP score = fitFormula + annSim}
#'
#' With: \itemize{
#'
#'   \item \code{annSim}: the \link[=id-conf]{annotation similarity}
#'
#'   \item \code{fitFormula}: the common element count divided by the total element count for the formulae of the
#'   parent/TP or TP/parent (maximum is taken)
#'
#' }
#'
#'   To speed up the calculation process, several thresholds are applied to rule out unlikely candidates. These
#'   thresholds are defaulted to those derived in \insertCite{Helmus2025}{patRoon}.
#'
#'   Unlike most other TP generation algorithms, no additional suspect screening step is required.
#'
#' @param formulas The \code{\link{formulas}} object containing the formula candidates.
#' @param minFitFormula Minimum \code{fitFormula} (see Details sections) to filter out unlikely candidates.
#'
#' @templateVar req formula
#' @template tp_gen-scr
#'
#' @template parallel-arg
#'
#' @return \code{generateTPsAnnForm} returns an object of the class \code{\link{transformationProductsAnnForm}}. Please
#'   see its documentation for \emph{e.g.} filtering steps that can be performed on this object.
#'
#' @note Setting \code{parallel=TRUE} may speed up calculations, but is only favorable for long calculations due to the
#'   overhead of setting up multiple \R processes.
#'
#' @references \insertAllCited{}
#'
#' @export
generateTPsAnnForm <- function(parents, formulas, minFitFormula = 0, skipInvalid = TRUE, prefCalcChemProps = TRUE,
                               neutralChemProps = FALSE, parallel = TRUE)
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
    checkmate::assertClass(formulas, "formulas", add = ac)
    checkmate::assertNumber(minFitFormula, add = ac)
    aapply(checkmate::assertFlag, . ~ skipInvalid + prefCalcChemProps + neutralChemProps + parallel,
           fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    parents <- getTPParents(parents, skipInvalid, prefCalcChemProps, neutralChemProps, "formula")
    
    if (nrow(parents) == 0)
        results <- list()
    else
    {
        cacheDB <- openCacheDBScope()
        annTable <- as.data.table(formulas)
        
        parsSplit <- split(parents, seq_len(nrow(parents)))
        names(parsSplit) <- parents$name
        
        baseHash <- makeHash(formulas, minFitFormula)
        setHash <- makeHash(parents, baseHash)
        cachedSet <- loadCacheSet("TPsAnnForm", setHash, cacheDB)
        hashes <- sapply(parsSplit, function(par) makeHash(baseHash, par[, c("name", "formula")],
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
