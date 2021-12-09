#' @include main.R
#' @include TP.R
NULL

#' @export
transformationProductsStructure <- setClass("transformationProductsStructure", contains = "transformationProducts")

#' @templateVar class transformationProductsStructure
#' @template convertToMFDB
#' @export
setMethod("convertToMFDB", "transformationProductsStructure", function(TPs, out, includeParents = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertPathForOutput(out, overwrite = TRUE, add = ac) # NOTE: assert doesn't work on Windows...
    checkmate::assertFlag(includeParents, add = ac)
    checkmate::reportAssertions(ac)
    
    # UNDONE: update collapse function and move to doConvertToMFDB(), where it is optional and not used for TP-library
    cat("Collapsing results... ")
    prodAll <- if (length(TPs) > 0) collapseBTResults(TPs@products) else data.table()
    cat("Done!\n")
    
    doConvertToMFDB(prodAll, parents(TPs), out, includeParents)
})

setMethod("linkParentsToFGroups", "transformationProductsStructure", function(TPs, fGroups)
{
    return(screenInfo(fGroups)[name %in% names(TPs), c("name", "group"), with = FALSE])
})

#' @describeIn transformationProductsStructure Performs rule-based filtering of the \command{BioTransformer} predictions.
#'   Useful to simplify and clean-up the data.
#'
#' @param removeDuplicates If \code{TRUE} then the TPs of a parent with duplicate structures (\acronym{SMILES}) are
#'   removed. Such duplicates may occur when different transformation pathways yield the same TPs. The first TP
#'   candidate with duplicate structure will be kept.
#' @param removeParentIsomers If \code{TRUE} then TPs with an equal formula as their parent (isomers) are removed.
#' @param removeTPIsomers If \code{TRUE} then all TPs with equal formula as any sibling TPs (isomers) are removed.
#'   Unlike \code{removeDuplicates}, \emph{all} TP candidates are removed (including the first match). This filter
#'   automatically sets \code{removeDuplicates=TRUE} to avoid complete removal of TPs with equal structure.
#' @param minSimilarity Minimum structure similarity (\samp{0-1}) that a TP should have relative to its parent. For
#'   details on how these similarities are calculated, see the \code{\link{generateTPsBioTransformer}} function. May be
#'   useful under the assumption that parents and TPs who have a high structural similarity, also likely have a high
#'   MS/MS spectral similarity (which can be evaluated after componentization with \code{\link{generateComponentsTPs}}.
#' @param negate If \code{TRUE} then filters are performed in opposite manner.
#'
#' @return \code{filter} returns a filtered \code{transformationProductsStructure} object.
#'
#' @export
setMethod("filter", "transformationProductsStructure", function(obj, removeParentIsomers = FALSE,
                                                                removeTPIsomers = FALSE, removeDuplicates = FALSE,
                                                                minSimilarity = NULL, negate = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertFlag, . ~ removeParentIsomers + removeTPIsomers + removeDuplicates + negate,
           fixed = list(add = ac))
    checkmate::assertNumber(minSimilarity, lower = 0, finite = TRUE, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (length(obj) == 0)
        return(obj)
    
    if (removeTPIsomers)
        removeDuplicates <- TRUE
    
    oldn <- length(obj)
    
    hash <- makeHash(obj, removeParentIsomers, removeTPIsomers, removeDuplicates, minSimilarity, negate)
    cache <- loadCacheData("filterTPs", hash)
    if (!is.null(cache))
        obj <- cache
    else
    {
        if (removeParentIsomers || removeTPIsomers || removeDuplicates)
        {
            # NOTE: obj@products should be first arg to Map to keep names...
            obj@products <- Map(obj@products, parents(obj)$formula, f = function(prod, pform)
            {
                if (removeDuplicates)
                    prod <- if (negate) prod[duplicated(SMILES)] else prod[!duplicated(SMILES)]
                if (removeParentIsomers)
                    prod <- if (negate) prod[formula == pform] else prod[formula != pform]
                if (removeTPIsomers)
                {
                    df <- getDuplicatedStrings(prod$formula)
                    prod <- if (negate) prod[formula %chin% df] else prod[!formula %chin% df]
                }
                return(prod)
            })
        }
        
        if (!is.null(minSimilarity))
        {
            pred <- if (negate) function(x) x < minSimilarity else function(x) numGTE(x, minSimilarity)
            obj@products <- lapply(obj@products, function(p) p[pred(similarity)])
        }
        
        obj@products <- pruneList(obj@products, checkZeroRows = TRUE)
        obj@parents <- obj@parents[name %in% names(obj@products)]
        
        saveCacheData("filterTPs", obj, hash)
    }
    
    newn <- length(obj)
    printf("Done! Filtered %d (%.2f%%) TPs. Remaining: %d\n", oldn - newn, if (oldn == 0) 0 else (1-(newn/oldn))*100, newn)
    
    return(obj)
})
