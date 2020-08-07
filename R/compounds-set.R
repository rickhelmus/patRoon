#' @include main.R
#' @include compounds.R
NULL

syncCompoundsSetObjects <- function(compoundsSet)
{
    # re-generate
    # UNDONE
    compoundsSet@featureFormulas <- Reduce(modifyList, lapply(compoundsSet@setObjects, formulaTable, features = TRUE))
    groupFormsList <- sapply(compoundsSet@setObjects, formulaTable, features = FALSE, simplify = FALSE)
    compoundsSet@formulas <- generateGroupFormulasByConsensus(groupFormsList, compoundsSet@setThreshold,
                                                             compoundsSet@origFGNames, "set", "setCoverage")
    
    compoundsSet@scoreRanges <- compoundsSet@scoreRanges[groupNames(compoundsSet)]
    compoundsSet@adducts <- compoundsSet@adducts[names(compoundsSet@setObjects)]
    
    return(compoundsSet)
}

compoundsSet <- setClass("compoundsSet", slots = c(adducts = "list", setObjects = "list",
                                                   setThreshold = "numeric"),
                        contains = "compounds")

setMethod("initialize", "compoundsSet",
          function(.Object, ...) callNextMethod(.Object, algorithm = "set", ...))

setMethod("sets", "compoundsSet", function(obj) names(obj@setObjects))
setMethod("adducts", "compoundsSet", function(obj) obj@adducts)

#' @describeIn compoundsSet Shows summary information for this object.
#' @export
setMethod("show", "compoundsSet", function(object)
{
    callNextMethod(object)
    printf("Sets: %s\n", paste0(sets(object), collapse = ", "))
    printf("Adducts: %s\n", paste0(sapply(adducts(object), as.character), collapse = ", "))
    if (length(object@setObjects[[1]]) > 0)
        printf("Original algorithm: %s\n", algorithm(object@setObjects[[1]]))
})

setMethod("[", c("compoundsSet", "ANY", "missing", "missing"), function(x, i, j, ..., sets = NULL, drop = TRUE)
{
    assertSets(x, sets)
    
    if (!is.null(sets))
        x@setObjects <- x@setObjects[sets]
    
    if (!missing(i))
    {
        # NOTE: assume that subsetting with non-existing i will not result in errors
        x@setObjects <- lapply(x@setObjects, "[", i = i)
        x@setObjects <- pruneList(x@setObjects, checkEmptyElements = TRUE)
        
        x <- synccompoundsSetObjects(x)
    }
    
    return(x)
})

setMethod("as.data.table", "compoundsSet", function(x, fGroups = NULL, average = FALSE, ...)
{
    ret <- callNextMethod(x, fGroups = fGroups, average = average, ...)
    if (average)
        ret[, formula := NULL] # doesn't make a lot of sense anymore with different adducts
    return(ret[])
})

#' @export
setMethod("filter", "compoundsSet", function(obj, ..., negate = FALSE, sets = NULL)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFlag(negate, add = ac)
    assertSets(obj, sets, add = ac)
    checkmate::reportAssertions(ac)
    
    if (!is.null(sets) && length(sets) > 0)
    {
        if (negate)
            sets <- setdiff(get("sets", pos = 2)(obj), sets)
        obj <- obj[, sets = sets]
    }
    
    if (...length() > 0)
    {
        # filter set objects and re-generate annotation consensus
        
        obj@setObjects <- lapply(obj@setObjects, filter, ..., negate = negate)
        obj@setObjects <- pruneList(obj@setObjects, checkEmptyElements = TRUE)
        
        # synchronize other objects
        cat("Synchronizing set objects...\n")
        obj <- synccompoundsSetObjects(obj)
        cat("Done!\n")
    }
    
    return(obj)
})

#' @export
setMethod("plotSpec", "compoundsSet", function(obj, precursor, groupName, analysis = NULL, MSPeakLists,
                                              title = NULL, useGGPlot2 = FALSE, xlim = NULL, ylim = NULL,
                                              perSet = FALSE, mirror = TRUE, ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertString(precursor, min.chars = 1, add = ac)
    checkmate::assertString(groupName, min.chars = 1, add = ac)
    checkmate::assertString(analysis, min.chars = 1, null.ok = TRUE, add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertString(title, null.ok = TRUE, add = ac)
    assertXYLim(xlim, ylim, add = ac)
    aapply(checkmate::assertFlag, . ~ useGGPlot2 + perSet, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    if (!perSet || length(sets(obj)) == 1 || !is.null(analysis))
        return(callNextMethod(obj, precursor, groupName, analysis, MSPeakLists, title,
                              useGGPlot2, xlim, ylim, ...))
    
    spec <- annotatedPeakList(obj, precursor, groupName, analysis, MSPeakLists)
    if (is.null(spec))
        return(NULL)
    
    if (is.null(title))
        title <- subscriptFormula(precursor)
    
    return(makeMSPlotSets(spec, title, mirror, sets(obj), xlim, ylim, useGGPlot2, ...))
})


generateCompoundsSet <- function(fGroupsSet, MSPeakListsSet, generator, ..., setThreshold)
{
    checkmate::assertNumber(setThreshold, lower = 0, finite = TRUE)
    msplArgs <- assertAndGetMSPLSetsArgs(fGroupsSet, MSPeakListsSet)
    
    ionizedFGroupsList <- sapply(sets(fGroupsSet), ionize, obj = fGroupsSet, simplify = FALSE)
    ionizedCompoundsList <- mapply(ionizedFGroupsList, msplArgs, adducts(fGroupsSet),
                                   FUN = function(fg, mspl, a) generator(fGroups = fg, MSPeakLists = mspl[[1]], adduct = a, ...),
                                   SIMPLIFY = FALSE)
    
    # generate consensus by...
    # - checking setThreshold
    # - merging by IK1 (or identifier?)
    # - average scores
    # - merge fragInfos and update PLIndex
    
    # get all annotated fGroup names with original order
    allAnnGNames <- intersect(names(fGroupsSet), unique(unlist(lapply(ionizedCompoundsList, groupNames))))
    
    cons <- sapply(allAnnGNames, function(gName)
    {
        allResults <- pruneList(sapply(ionizedCompoundsList, "[[", gName, simplify = FALSE))
        # allResults <- lapply(ionizedCompoundsList, function(cmp)
        # {
        #     ret <- compoundTable(cmp)
        #     ret[setdiff(allAnnGNames, groupNames(cmp))] <- list(data.table()) # add empty table for missing results
        #     return(ret[allAnnGNames]) # sync order
        # })
        
        Reduce(x = allResults, f = function(left, right)
        {
            if (nrow(left) == 0)
                return(right)
            if (nrow(right) == 0)
                return(left)
            merged <- copy(left)
            
            # UNDONE: assume score cols are same for left/right, should always be the case?
            
            # average scores for overlapping candidates
            scoreCols <- getAllCompCols(getCompScoreColNames(), names(left), NULL)
            merged[right, (scoreCols) := lapply(scoreCols,
                                                function(sc) .rowMeans(unlist(mget(c(sc, paste0("i.", sc)), inherits = TRUE)), na.rm = TRUE, m = .N, n = 2)),
                   on = "identifier"]
            
            # add missing candidates from right
            merged <- rbind(merged, right[!identifier %in% merged$identifier])
            
            # re-sort
            setorderv(merged, "score")
            
            return(merged)
        })
    }, simplify = FALSE)
    
    
browser()    

    ret <- compoundsSet(adducts = adducts(fGroupsSet), setObjects = ionizedFormulasList,
                        setThreshold = setThreshold,
                        formulas = groupForms, featureFormulas = combFormulas)
    
    return(ret)
}
