#' @include main.R
#' @include formulas.R
NULL

syncFormulasSetObjects <- function(formulasSet)
{
    # re-generate
    formulasSet@featureFormulas <- Reduce(modifyList, lapply(formulasSet@setObjects, formulaTable, features = TRUE))
    groupFormsList <- sapply(formulasSet@setObjects, formulaTable, features = FALSE, simplify = FALSE)
    formulasSet@formulas <- generateGroupFormulasByConsensus(groupFormsList, formulasSet@setThreshold,
                                                             formulasSet@origFGNames, "set", "setCoverage")

    formulasSet@scoreRanges <- formulasSet@scoreRanges[groupNames(formulasSet)]
    formulasSet@adducts <- formulasSet@adducts[names(formulasSet@setObjects)]
    
    return(formulasSet)
}

formulasSet <- setClass("formulasSet", slots = c(adducts = "list", setObjects = "list",
                                                 setThreshold = "numeric",
                                                 origFGNames = "character"),
                        contains = "formulas")

setMethod("initialize", "formulasSet",
          function(.Object, ...) callNextMethod(.Object, algorithm = "set", ...))

setMethod("sets", "formulasSet", function(obj) names(obj@setObjects))
setMethod("adducts", "formulasSet", function(obj) obj@adducts)

#' @describeIn formulasSet Shows summary information for this object.
#' @export
setMethod("show", "formulasSet", function(object)
{
    callNextMethod(object)
    printf("Sets: %s\n", paste0(sets(object), collapse = ", "))
    printf("Adducts: %s\n", paste0(sapply(adducts(object), as.character), collapse = ", "))
    if (length(object@setObjects[[1]]) > 0)
        printf("Original algorithm: %s\n", algorithm(object@setObjects[[1]]))
})

setMethod("[", c("formulasSet", "ANY", "missing", "missing"), function(x, i, j, ..., sets = NULL, drop = TRUE)
{
    assertSets(x, sets)
    
    if (!is.null(sets))
        x@setObjects <- x@setObjects[sets]
    
    if (!missing(i))
    {
        # NOTE: assume that subsetting with non-existing i will not result in errors
        x@setObjects <- lapply(x@setObjects, "[", i = i)
        x@setObjects <- pruneList(x@setObjects, checkEmptyElements = TRUE)
        
        x <- syncFormulasSetObjects(x)
    }
    
    return(x)
})

setMethod("as.data.table", "formulasSet", function(x, fGroups = NULL, average = FALSE, ...)
{
    ret <- callNextMethod(x, fGroups = fGroups, average = average, ...)
    if (average)
        ret[, formula := NULL] # doesn't make a lot of sense anymore with different adducts
    return(ret[])
})

#' @export
setMethod("filter", "formulasSet", function(obj, ..., negate = FALSE, sets = NULL)
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
        obj <- syncFormulasSetObjects(obj)
        cat("Done!\n")
    }
    
    return(obj)
})

generateFormulasSet <- function(fGroupsSet, MSPeakListsSet, generator, ..., setThreshold)
{
    ionizedFGroupsList <- sapply(sets(fGroupsSet), ionize, obj = fGroupsSet, simplify = FALSE)
    ionizedMSPeaksList <- sapply(sets(MSPeakListsSet), ionize, obj = MSPeakListsSet, simplify = FALSE)
    # UNDONE: sync sets
    ionizedFormulasList <- mapply(ionizedFGroupsList, ionizedMSPeaksList, adducts(fGroupsSet),
                                  FUN = function(fg, mspl, a) generator(fGroups = fg, MSPeakLists = mspl, adduct = a, ...),
                                  SIMPLIFY = FALSE)
    
    combFormulas <- Reduce(modifyList, lapply(ionizedFormulasList, formulaTable, features = TRUE))
    
    groupFormsList <- sapply(ionizedFormulasList, formulaTable, features = FALSE, simplify = FALSE)
    groupForms <- generateGroupFormulasByConsensus(groupFormsList, setThreshold, names(fGroupsSet),
                                                   "set", "setCoverage")
    
    ret <- formulasSet(adducts = adducts(fGroupsSet), setObjects = ionizedFormulasList,
                       origFGNames = names(fGroupsSet), setThreshold = setThreshold,
                       formulas = groupForms, featureFormulas = combFormulas)
    
    return(ret)
}
