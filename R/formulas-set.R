#' @include main.R
#' @include formulas.R
NULL

syncFormulasSetObjects <- function(formulasSet)
{
    formulasSet@adducts <- formulasSet@adducts[names(formulasSet@setObjects)] # in case sets were removed
    
    # re-generate
    formulasSet@featureFormulas <- Reduce(modifyList, lapply(formulasSet@setObjects, formulaTable, features = TRUE))
    groupFormsList <- sapply(formulasSet@setObjects, formulaTable, features = FALSE, simplify = FALSE)
    formulasSet@formulas <- generateGroupFormulasByConsensus(groupFormsList, formulasSet@setThreshold,
                                                             formulasSet@origFGNames, "set", "setCoverage")
    
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
    
    # parent method is not so useful currently, as we need to re-generate
    # formulas/featureFormulas slot data
    # x <- callNextMethod(x, i, j, ...)

    if (!is.null(sets))
        x@setObjects <- x@setObjects[sets]
    
    if (!missing(i))
    {
        # from parent method
        i <- assertSubsetArgAndToChr(i, groupNames(x))
        x@scoreRanges <- x@scoreRanges[i]
        
        # NOTE: assume that subsetting with non-existing i will not result in errors
        x@setObjects <- lapply(x@setObjects, "[", i = i)
        x@setObjects <- pruneList(x@setObjects, checkEmptyElements = TRUE)
    }
    
    if (!is.null(sets) || !missing(i))
        x <- syncFormulasSetObjects(x)
    
    return(x)
})

# UNDONE: test formula averaging of as.data.table()

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
