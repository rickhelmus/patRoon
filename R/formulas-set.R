#' @include main.R
#' @include formulas.R
NULL

syncFormulasSetObjects <- function(formulasSet, j)
{
    # NOTE: assume that subsetting with non-existing j will not result in errors
    formulasSet@setObjects <- lapply(formulasSet@setObjects, "[", j = j)
    formulasSet@setObjects <- pruneList(formulasSet@setObjects, checkEmptyElements = TRUE)
    formulasSet@adducts <- formulasSet@adducts[names(formulasSet@setObjects)] # in case sets were removed
    
    return(formulasSet)
}

formulasSet <- setClass("formulasSet", slots = c(adducts = "list", setObjects = "list"),
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

setMethod("[", c("formulasSet", "ANY", "missing", "missing"), function(x, i, j, ...)
{
    x <- callNextMethod(x, i, j, ...)
    
    if (!missing(i))
        x <- syncFormulasSetObjects(x, i)
    
    return(x)
})

# UNDONE: test formula averaging of as.data.table()
# UNDONE: don't average all frag errors across sets, check if averaging and normalization of scorings are OK

generateFormulasSet <- function(fGroupsSet, MSPeakListsSet, generator, ..., setThreshold)
{
    ionizedFGroupsList <- sapply(sets(fGroupsSet), ionize, obj = fGroupsSet, simplify = FALSE)
    ionizedMSPeaksList <- sapply(sets(MSPeakListsSet), ionize, obj = MSPeakListsSet, simplify = FALSE)
    # UNDONE: sync sets
    ionizedFormulasList <- mapply(ionizedFGroupsList, ionizedMSPeaksList, adducts(fGroupsSet),
                                  FUN = function(fg, mspl, a) generator(fGroups = fg, MSPeakLists = mspl, adduct = a, ...),
                                  SIMPLIFY = FALSE)
    
    neutralizeFTable <- function(ft, add)
    {
        ft <- copy(ft)
        ft[, formula := neutral_formula]
        ft[, formula_mz := formula_mz - adductMZDelta(add)]
        # UNDONE: take multiple M into account?
        ft[, frag_formula_mz := frag_formula_mz - adductMZDelta(adduct(charge = add@charge))]
        # UNDONE: more? check SIRIUS/DA. neutral_loss?
        return(ft)
    }
    neutralFormulasList <- mapply(ionizedFormulasList, adducts(fGroupsSet), FUN = function(forms, add)
    {
        forms@featureFormulas <- lapply(forms@featureFormulas, lapply, neutralizeFTable, add = add)
        forms@formulas <- lapply(forms@formulas, neutralizeFTable, add = add)
        return(forms)
    }, SIMPLIFY = FALSE)
    
    combFormulas <- Reduce(modifyList, lapply(neutralFormulasList, formulaTable, features = TRUE))
    combIonFormulas <- Reduce(modifyList, lapply(ionizedFormulasList, formulaTable, features = TRUE))
    
    groupFormsList <- sapply(neutralFormulasList, formulaTable, features = FALSE, simplify = FALSE)
    groupForms <- generateGroupFormulasByConsensus(groupFormsList, setThreshold, names(fGroupsSet),
                                                   "set", "setCoverage")
    
    ret <- formulasSet(adducts = adducts(fGroupsSet), setObjects = ionizedFormulasList,
                       formulas = groupForms, featureFormulas = combFormulas)
    
    return(ret)
}
