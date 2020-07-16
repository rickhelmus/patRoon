#' @include main.R
#' @include formulas.R
NULL

syncFormulasSetObjects <- function(formulasSet, j)
{
    # NOTE: assume that subsetting with non-existing j will not result in errors
    formulasSet@setObjects <- lapply(formulasSet@setObjects, "[", j = j)
    formulasSet@setObjects <- pruneList(formulasSet@setObjects, checkEmptyElements = TRUE)
    
    # re-generate
    formulasSet@ionizedFeatureFormulas <- Reduce(modifyList, lapply(formulasSet@setObjects, formulaTable, features = TRUE))
    
    formulasSet@adducts <- formulasSet@adducts[names(formulasSet@setObjects)] # in case sets were removed
    
    # average ionized if (now) possible
    # if (allSame(adducts(formulasSet)))
    #     formulasSet@ionizedAveragedPeakLists <- do.call(averageMSPeakLists,
    #                                                        c(list(formulasSet@ionizedPeakLists,
    #                                                               formulasSet@origFGNames),
    #                                                          formulasSet@avgPeakListArgs))
    
    return(formulasSet)
}

formulasSet <- setClass("formulasSet", slots = c(adducts = "list", setObjects = "list",
                                                 ionizedFeatureFormulas = "list"),
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



generateFormulasSet <- function(fGroupsSet, MSPeakListsSet, generator, ..., setThreshold)
{
    ionizedFGroupsList <- sapply(sets(fGroupsSet), ionize, obj = fGroupsSet, simplify = FALSE)
    ionizedMSPeaksList <- sapply(sets(fGroupsSet), ionize, obj = MSPeakListsSet, simplify = FALSE)
    ionizedFormulasList <- mapply(ionizedFGroupsList, ionizedMSPeaksList, adducts(fGroupsSet),
                                  FUN = function(fg, mspl, a) generator(fGroups = fg, MSPeakLists = mspl, adduct = a, ...),
                                  SIMPLIFY = FALSE)
    
    neutralizeFTable <- function(ft)
    {
        ft <- copy(ft)
        ft[, formula := neutral_formula]
        # UNDONE: more? (check SIRIUS/DA)
        return(ft)
    }
    neutralFormulasList <- sapply(ionizedFormulasList, function(forms)
    {
        forms@featureFormulas <- lapply(forms@featureFormulas, lapply, neutralizeFTable)
        forms@formulas <- lapply(forms@formulas, neutralizeFTable)
        return(forms)
    }, simplify = FALSE)
    
    combFormulas <- Reduce(modifyList, lapply(neutralFormulasList, formulaTable, features = TRUE))
    combIonFormulas <- Reduce(modifyList, lapply(ionizedFormulasList, formulaTable, features = TRUE))
    
    groupFormsList <- sapply(neutralFormulasList, formulaTable, features = FALSE, simplify = FALSE)
    groupForms <- generateGroupFormulasByConsensus(groupFormsList, setThreshold, names(fGroupsSet),
                                                   "set", "setCoverage")
    
    ret <- formulasSet(adducts = adducts(fGroupsSet), setObjects = ionizedFormulasList,
                       ionizedFeatureFormulas = combIonFormulas, formulas = groupForms,
                       featureFormulas = combFormulas)
    
    # if (allSame(adducts(ret)))
    # {
    #     # ionize averaged combined spectra if all adducts are the same
    #     ret@ionizedAveragedPeakLists <- sapply(ret@averagedPeakLists, ionizeMSPeakList,
    #                                            adduct = adducts(fGroupsSet)[[1]], ionize = TRUE,
    #                                            simplify = FALSE)
    # }
    
    return(ret)
}
