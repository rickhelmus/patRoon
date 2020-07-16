#' @include main.R
#' @include formulas.R
NULL

formulasSet <- setClass("formulasSet", slots = c(adducts = "list", setObjects = "list",
                                                 ionizedFeatureFormulas = "list"),
                        contains = "formulas")

setMethod("initialize", "formulasSet",
          function(.Object, ...) callNextMethod(.Object, algorithm = "set", ...))


generateFormulasSet <- function(fGroupsSet, MSPeakListsSet, generator, ..., setThreshold)
{
    ionizedFGroupsList <- sapply(sets(fGroupsSet), ionize, obj = fGroupsSet, simplify = FALSE)
    ionizedMSPeaksList <- sapply(sets(fGroupsSet), ionize, obj = MSPeakListsSet, simplify = FALSE)
    ionizedFormulasList <- mapply(ionizedFGroupsList, ionizedMSPeaksList,
                                  FUN = function(fg, mspl) generator(fGroups = fg, MSPeakLists = mspl, ...),
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
