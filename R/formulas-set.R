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
        i <- assertSubsetArgAndToChr(i, groupNames(x))
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

#' @export
setMethod("plotSpec", "formulasSet", function(obj, precursor, groupName, analysis = NULL, MSPeakLists,
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
    aapply(checkmate::assertFlag, . ~ useGGPlot2 + perSet + mirror, fixed = list(add = ac))
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

setMethod("plotSpecHash", "formulasSet", function(obj, precursor, groupName, analysis = NULL, MSPeakLists,
                                                  title = NULL, useGGPlot2 = FALSE, xlim = NULL, ylim = NULL,
                                                  perSet = FALSE, mirror = TRUE, ...)
{
    return(makeHash(callNextMethod(obj, precursor, groupName, analysis, MSPeakLists,
                                   title, useGGPlot2, xlim, ylim),
                    perSet, mirror))
})


generateFormulasSet <- function(fGroupsSet, generator, ..., setArgs, setThreshold)
{
    checkmate::assertNumber(setThreshold, lower = 0, finite = TRUE)
    
    # UNDONE: mention that adduct argument is automatically set

    ionizedFGroupsList <- sapply(sets(fGroupsSet), ionize, obj = fGroupsSet, simplify = FALSE)
    ionizedFormulasList <- mapply(ionizedFGroupsList, adducts(fGroupsSet), setArgs,
                                  FUN = function(fg, a, sa) do.call(generator, c(list(fGroups = fg, adduct = a, ...), sa)),
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

formulasSetIonized <- setClass("formulasSetIonized", contains = "formulas")
setMethod("initialize", "formulasSetIonized",
          function(.Object, ...) callNextMethod(.Object, algorithm = "set_ionized", ...))
setMethod("ionize", "formulasSetIonized", function(obj, sets)
{
    if (!is.null(sets) && length(sets) > 0)
        obj <- obj[, sets = sets]
    
    assertEqualAdducts(adducts(obj))
    
    groupForms <- copy(formulaTable(obj))
    groupForms <- lapply(groupForms, set, j = c("set", "setCoverage"), value = NULL)
    
    return(formulasSetIonized(formulas = groupForms, featureFormulas = formulaTable(obj, features = TRUE)))
})
