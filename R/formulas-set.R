#' @include main.R
#' @include formulas.R
#' @include workflow-step-set.R
NULL

syncFormulasSetObjects <- function(formulasSet, makeCons)
{
    # re-generate
    formulasSet@featureFormulas <- Reduce(modifyList, lapply(formulasSet@setObjects, formulaTable, features = TRUE))
    
    if (makeCons)
    {
        groupFormsList <- sapply(formulasSet@setObjects, formulaTable, features = FALSE, simplify = FALSE)
        gNames <- groupNames(formulasSet)
        mc <- setNames(rep(length(groupFormsList), length(gNames)), gNames)
        formulasSet@formulas <- generateGroupFormulasByConsensus(groupFormsList, mc, formulasSet@setThreshold,
                                                                 formulasSet@setThresholdAnn,
                                                                 formulasSet@origFGNames, "set", "setCoverage",
                                                                 "setCoverageAnn")
    }
    else
    {
        # sync available feature groups
        allFGroups <- unique(sapply(setObjects(formulasSet), groupNames))
        formulasSet@formulas <- formulasSet@formulas[intersect(groupNames(formulasSet), allFGroups)]
        
        # only keep results from sets still present
        formulasSet@formulas <- lapply(formulasSet@formulas, function(ft) ft[set %in% sets(formulasSet)])
    }
    
    formulasSet@scoreRanges <- formulasSet@scoreRanges[groupNames(formulasSet)]
    
    return(formulasSet)
}

formulasSet <- setClass("formulasSet", slots = c(setThreshold = "numeric",
                                                 setThresholdAnn = "numeric",
                                                 origFGNames = "character"),
                        contains = c("formulas", "workflowStepSet"))

#' @describeIn formulasSet Shows summary information for this object.
#' @export
setMethod("show", "formulasSet", function(object)
{
    callAllNextMethods(object, show, firstClass = "formulas", startFrom = "formulasSet")
})

setMethod("[", c("formulasSet", "ANY", "missing", "missing"), function(x, i, j, ..., sets = NULL, drop = TRUE)
{
    assertSets(x, sets, TRUE)
    
    if (!is.null(sets))
        x@setObjects <- x@setObjects[sets]
    
    if (!missing(i))
    {
        i <- assertSubsetArgAndToChr(i, groupNames(x))
        # NOTE: assume that subsetting with non-existing i will not result in errors
        x@setObjects <- lapply(x@setObjects, "[", i = i)
        x@setObjects <- pruneList(x@setObjects, checkEmptyElements = TRUE)
    }

    if (!is.null(sets) || !missing(i))
        x <- syncFormulasSetObjects(x, FALSE)
    
    return(x)
})

setMethod("as.data.table", "formulasSet", function(x, fGroups = NULL, average = FALSE, ...)
{
    ret <- callNextMethod(x, fGroups = fGroups, average = average, ...)
    if (average)
        ret[, c("formula", "set") := NULL] # formula column doesn't make sense anymore, set column is left-over
    return(ret[])
})

#' @export
setMethod("filter", "formulasSet", function(obj, ..., negate = FALSE, sets = NULL)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFlag(negate, add = ac)
    assertSets(obj, sets, TRUE, add = ac)
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
        obj <- syncFormulasSetObjects(obj, TRUE)
        cat("Done!\n")
    }
    
    return(obj)
})

#' @export
setMethod("plotSpectrum", "formulasSet", function(obj, precursor, groupName, analysis = NULL, MSPeakLists,
                                                  title = NULL, useGGPlot2 = FALSE, mincex = 0.9, xlim = NULL,
                                                  ylim = NULL, perSet = TRUE, mirror = TRUE, ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertString(precursor, min.chars = 1, add = ac)
    checkmate::assertString(groupName, min.chars = 1, add = ac)
    checkmate::assertString(analysis, min.chars = 1, null.ok = TRUE, add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertString(title, null.ok = TRUE, add = ac)
    checkmate::assertNumber(mincex, lower = 0, finite = TRUE, add = ac)
    assertXYLim(xlim, ylim, add = ac)
    aapply(checkmate::assertFlag, . ~ useGGPlot2 + perSet + mirror, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    if (!perSet || length(sets(obj)) == 1 || !is.null(analysis))
        return(callNextMethod(obj, precursor, groupName, analysis, MSPeakLists, title,
                              useGGPlot2, mincex, xlim, ylim, ...))
    
    spec <- annotatedPeakList(obj, precursor, groupName, analysis, MSPeakLists)
    if (is.null(spec))
        return(NULL)

    if (is.null(title))
        title <- subscriptFormula(precursor)
    
    return(makeMSPlotSets(spec, title, mirror, sets(obj), mincex, xlim, ylim, useGGPlot2, ...))
})

setMethod("plotSpectrumHash", "formulasSet", function(obj, precursor, groupName, analysis = NULL, MSPeakLists,
                                                      title = NULL, useGGPlot2 = FALSE, mincex = 0.9,
                                                      xlim = NULL, ylim = NULL, perSet = TRUE, mirror = TRUE, ...)
{
    return(makeHash(callNextMethod(obj, precursor, groupName, analysis, MSPeakLists,
                                   title, useGGPlot2, mincex, xlim, ylim, ...),
                    perSet, mirror))
})


generateFormulasSet <- function(fGroupsSet, generator, ..., setArgs, setThreshold, setThresholdAnn)
{
    aapply(checkmate::assertNumber, . ~ setThreshold + setThresholdAnn, lower = 0, upper = 1, finite = TRUE)
    
    # UNDONE: mention that adduct argument is automatically set

    unsetFGroupsList <- sapply(sets(fGroupsSet), unset, obj = fGroupsSet, simplify = FALSE)
    setObjects <- Map(unsetFGroupsList, setArgs,
                      f = function(fg, sa) do.call(generator, c(list(fGroups = fg, ...), sa)))
    
    combFormulas <- Reduce(modifyList, lapply(setObjects, formulaTable, features = TRUE))
    
    groupFormsList <- sapply(setObjects, formulaTable, features = FALSE, simplify = FALSE)
    mc <- setNames(rep(length(setObjects), length(fGroupsSet)), names(fGroupsSet))
    groupForms <- generateGroupFormulasByConsensus(groupFormsList, mc, setThreshold, setThresholdAnn, names(fGroupsSet),
                                                   "set", "setCoverage", "setCoverageAnn")
    
    return(formulasSet(setObjects = setObjects, origFGNames = names(fGroupsSet), setThreshold = setThreshold,
                       setThresholdAnn = setThresholdAnn, formulas = groupForms, featureFormulas = combFormulas,
                       algorithm = makeSetAlgorithm(setObjects)))
}

formulasUnset <- setClass("formulasUnset", contains = "formulas")
setMethod("unset", "formulasSet", function(obj, set)
{
    assertSets(obj, set, FALSE)
    obj <- obj[, sets = set]
    
    groupForms <- lapply(formulaTable(obj), copy)
    groupForms <- lapply(groupForms, data.table::set, j = c("set", "setCoverage"), value = NULL)
    
    return(formulasUnset(formulas = groupForms, featureFormulas = formulaTable(obj, features = TRUE),
                         algorithm = paste0(algorithm(obj), "_unset")))
})
