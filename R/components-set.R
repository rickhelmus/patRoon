#' @include main.R
#' @include components.R
#' @include workflow-step-set.R
NULL

syncComponentsSetObjects <- function(componentsSet)
{
    # re-generate
    if (length(componentsSet@setObjects) > 1)
    {
        mcmp <- mergeComponents(componentsSet@setObjects, sets(componentsSet), "set")
        componentsSet@components <- mcmp$components
        componentsSet@componentInfo <- mcmp$componentInfo
    }    
    return(componentsSet)
}

componentsSet <- setClass("componentsSet", slots = c(setObjects = "list"),
                          contains = c("components", "workflowStepSet"))

#' @describeIn componentsSet Shows summary information for this object.
#' @export
setMethod("show", "componentsSet", function(object)
{
    callAllNextMethods(object, show, firstClass = "components", startFrom = "componentsSet")
})

setMethod("[", c("componentsSet", "ANY", "ANY", "missing"), function(x, i, j, ..., sets = NULL, drop = TRUE)
{
    assertSets(x, sets, TRUE)
    if (!is.null(sets))
        x@setObjects <- x@setObjects[sets]

    if (!missing(i) || !missing(j))
    {
        args <- list()
        if (!missing(i))
            args <- c(args, list(i = assertSubsetArgAndToChr(i, names(x))))
        if (!missing(j))
            args <- c(args, list(j = assertSubsetArgAndToChr(j, groupNames(x))))
        
        # NOTE: assume that subsetting with non-existing i/j will not result in errors
        x@setObjects <- lapply(x@setObjects, function(o) do.call("[", args = c(list(x = o), args)))
        x@setObjects <- pruneList(x@setObjects, checkEmptyElements = TRUE)
    }
    
    if (!is.null(sets) || !missing(i) || !missing(j))
        x <- syncComponentsSetObjects(x)
    
    return(x)
})

#' @export
setMethod("filter", "componentsSet", function(obj, ..., negate = FALSE, sets = NULL)
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
        obj@setObjects <- lapply(obj@setObjects, filter, ..., negate = negate)
        obj@setObjects <- pruneList(obj@setObjects, checkEmptyElements = TRUE)
        
        # synchronize other objects
        cat("Synchronizing set objects...\n")
        obj <- syncComponentsSetObjects(obj)
        cat("Done!\n")
    }
    
    return(obj)
})


generateComponentsSet <- function(fGroupsSet, generator, ..., classGenerator = componentsSet)
{
    annTable <- annotations(fGroupsSet)
    ionizations <- sapply(sets(fGroupsSet), function(s)
    {
        at <- annTable[set == s]
        if (nrow(at) == 0)
            return("positive") # UNDONE: relevant?
        return(if (as.adduct(at$adduct[1])@charge < 0) "negative" else "positive")
    })
    
    unsetFGroupsList <- sapply(sets(fGroupsSet), unset, obj = fGroupsSet, simplify = FALSE)
    setObjects <- Map(unsetFGroupsList, ionizations,
                      f = function(fg, i) do.call(generator, list(fGroups = fg, ionization = i, ...)))
    
    mcmp <- mergeComponents(setObjects, sets(fGroupsSet), "set")
    
    return(classGenerator(setObjects = setObjects, components = mcmp$components, componentInfo = mcmp$componentInfo,
                          algorithm = makeSetAlgorithm(setObjects)))
}


componentsUnset <- setClass("componentsUnset", contains = "components")
setMethod("unset", "componentsSet", function(obj, set)
{
    assertSets(obj, set, FALSE)
    obj <- obj[, sets = set]
    return(componentsUnset(components = componentTable(obj), componentInfo = componentInfo(obj),
                           algorithm = paste0(algorithm(obj), "_unset")))
})

