#' @include main.R
#' @include components.R
#' @include workflow-step-set.R
NULL

syncComponentsSetObjects <- function(componentsSet)
{
    # re-generate
    mcmp <- mergeComponents(componentsSet@setObjects, sets(componentsSet), "set")
    componentsSet@components <- mcmp$components
    componentsSet@componentInfo <- mcmp$componentInfo
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
    if (length(x) == 0)
        return(x)
    
    assertSets(x, sets, TRUE)
    if (!is.null(sets))
    {
        x@setObjects <- x@setObjects[sets]
        x <- syncComponentsSetObjects(x)
    }
    
    if (!missing(i) || !missing(j))
        x <- callNextMethod(x, i, j, ..., drop = drop)

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

#' @export
setMethod("delete", "componentsSet", function(obj, i = NULL, j = NULL, ...)
{
    i <- assertDeleteArgAndToChr(i, names(obj))
    
    # figure out corresponding sets and revert names of i to non-set names
    isets <- componentInfo(obj)[match(i, name)]$set
    i <- mapply(i, isets, FUN = function(x, s) sub(paste0("-", s, "$"), "", x))
    
    unisets <- unique(isets)
    obj@setObjects[unisets] <- Map(obj@setObjects[unisets], unisets,
                                   f = function(o, s) delete(o, i = i[isets == s], j = j, ...))
    
    obj <- syncComponentsSetObjects(obj)
    
    return(obj)
})

generateComponentsSet <- function(fGroupsSet, generator, setIonization, ..., setArgs = list(),
                                  classGenerator = componentsSet)
{
    annTable <- annotations(fGroupsSet)
    unsetFGroupsList <- sapply(sets(fGroupsSet), unset, obj = fGroupsSet, simplify = FALSE)
    
    if (length(setArgs) == 0)
        setArgs <- vector("list", length(unsetFGroupsList))
    
    if (setIonization)
    {
        for (i in seq_along(unsetFGroupsList))
        {
            s <- sets(fGroupsSet)[i]
            at <- annTable[set == s]
            # UNDONE: default to positive OK?
            ion <- if (nrow(at) == 0 || as.adduct(at$adduct[1])@charge >= 0) "positive" else "negative"
            setArgs[[i]]$ionization <- ion
        }
    }
        
    setObjects <- Map(unsetFGroupsList, setArgs,
                      f = function(fg, sa) do.call(generator, c(list(fGroups = fg, ...), sa)))

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

