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

setMethod("initialize", "componentsSet",
          function(.Object, ...) callNextMethod(.Object, algorithm = "set", ...))

#' @describeIn componentsSet Shows summary information for this object.
#' @export
setMethod("show", "componentsSet", function(object)
{
    callAllNextMethods(object, show, firstClass = "components")
})

setMethod("[", c("componentsSet", "ANY", "ANY", "missing"), function(x, i, j, ..., sets = NULL, drop = TRUE)
{
    assertSets(x, sets)
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
        obj@setObjects <- lapply(obj@setObjects, filter, ..., negate = negate)
        obj@setObjects <- pruneList(obj@setObjects, checkEmptyElements = TRUE)
        
        # synchronize other objects
        cat("Synchronizing set objects...\n")
        obj <- syncComponentsSetObjects(obj)
        cat("Done!\n")
    }
    
    return(obj)
})


generateComponentsSet <- function(fGroupsSet, generator, ...)
{
    posneg <- function(add) if (add@charge < 0) "negative" else "positive"
    
    ionizedFGroupsList <- sapply(sets(fGroupsSet), ionize, obj = fGroupsSet, simplify = FALSE)
    ionizedComponentsList <- mapply(ionizedFGroupsList, adducts(fGroupsSet),
                                    FUN = function(fg, a) do.call(generator, list(fGroups = fg, ionization = posneg(a), ...)),
                                    SIMPLIFY = FALSE)
    
    mcmp <- mergeComponents(ionizedComponentsList, sets(fGroupsSet), "set")
    
    return(componentsSet(adducts = adducts(fGroupsSet), setObjects = ionizedComponentsList,
                         components = mcmp$components, componentInfo = mcmp$componentInfo))
}

