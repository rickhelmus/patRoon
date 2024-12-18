# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

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

#' @param set \setsWF The name of the set.
#' @param sets \setsWF A \code{character} with name(s) of the sets to keep (or remove if \code{negate=TRUE}).
#'
#' @section Sets workflows: \setsWFClass{componentsSet}{components}
#'
#'   \setsWFNewMethodsSO{componentsUnset}{Only the components in the specified set are kept.}
#'
#'   \setsWFChangedMethods{
#'
#'   \item \code{filter} and the subset operator (\code{[}) Can be used to select components that are only present for
#'   selected sets.
#'
#'   }
#'
#' @rdname components-class
#' @export
componentsSet <- setClass("componentsSet", contains = c("components", "workflowStepSet"))


#' @rdname components-class
#' @export
setMethod("show", "componentsSet", function(object)
{
    callAllNextMethods(object, show, firstClass = "components", startFrom = "componentsSet")
})

#' @rdname components-class
#' @export
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

#' @rdname components-class
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
        obj <- callNextMethod(obj, ..., negate = negate)

    return(obj)
})

#' @rdname components-class
#' @export
setMethod("delete", "componentsSet", function(obj, i = NULL, j = NULL, ...)
{
    old <- obj
    obj <- callNextMethod()

    # sync setObjects
    cmpTab <- componentTable(obj); cmpTabOld <- componentTable(old)
    obj@setObjects <- Map(obj@setObjects, sets(obj), f = function(so, soName)
    {
        sNames <- paste0(names(so), "-", soName)
        removed <- !sNames %chin% names(obj)
        so <- delete(so, i = removed)
        
        if (length(so) > 0)
        {
            # sync the rest
            sNames <- sNames[!removed] # update
            so@components <- Map(so@components, sNames, f = function(socmp, sn) fintersect(obj[[sn]], socmp))
        }
        
        return(so)
    })
    
    return(obj)
})

#' @rdname components-class
#' @export
setMethod("consensus", "componentsSet", function(obj, ...)
{
    allComponents <- c(list(obj), list(...))
    
    checkmate::assertList(allComponents, types = "componentsSet", min.len = 2, any.missing = FALSE,
                          unique = TRUE, .var.name = "...")
    
    if (!allSame(lapply(allComponents, sets)))
        stop("All objects must have the same sets.")
    
    setObjects <- sapply(sets(obj), function(set)
    {
        return(do.call(consensus, lapply(lapply(allComponents, setObjects), "[[", set)))
    }, simplify = FALSE)
    
    mcmp <- mergeComponents(setObjects, names(setObjects), "set")
    
    return(componentsSet(setObjects = setObjects, components = mcmp$components, componentInfo = mcmp$componentInfo,
                         algorithm = paste0(unique(sapply(allComponents, algorithm)), collapse = ",")))
})

generateComponentsSet <- function(fGroupsSet, ionization, generator, setIonization, ..., setArgs = list(),
                                  classGenerator = componentsSet)
{
    annTable <- annotations(fGroupsSet)
    unsetFGroupsList <- sapply(sets(fGroupsSet), unset, obj = fGroupsSet, simplify = FALSE)
    verifyNoAdductIonizationArg(ionization)
    
    if (length(setArgs) == 0)
        setArgs <- vector("list", length(unsetFGroupsList))
    
    if (setIonization)
    {
        for (i in seq_along(unsetFGroupsList))
        {
            s <- sets(fGroupsSet)[i]
            at <- annTable[set == s]
            setArgs[[i]]$ionization <- getIonizationFromAnnTable(at)
        }
    }
        
    setObjects <- Map(unsetFGroupsList, setArgs,
                      f = function(fg, sa) do.call(generator, c(list(fGroups = fg, ...), sa)))

    mcmp <- mergeComponents(setObjects, sets(fGroupsSet), "set")
    
    return(classGenerator(setObjects = setObjects, components = mcmp$components, componentInfo = mcmp$componentInfo,
                          algorithm = makeSetAlgorithm(setObjects)))
}


#' @rdname components-class
#' @export
componentsUnset <- setClass("componentsUnset", contains = "components")

#' @rdname components-class
#' @export
setMethod("unset", "componentsSet", function(obj, set)
{
    assertSets(obj, set, FALSE)
    obj <- obj[, sets = set]
    cInfo <- copy(componentInfo(obj))
    cTable <- componentTable(obj)
    if (nrow(cInfo) > 0)
    {
        cInfo[, set := NULL]
        cInfo[, name := sub(paste0("\\-", set), "", name)][] # get rid of set specific names
        cTable <- setNames(cTable, cInfo$name)
    }
    return(componentsUnset(components = cTable, componentInfo = cInfo, algorithm = paste0(algorithm(obj), "_unset")))
})

