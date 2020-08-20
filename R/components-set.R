#' @include main.R
#' @include components.R
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
                          contains = "components")

setMethod("initialize", "componentsSet",
          function(.Object, ...) callNextMethod(.Object, algorithm = "set", ...))

setMethod("sets", "componentsSet", function(obj) names(obj@setObjects))

#' @describeIn componentsSet Shows summary information for this object.
#' @export
setMethod("show", "componentsSet", function(object)
{
    callNextMethod(object)
    printf("Sets: %s\n", paste0(sets(object), collapse = ", "))
    if (length(object@setObjects[[1]]) > 0)
        printf("Original algorithm: %s\n", algorithm(object@setObjects[[1]]))
})

setMethod("[", c("componentsSet", "ANY", "ANY", "missing"), function(x, i, j, ..., sets = NULL)
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

#' @export
setMethod("plotSpec", "componentsSet", function(obj, index, groupName, MSPeakLists, formulas = NULL,
                                               plotStruct = TRUE, title = NULL, useGGPlot2 = FALSE, xlim = NULL,
                                               ylim = NULL, perSet = FALSE, mirror = TRUE, ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertCount(index, positive = TRUE, add = ac)
    checkmate::assertString(groupName, min.chars = 1, add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertClass(formulas, "formulas", null.ok = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ plotStruct + useGGPlot2 + perSet + mirror, fixed = list(add = ac))
    assertXYLim(xlim, ylim, add = ac)
    checkmate::reportAssertions(ac)
    
    if (!perSet || length(sets(obj)) == 1)
        return(callNextMethod(obj, index, groupName, MSPeakLists, formulas,
                              plotStruct, title, useGGPlot2, xlim, ylim, ...))
    
    spec <- annotatedPeakList(obj, index, groupName, MSPeakLists, formulas)
    if (is.null(spec))
        return(NULL)
    
    compr <- obj[[groupName]][index, ]    
    mol <- NULL
    if (plotStruct)
    {
        mol <- getMoleculesFromSMILES(compr$SMILES)
        if (!isValidMol(mol))
            mol <- NULL
    }
    
    if (is.null(title))
        title <- getComponentsSpecPlotTitle(compr$compoundName, compr$formula)
    
    return(makeMSPlotSets(spec, title, mirror, sets(obj), xlim, ylim, useGGPlot2, ..., mol = mol))
})

setMethod("plotSpecHash", "componentsSet", function(obj, index, groupName, MSPeakLists, formulas = NULL,
                                                   plotStruct = TRUE, title = NULL, useGGPlot2 = FALSE, xlim = NULL,
                                                   ylim = NULL, perSet = FALSE, mirror = TRUE, ...)
{
    return(makeHash(callNextMethod(obj, index, groupName, MSPeakLists, formulas,
                                   plotStruct, title, useGGPlot2, xlim,
                                   ylim, ...),
                    perSet, mirror))
})

generateComponentsSet <- function(fGroupsSet, generator, ..., setArgs)
{
    ionizedFGroupsList <- sapply(sets(fGroupsSet), ionize, obj = fGroupsSet, simplify = FALSE)
    ionizedComponentsList <- mapply(ionizedFGroupsList, setArgs,
                                    FUN = function(fg, sa) do.call(generator, c(list(fGroups = fg, ...), sa)),
                                    SIMPLIFY = FALSE)
    
    mcmp <- mergeComponents(ionizedComponentsList, sets(fGroupsSet), "set")
    
    return(componentsSet(setObjects = ionizedComponentsList, components = mcmp$components, componentInfo = mcmp$componentInfo))
}

