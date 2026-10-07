# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include workflow.R
NULL


workflowSet <- setClass("workflowSet", slots = c(setObjects = "list"), contains = "workflow")

getAllWFSetsObject <- function(obj, sn)
{
    sapply(setObjects(obj), "[[", sn, simplify = FALSE)
}

setMethod("initialize", "workflowSet", function(.Object, ...)
{
    # args can be objects that should be passed to the parent constructor (e..g pre-made features etc) or set specific
    # analysisInfo
    args <- list(...)
    
    if (checkmate::testNames(names(args)))
    {
        for (an in names(args))
        {
            if (an != "analysisInfo" && checkmate::testDataFrame(args[[an]]))
            {
                .Object@setObjects[[an]] <- list(analysisInfo = assertAndPrepareAnaInfo(args[[an]], .var.name = sprintf("...[[%s]]", an)))
                args[[an]] <- NULL
            }
        }
    }
    
    setNames <- names(.Object@setObjects)
    if (checkmate::testNames(setNames))
    {
        if (!is.null(args[["analysisInfo"]]))
            stop("Cannot specify analysisInfo to constructor when setting set-specific analysisInfo objects",
                 call. = FALSE)
        anaInfoSets <- getAllWFSetsObject(.Object, "analysisInfo")
        args[["analysisInfo"]] <- getSetsAnaInfo(anaInfoSets)
    }
    
    do.call(callNextMethod, c(list(.Object = .Object), args))
})

#' @export
setMethod("setObjects", "workflowSet", function(obj) obj@setObjects)
    
#' @export
setMethod("sets", "workflowSet", function(obj)
{
    if (checkmate::testNames(names(setObjects(obj))))
        return(names(setObjects(obj)))
    return(character())
})

#' @export
setMethod("makeSet", "workflowSet", function(obj, ...)
{
    args <- list(...)
    if (!is.null(args[["labels"]]))
        stop("Customization of labels should be done via workflowStep()", call. = FALSE)
    args$labels <- sets(obj)
    
    checkIfAllSOPresent <- function(objs, wh)
    {
        missing <- setdiff(sets(obj), names(objs))
        if (length(objs) > 0 && length(missing) > 0)
            stop(sprintf("Non-sets %s objects are missing for sets: %s", wh, paste0(missing, collapse = ", ")),
                 call. = FALSE)
    }
    
    featsSO <- pruneList(getAllWFSetsObject(obj, "features"))
    checkIfAllSOPresent(featsSO, "features")
    fGroupsSO <- pruneList(getAllWFSetsObject(obj, "fGroups"))
    checkIfAllSOPresent(fGroupsSO, "fGroups")
    
    if (length(pruneList(fGroupsSO)) == 0)
        obj@features <- do.call(makeSet, c(unname(featsSO), args))
    else
        obj@fGroups <- do.call(makeSet, c(unname(fGroupsSO), args))
    
    return(obj)
})

#' @export
setMethod("unset", "workflowSet", function(obj, set)
{
    assertSets(obj, set, FALSE)
    
    so <- setObjects(obj, set)
    # UNDONE: make and return a workflowUnSet class for consistency?
    ret <- workflow(analysisInfo = so$analysisInfo)
    
    for (n in names(obj))
    {
        if (!is.null(slot(obj, n)))
            slot(ret, n) <- unset(slot(obj, n), set)
    }
    
    return(obj)
})
