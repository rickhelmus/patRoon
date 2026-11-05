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
    else
        stop("All arguments to workflowSet constructor must be named", call. = FALSE)
    
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
setMethod("[", c("workflowSet", "ANY", "ANY", "missing"), function(x, i, j, ..., sets = NULL, reorder = FALSE,
                                                                   drop = TRUE)
{
    assertSets(x, sets, TRUE)

    if (!is.null(sets))
    {
        if (is.null(x@fGroups))
        {
            if (!is.null(x@features))
            {
                x@features <- x@features[, sets = sets, reorder = reorder]
                x@analysisInfo <- analysisInfo(x@features)
            }
            else
            {
                x@analysisInfo <- x@analysisInfo[set %in% sets]
                if (reorder)
                    x@analysisInfo <- x@analysisInfo[order(match(set, sets))]
            }
        }

        selectedSets <- if (reorder) sets else intersect(names(x@setObjects), sets)
        x@setObjects <- x@setObjects[selectedSets]
    }

    return(callNextMethod(x, i, j, ..., sets = sets, reorder = reorder))
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
    {
        if (!is.null(args[["groupAlgo"]]))
            stop("Cannot specify groupAlgo when calling makeSet() on features", call. = FALSE)
        obj@features <- do.call(makeSet, c(unname(featsSO), args))
    }
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
