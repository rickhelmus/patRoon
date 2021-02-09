#' @include main.R
#' @include mspeaklists.R
NULL

syncMSPeakListsSetObjects <- function(MSPeakListsSet)
{
    # update/initialize from setObjects
    if (length(setObjects(MSPeakListsSet)) >= 2)
    {
        MSPeakListsSet@peakLists <- Reduce(modifyList, lapply(MSPeakListsSet@setObjects, peakLists))
        MSPeakListsSet@averagedPeakLists <- averageMSPeakLists(MSPeakListsSet)
    }
    else if (length(setObjects(MSPeakListsSet)) == 1)
    {
        MSPeakListsSet@peakLists <- peakLists(MSPeakListsSet@setObjects[[1]])
        MSPeakListsSet@averagedPeakLists <- averagedPeakLists(MSPeakListsSet@setObjects[[1]])
    }
    else
        MSPeakListsSet@peakLists <- MSPeakListsSet@averagedPeakLists <- list()
    
    MSPeakListsSet@analysisInfo <-
        MSPeakListsSet@analysisInfo[MSPeakListsSet@analysisInfo$analysis %in% names(peakLists(MSPeakListsSet)), ]
    
    return(MSPeakListsSet)
}

#' @export
MSPeakListsSet <- setClass("MSPeakListsSet",
                           slots = c(analysisInfo = "data.frame"),
                           contains = c("MSPeakLists", "workflowStepSet"))

setMethod("averageMSPeakLists", "MSPeakListsSet", function(obj)
{
    # create 'averaged' peak lists by simply merging the averaged lists from the setObjects
    
    cat("Merging set-averaged peak lists... ")

    hash <- makeHash(lapply(obj@setObjects, averagedPeakLists))
    avgPLists <- loadCacheData("MSPeakListsSetAvg", hash)
    
    gNames <- unique(unlist(sapply(obj@setObjects, groupNames, simplify = FALSE), use.names = FALSE))
    gNames <- intersect(obj@origFGNames, gNames) # sort to original order
    
    if (length(gNames) == 0)
        avgPLists <- list()
    else if (is.null(avgPLists))
    {
        avgPLists <- sapply(gNames, function(gName)
        {
            PLMS <- rbindlist(pruneList(lapply(obj@setObjects, function(mspl) mspl[[gName]][["MS"]])), idcol = "set")
            PLMSMS <- rbindlist(pruneList(lapply(obj@setObjects, function(mspl) mspl[[gName]][["MSMS"]])), idcol = "set")
            
            return(pruneList(list(MS = if (nrow(PLMS) > 0) PLMS else NULL,
                                  MSMS = if (nrow(PLMSMS) > 0) PLMSMS else NULL)))
        }, simplify = FALSE)
        avgPLists <- pruneList(avgPLists, checkEmptyElements = TRUE)
        saveCacheData("MSPeakListsSetAvg", avgPLists, hash)
    }

    cat("Done!\n")
    
    return(avgPLists)
})


#' @describeIn MSPeakListsSet Get analysis information
#' @return \code{analysisInfo}: A \code{data.frame} containing a column with
#'   analysis name (\code{analysis}), its path (\code{path}), and other columns
#'   such as replicate group name (\code{group}) and blank reference
#'   (\code{blank}).
#' @export
setMethod("analysisInfo", "MSPeakListsSet", function(obj) obj@analysisInfo)

#' @describeIn MSPeakListsSet Shows summary information for this object.
#' @export
setMethod("show", "MSPeakListsSet", function(object)
{
    callAllNextMethods(object, show, firstClass = "MSPeakLists", startFrom = "MSPeakListsSet")
})

#' @describeIn MSPeakListsSet Subset on analyses/feature groups.
#' @param reAverage Set to \code{TRUE} to regenerate averaged MS peak lists
#'   after subsetting analyses.
#' @export
setMethod("[", c("MSPeakListsSet", "ANY", "ANY", "missing"), function(x, i, j, ..., reAverage = FALSE,
                                                                      sets = NULL, drop = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    assertSets(x, sets, TRUE, add = ac)
    checkmate::assertFlag(reAverage, add = ac)
    checkmate::reportAssertions(ac)

    # NOTE: reAverage is ignored here, as syncMSPeakListsSetObjects() should
    # always be called and averaging for MSPeaksListSet actually only concerns
    # merging

    if (!is.null(sets) && length(sets) > 0)
        i <- mergeAnaSubsetArgWithSets(i, sets, analysisInfo(x))
    
    if (!missing(i) || !missing(j))
    {
        args <- list(reAverage = reAverage)
        if (!missing(i))
            args <- c(args, list(i = assertSubsetArgAndToChr(i, analyses(x))))
        if (!missing(j))
            args <- c(args, list(j = assertSubsetArgAndToChr(j, groupNames(x))))
        
        # NOTE: assume that subsetting with non-existing i/j will not result in errors
        x@setObjects <- lapply(x@setObjects, function(o) do.call("[", args = c(list(x = o), args)))
        x@setObjects <- pruneList(x@setObjects, checkEmptyElements = TRUE)
        
        x <- syncMSPeakListsSetObjects(x)
    }
    
    return(x)
})

#' @describeIn MSPeakListsSet Returns all MS peak list data in a table.
#'
#' @param averaged If \code{TRUE} then feature group averaged peak list data is
#'   used.
#'
#' @template as_data_table-args
#'
#' @export
setMethod("as.data.table", "MSPeakListsSet", function(x, fGroups = NULL, averaged = TRUE, sets = NULL)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroupsSet", null.ok = TRUE, add = ac)
    checkmate::assertFlag(averaged, add = ac)
    checkmate::reportAssertions(ac)
    
    ret <- callNextMethod(x, fGroups = fGroups, averaged = averaged)
    
    if (!averaged) # add set column
    {
        anaInfo <- analysisInfo(x)
        ret[, set := anaInfo[match(analysis, anaInfo$analysis), "set"]]
        setcolorder(ret, "set")
    }
    
    return(ret[])
})

#' @export
setMethod("filter", "MSPeakListsSet", function(obj, ..., negate = FALSE, sets = NULL)
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
        obj <- syncMSPeakListsSetObjects(obj)
        cat("Done!\n")
    }
    
    return(obj)
})

#' @export
setMethod("plotSpectrum", "MSPeakListsSet", function(obj, groupName, analysis = NULL, MSLevel = 1, title = NULL,
                                                     useGGPlot2 = FALSE, mincex = 0.9, xlim = NULL,
                                                     ylim = NULL, perSet = TRUE, mirror = TRUE, ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertString(groupName, min.chars = 1, add = ac)
    checkmate::assertString(analysis, min.chars = 1, null.ok = TRUE, add = ac)
    checkmate::assertChoice(MSLevel, 1:2, add = ac)
    checkmate::assertNumber(mincex, lower = 0, finite = TRUE, add = ac)
    assertXYLim(xlim, ylim, add = ac)
    aapply(checkmate::assertFlag, . ~ useGGPlot2 + perSet, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    if (!perSet || length(sets(obj)) == 1 || !is.null(analysis))
        return(callNextMethod(obj, groupName, analysis, MSLevel, title, useGGPlot2, mincex, xlim, ylim, ...))
    
    spec <- getSpec(obj, groupName, MSLevel, NULL)
    if (is.null(spec))
        return(NULL)

    if (is.null(title))
        title <- getMSPeakListPlotTitle(MSLevel, analysis, groupName)
    
    return(makeMSPlotSets(spec, title, mirror, sets(obj), mincex, xlim, ylim, useGGPlot2, ...))
})

generateMSPeakListsSet <- function(fGroupsSet, generator, ...)
{
    # unset all fGroups sets, calculate MS peak lists for each set and store in setObjects
    # store combined setObject results in peakLists
    # store merged averaged peak lists in averagedPeakLists
    
    unsetFGroupsList <- sapply(sets(fGroupsSet), unset, obj = fGroupsSet, simplify = FALSE)
    setObjects <- sapply(unsetFGroupsList, generator, ..., simplify = FALSE)
    
    # combine non averaged (per analysis) MSPeakLists
    combPL <- Reduce(modifyList, lapply(setObjects, peakLists))

    # UNDONE: set metadata?
    ret <- MSPeakListsSet(setObjects = setObjects,
                          analysisInfo = analysisInfo(fGroupsSet),
                          peakLists = combPL, metadata = list(),
                          origFGNames = names(fGroupsSet),
                          algorithm = makeSetAlgorithm(setObjects))
    
    return(ret)
}

MSPeakListsUnset <- setClass("MSPeakListsUnset", contains = "MSPeakLists")
setMethod("unset", "MSPeakListsSet", function(obj, set)
{
    assertSets(obj, set, FALSE)
    obj <- obj[, sets = set]
    
    avArgs <- if (length(obj@setObjects) > 0) obj@setObjects[[1]]@avgPeakListArgs else list()
    
    # # only re-average if >1 sets, otherwise just copy the setObject
    # UNDONE: re-enable if unset will support >1 sets again...
    # if (length(obj@setObjects) > 1)
    #     return(MSPeakListsUnset(peakLists = obj@peakLists, metadata = list(), avgPeakListArgs = avArgs,
    #                             origFGNames = obj@origFGNames, algorithm = paste0(algorithm(obj), "_unset")))
    
    return(MSPeakListsUnset(peakLists = obj@setObjects[[1]]@peakLists,
                            averagedPeakLists = averagedPeakLists(obj@setObjects[[1]]),
                            metadata = list(), avgPeakListArgs = avArgs,
                            origFGNames = obj@origFGNames, algorithm = paste0(algorithm(obj), "_unset")))
})

