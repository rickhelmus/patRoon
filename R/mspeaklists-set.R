#' @include main.R
#' @include mspeaklists.R
NULL

ionizeMSPeakList <- function(pl, adduct, ionize)
{
    adductMZ <- adductMZDelta(adduct)

    adjPL <- function(x)
    {
        x <- copy(x)
        if (ionize)
            x[, mz := mz + adductMZ]
        else
            x[, mz := mz - adductMZ]
        return(x)
    }
    
    if (!is.null(pl[["MS"]]))
        pl$MS <- adjPL(pl$MS)
    if (!is.null(pl[["MSMS"]]))
        pl$MSMS <- adjPL(pl$MSMS)
    
    return(pl)
}

#' @export
MSPeakListsSet <- setClass("MSPeakListsSet",
                           slots = c(adducts = "list", setObjects = "list", ionizedPeakLists = "list",
                                     ionizedAveragedPeakLists = "list", analysisInfo = "data.frame"),
                           contains = "MSPeakLists")

setMethod("initialize", "MSPeakListsSet",
          function(.Object, ...) callNextMethod(.Object, algorithm = "set", ...))

setMethod("sets", "MSPeakListsSet", function(obj) names(obj@setObjects))
setMethod("adducts", "MSPeakListsSet", function(obj) obj@adducts)

#' @describeIn features Get analysis information
#' @return \code{analysisInfo}: A \code{data.frame} containing a column with
#'   analysis name (\code{analysis}), its path (\code{path}), and other columns
#'   such as replicate group name (\code{group}) and blank reference
#'   (\code{blank}).
#' @export
setMethod("analysisInfo", "MSPeakListsSet", function(obj) obj@analysisInfo)

#' @describeIn MSPeakListsSet Subset on analyses/feature groups.
#' @param reAverage Set to \code{TRUE} to regenerate averaged MS peak lists
#'   after subsetting analyses.
#' @export
setMethod("[", c("MSPeakListsSet", "ANY", "ANY", "missing"), function(x, i, j, ..., reAverage = TRUE,
                                                                      sets = NULL, drop = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    assertSets(x, sets, add = ac)
    checkmate::assertFlag(reAverage, add = ac)
    checkmate::reportAssertions(ac)
    
    if (!is.null(sets) && length(sets) > 0)
        i <- mergeAnaSubsetArgWithSets(i, sets, analysisInfo(x))
    
    x <- callNextMethod(x, i, j, ..., reAverage = TRUE)

    if (!missing(i))
        x@analysisInfo <- x@analysisInfo[x@analysisInfo$analysis %in% names(peakLists(x)), ]
    
    if (!missing(i) || !missing(j))
    {
        args <- list(reAverage = reAverage)
        if (!missing(i))
            args <- c(args, list(i = analyses(x)))
        if (!missing(j))
            args <- c(args, list(j = groupNames(x)))
        
        # NOTE: assume that subsetting with non-existing i/j will not result in errors
        x@setObjects <- lapply(x@setObjects, function(o) do.call("[", args = c(list(x = o), args)))
        x@setObjects <- pruneList(x@setObjects, checkEmptyElements = TRUE)
        
        # re-generate
        x@ionizedPeakLists <- Reduce(modifyList, lapply(x@setObjects, peakLists))

        x@adducts <- x@adducts[names(x@setObjects)] # in case sets were removed
        
        # average ionized if (now) possible
        if (allSame(adducts(x)))
            x@ionizedAveragedPeakLists <- do.call(averageMSPeakLists,
                                                  c(list(x@ionizedPeakLists, x@origFGNames), x@avgPeakListArgs))
    }
    
    return(x)
})

#' @describeIn MSPeakListsSet Extract a list with MS and MS/MS (if available) peak
#'   lists. If the second argument (\code{j}) is not specified the averaged peak
#'   lists for the group specified by the first argument (\code{i}) will be
#'   returned.
#' @export
setMethod("[[", c("MSPeakListsSet", "ANY", "ANY"), function(x, i, j, neutralized = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    assertExtractArg(i, add = ac)
    if (!missing(j))
        assertExtractArg(j, add = ac)
    checkmate::assertFlag(neutralized, add = ac)
    checkmate::reportAssertions(ac)
    
    if (!neutralized && missing(j) && !allSame(adducts(x)))
        stop("No averaged ionized peak lists available because not all set adducts are equal")
    
    if (neutralized)
        return(callNextMethod(x, i, j))
    
    if (!missing(j))
        return(ionize(x)[[i, j]])
    return(ionize(x)[[i]])
})



generateMSPeakListsSet <- function(fGroupsSet, generator, ..., avgSetParams)
{
    assertAvgPListParams(avgSetParams) # UNDONE: move?
    
    # ionize all sets
    # calculate MS peak lists for each set
    # store set results in setObjects
    # store combined ionized results in ionizedPeakLists/ionizedAveragedPeakLists
    # store combined neutralized results in peakLists/averagedPeakLists
    
    ionizedFGroupsList <- sapply(sets(fGroupsSet), ionize, obj = fGroupsSet, simplify = FALSE)
    ionizedMSPeakLists <- sapply(ionizedFGroupsList, generator, ..., simplify = FALSE)
    
    neutralizedMSPL <- mapply(ionizedMSPeakLists, adducts(fGroupsSet), FUN = function(pl, add)
    {
        pl@peakLists <- lapply(pl@peakLists, lapply, ionizeMSPeakList, adduct = add, ionize = FALSE)
        pl@averagedPeakLists <- lapply(pl@averagedPeakLists, ionizeMSPeakList, adduct = add, ionize = FALSE)
        return(pl)
    }, SIMPLIFY = FALSE, USE.NAMES = TRUE)
    
    # combine (neutralized) MSPeakLists
    combPL <- Reduce(modifyList, lapply(neutralizedMSPL, peakLists))
    combPLIon <- Reduce(modifyList, lapply(ionizedMSPeakLists, peakLists))
    
    # UNDONE: set metadata?
    ret <- MSPeakListsSet(adducts = adducts(fGroupsSet), setObjects = ionizedMSPeakLists,
                          ionizedPeakLists = combPLIon, analysisInfo = analysisInfo(fGroupsSet),
                          peakLists = combPL, metadata = list(), avgPeakListArgs = avgSetParams,
                          origFGNames = names(fGroupsSet))
    
    if (allSame(adducts(ret)))
    {
        # ionize averaged combined spectra if all adducts are the same
        ret@ionizedAveragedPeakLists <- sapply(ret@averagedPeakLists, ionizeMSPeakList,
                                               adduct = adducts(fGroupsSet)[[1]], ionize = TRUE,
                                               simplify = FALSE)
    }
    
    return(ret)
}

MSPeakListsSetIonized <- setClass("MSPeakListsSetIonized", contains = "MSPeakLists")
setMethod("initialize", "MSPeakListsSetIonized",
          function(.Object, ...) callNextMethod(.Object, algorithm = "set_ionized", ...))
setMethod("ionize", "MSPeakListsSet", function(obj, sets)
{
    if (!is.null(sets) && length(sets) > 0)
        obj <- obj[, sets = sets]
    
    if (!allSame(adducts(obj)))
        stop("Selected sets for conversion must have have equal adducts")
    
    return(MSPeakListsSetIonized(peakLists = obj@ionizedPeakLists,
                                 averagedPeakLists = obj@ionizedAveragedPeakLists,
                                 metadata = list(), avgPeakListArgs = obj@avgPeakListArgs,
                                 origFGNames = obj@origFGNames))
})

