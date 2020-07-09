#' @include main.R
#' @include mspeaklists.R
NULL

neutralizeMSPeakLists <- function(MSPeakLists, adduct)
{
    adductMZ <- adductMZDelta(adduct)
    
    neutralizePL <- function(pl)
    {
        if (!is.null(pl[["MS"]]))
        {
            pl$MS <- copy(pl$MS)
            pl$MS[, mz := mz - adductMZ]
        }
        if (!is.null(pl[["MSMS"]]))
        {
            pl$MSMS <- copy(pl$MSMS)
            pl$MSMS[, mz := mz - adductMZ]
        }
        return(pl)
    }
    
    MSPeakLists@peakLists <- lapply(MSPeakLists@peakLists, lapply, neutralizePL)
    MSPeakLists@averagedPeakLists <- lapply(MSPeakLists@averagedPeakLists, neutralizePL)
    
    return(MSPeakLists)
}

#' @export
MSPeakListsSet <- setClass("MSPeakListsSet",
                           slots = c(setObjects = "list", ionizedPeakLists = "list",
                                     analysisInfo = "data.frame"),
                           contains = "MSPeakLists")

setMethod("initialize", "MSPeakListsSet",
          function(.Object, ...) callNextMethod(.Object, algorithm = "set", ...))

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
    }
    
    return(x)
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
    
    neutralizedMSPL <- mapply(ionizedMSPeakLists, adducts(fGroupsSet), FUN = neutralizeMSPeakLists,
                              SIMPLIFY = FALSE, USE.NAMES = TRUE)
    
    # combine (neutralized) MSPeakLists
    combPL <- Reduce(modifyList, lapply(neutralizedMSPL, peakLists))
    combPLIon <- Reduce(modifyList, lapply(ionizedMSPeakLists, peakLists))
    
    # UNDONE: set metadata?
    return(MSPeakListsSet(setObjects = ionizedMSPeakLists, ionizedPeakLists = combPLIon,
                          analysisInfo = analysisInfo(fGroupsSet),
                          peakLists = combPL, metadata = list(),
                          avgPeakListArgs = avgSetParams, origFGNames = names(fGroupsSet)))
}
