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
                           slots = c(setObjects = "list"),
                           contains = "MSPeakLists")

setMethod("initialize", "MSPeakListsSet",
          function(.Object, ...) callNextMethod(.Object, algorithm = "set", ...))

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
    return(MSPeakListsSet(setObjects = ionizedMSPeakLists, peakLists = combPL, metadata = list(),
                          avgPeakListArgs = avgSetParams, origFGNames = names(fGroupsSet)))
}
