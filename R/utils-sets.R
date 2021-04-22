makeSetAlgorithm <- function(setObjects)
{
    if (length(setObjects) == 0)
        setAlgo <- "empty"
    else
        setAlgo <- paste0(unique(sapply(setObjects, algorithm)), collapse = ",")
    return(paste(setAlgo, "set", sep = "-"))
}

# checks if NULL, verifies if all sets are present and syncs sets order while unsetting
checkAndUnSetOther <- function(targetSets, obj, objName, handleNULL = FALSE)
{
    if (handleNULL && is.null(obj))
        return(setNames(as.list(rep(list(NULL), length(targetSets))), targetSets))
    
    missingSets <- setdiff(targetSets, sets(obj))
    if (length(missingSets) > 0)
        stop("Missing sets in ", objName, ": ", paste0(missingSets, collapse = ", "))
    return(sapply(targetSets, unset, obj = obj, simplify = FALSE))
}

assertAndGetMSPLSetsArgs <- function(fGroupsSet, MSPeakListsSet)
{
    checkmate::assertClass(MSPeakListsSet, "MSPeakListsSet")
    unsetMSPeaksList <- checkAndUnSetOther(sets(fGroupsSet), MSPeakListsSet, "MSPeakLists")
    return(lapply(unsetMSPeaksList, function(x) list(MSPeakLists = x)))
}

prepareMakeSetAdducts <- function(objects, adducts, labels)
{
    adducts <- lapply(adducts, checkAndToAdduct, .var.name = "adducts")
    adducts <- rep(adducts, length.out = length(objects))
    
    if (!is.null(labels))
        names(adducts) <- labels
    else
        names(adducts) <- make.unique(ifelse(sapply(adducts, "slot", "charge") < 0, "negative", "positive"))
    
    return(adducts)
}
