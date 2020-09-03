#' @include main.R
#' @include feature_groups-set.R
#' @include feature_groups-screening.R
NULL

# merges screening info from screenInfo slots or as.data.table() tables
mergeScreeningSetInfos <- function(setObjects, sInfos = lapply(setObjects, screenInfo))
{
    sets <- names(setObjects)
    
    scrInfo <- ReduceWithArgs(x = sInfos, paste0("-", names(setObjects)),
                              f = function(l, r, sl, sr) merge(l, r, suffixes = c(sl, sr),
                                                               by = c("name", "group"), all = TRUE))
    
    unCols <- c("rt", "formula", "d_rt", "d_mz", "fragments_formula")
    rmCols <- c("mz", "fragments_mz")
    
    getAllCols <- function(cols)
    {
        cols <- unlist(lapply(cols, paste0, "-", sets(fGroups)))
        return(cols[sapply(cols, function(x) !is.null(scrInfo[[x]]))])
    }
    
    for (col in unCols)
    {
        allCols <- getAllCols(col)
        if (length(allCols) > 0)
        {
            scrInfo[, (col) := get(allCols[1])] # just take first
            scrInfo[, (allCols) := NULL]
        }
    }
    
    scrInfo[, (getAllCols(rmCols)) := NULL]
    
    return(scrInfo[])
}

syncScreeningSetObjects <- function(obj)
{
    # BUG? can't call "[" directly here to subset??
    # obj@setObjects <- lapply(obj@setObjects, "[", i = analyses(obj), j = groupNames(obj))
    obj@setObjects <- lapply(obj@setObjects, function(x) x[analyses(obj), groupNames(obj)])
    obj@screenInfo <- mergeScreeningSetInfos(obj@setObjects)
    return(obj)
}

featureGroupsScreeningSet <- setClass("featureGroupsScreeningSet",
                                      slots = c(screenInfo = "data.table"),
                                      contains = c("featureGroupsSet", "workflowStepSet"))

setMethod("initialize", "featureGroupsScreeningSet",
          function(.Object, ...) callNextMethod(.Object, algorithm = "screening-set", ...))

setMethod("screenInfo", "featureGroupsScreeningSet", function(obj) obj@screenInfo)

setMethod("[", c("featureGroupsScreeningSet", "ANY", "ANY", "missing"), function(x, i, j, ..., sets = NULL, drop = TRUE)
{
    x <- callNextMethod()
    return(syncScreeningSetObjects(x))
})

setMethod("as.data.table", "featureGroupsScreeningSet", function(x, ..., collapseSuspects = TRUE,
                                                                 onlyHits = FALSE)
{
    # UNDONE: document that collapseSuspect==TRUE && features==TRUE will give lots of rows (per feature and per suspect)
    
    checkmate::assertFlag(onlyHits, add = ac)
    
    ret <- callNextMethod(x, ...)
    
    # get tables from setObjects. Note that we only want to get the screenInfo per group, so other arguments from
    # ... can be ignored.
    dtSets <- mergeScreeningSetInfos(x@setObjects, lapply(x@setObjects, as.data.table,
                                                          collapseSuspects = collapseSuspects,
                                                          onlyHits = onlyHits))
    dtSets <- dtSets[, c("group", setdiff(names(dtSets), names(ret))), with = FALSE] # only keep unique columns (and group)
    
    return(merge(ret, dtSets, by = "group", all.x = !onlyHits))
})              

setMethod("groupFeaturesScreening", "featureGroupsSet", function(fGroups, suspects, rtWindow, mzWindow,
                                                                 adduct, skipInvalid)
{
    # UNDONE: remove argument (and from generic?)
    if (!is.null(adduct))
        stop("adduct argument not supported for sets!")
    
    if (checkmate::testDataFrame(suspects))
    {
        assertSuspectList(suspects, adducts(fGroups)[1], skipInvalid)
        suspects <- sapply(sets(fGroups), function(s) suspects, simplify = FALSE) # same for all set
    }
    else
    {
        checkmate::assertList(suspects, "data.frame", any.missing = FALSE, all.missing = FALSE,
                              len = length(sets(fGroups)), names = "unique")
        checkmate::assertSubset(names(suspects), sets(fGroups), empty.ok = FALSE)
    }
    
    # sync order
    suspects <- suspects[sets(fGroups)]
    
    ionizedFGroupsList <- sapply(sets(fGroups), ionize, obj = fGroups, simplify = FALSE)
    ionizedFGScr <- mapply(ionizedFGroupsList, suspects, adducts(fGroups), SIMPLIFY = FALSE,
                           FUN = function(fg, s, a) groupFeaturesScreening(fg, s, rtWindow, mzWindow, a,
                                                                           skipInvalid, onlyHits))
    
    return(featureGroupsScreeningSet(screenInfo = mergeScreeningSetInfos(ionizedFGScr), setObjects = ionizedFGScr,
                                     groups = copy(groups(fGroups)), analysisInfo = analysisInfo(fGroups),
                                     groupInfo = groupInfo(fGroups), features = getFeatures(fGroups),
                                     ftindex = copy(groupFeatIndex(fGroups))))
})

