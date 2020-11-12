#' @include main.R
#' @include feature_groups-set.R
#' @include feature_groups-screening.R
NULL

# merges screening info from screenInfo slots or as.data.table() tables
mergeScreeningSetInfos <- function(setObjects, sInfos = lapply(setObjects, screenInfo),
                                   rmSetCols = TRUE)
{
    rmCols <- c("mz", "fragments_mz")

    if (length(setObjects) > 1)
    {
        sets <- names(setObjects)
        
        scrInfo <- ReduceWithArgs(x = sInfos, paste0("-", names(setObjects)),
                                  f = function(l, r, sl, sr) merge(l, r, suffixes = c(sl, sr),
                                                                   by = c("name", "group"), all = TRUE))
        
        unCols <- c("rt", "formula", "SMILES", "InChI", "InChIKey", "neutralMass", "d_rt", "d_mz", "fragments_formula")
        
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
        
        if (rmSetCols)
            scrInfo[, (getAllCols(rmCols)) := NULL]
    }
    else if (length(sInfos) == 1)
    {
        scrInfo <- copy(sInfos[[1]])
        if (rmSetCols)
        {
            rmCols <- intersect(rmCols, names(scrInfo))
            if (length(rmCols) > 0)
                scrInfo[, (intersect(rmCols, names(scrInfo))) := NULL]
        }
    }
    else
        scrInfo <- data.table()
    
    return(scrInfo[])
}

syncScreeningSetObjects <- function(obj)
{
    # BUG? can't call "[" directly here to subset??
    # obj@setObjects <- lapply(obj@setObjects, "[", i = analyses(obj), j = groupNames(obj))
    obj@setObjects <- lapply(obj@setObjects, function(x) x[analyses(obj), groupNames(obj)])
    obj@setObjects <- pruneList(obj@setObjects, checkEmptyElements = TRUE)
    obj@screenInfo <- mergeScreeningSetInfos(obj@setObjects)
    return(obj)
}

featureGroupsScreeningSet <- setClass("featureGroupsScreeningSet",
                                      slots = c(screenInfo = "data.table"),
                                      contains = c("featureGroupsSet", "workflowStepSet"))

setMethod("initialize", "featureGroupsScreeningSet",
          function(.Object, ...) callNextMethod(.Object, algorithm = "screening-set", ...))

setMethod("screenInfo", "featureGroupsScreeningSet", function(obj) obj@screenInfo)

setMethod("[", c("featureGroupsScreeningSet", "ANY", "ANY", "missing"), function(x, i, j, ..., rGroups, sets = NULL, drop = TRUE)
{
    x <- callNextMethod(x, i, j, ..., rGroups = rGroups, sets = sets, drop = drop)
    return(syncScreeningSetObjects(x))
})

setMethod("as.data.table", "featureGroupsScreeningSet", function(x, ..., collapseSuspects = ",",
                                                                 onlyHits = FALSE)
{
    # UNDONE: document that collapseSuspect!=NULL && features==TRUE will give lots of rows (per feature and per suspect)
    
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

setMethod("annotateSuspects", "featureGroupsScreeningSet", function(fGroups, MSPeakLists, formulas,
                                                                    compounds, ...)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertClass, . ~ MSPeakLists + formulas + compounds,
           c("MSPeakListsSet", "formulasSet", "compoundsSet"), null.ok = TRUE, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    ionOrNULL <- function(o) if (is.null(o)) rep(list(NULL), length(sets(fGroups))) else lapply(sets(fGroups), ionize, obj = o)
    
    ionizedMSPeakLists <- ionOrNULL(MSPeakLists)
    ionizedFormulas <- ionOrNULL(formulas)
    ionizedCompounds <- ionOrNULL(compounds)
    
    fGroups@setObjects <- mapply(setObjects(fGroups), ionizedMSPeakLists, ionizedFormulas, ionizedCompounds,
                                 FUN = annotateSuspects, MoreArgs = list(...), SIMPLIFY = FALSE)
    
    return(syncScreeningSetObjects(fGroups))
})

setMethod("filter", "featureGroupsScreeningSet", function(obj, ..., onlyHits = FALSE,
                                                          selectHitsBy = NULL, selectFGroupsBy = NULL,
                                                          maxLevel = NULL, maxFormRank = NULL, maxCompRank = NULL,
                                                          minAnnMSMSSim = NULL, minFragMatches = NULL, negate = FALSE)
{
    # filter functionality from fGroupsSet
    obj <- callNextMethod(obj, ..., negate = negate)
    obj <- syncScreeningSetObjects(obj)
    
    # filter functionality from screening (no need to pass ...)
    obj@setObjects <- lapply(obj@setObjects, filter, onlyHits = onlyHits, selectHitsBy = selectHitsBy,
                             selectFGroupsBy = selectFGroupsBy, maxLevel = maxLevel, maxFormRank = maxFormRank,
                             maxCompRank = maxCompRank, minAnnMSMSSim = minAnnMSMSSim, minFragMatches = minFragMatches,
                             negate = negate)
    # --> groups may have been removed
    obj <- obj[, unique(unlist(sapply(obj@setObjects, groupNames)))]
    obj <- syncScreeningSetObjects(obj)
    
    return(obj)
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


featureGroupsSetScreeningIonized <- setClass("featureGroupsSetScreeningIonized",
                                             contains = "featureGroupsScreening")
setMethod("ionize", "featureGroupsScreeningSet", function(obj, sets)
{
    iobj <- callNextMethod()
    
    if (!is.null(sets))
        obj <- obj[, sets = sets]
    sInfo <- mergeScreeningSetInfos(setObjects(obj), rmSetCols = FALSE)
    
    ret <- featureGroupsSetScreeningIonized(screenInfo = sInfo, groups = groups(iobj),
                                            groupInfo = groupInfo(iobj), analysisInfo = analysisInfo(iobj),
                                            features = getFeatures(iobj), ftindex = groupFeatIndex(iobj))
    # override after constructing: parent constructor already sets algorithm,
    # which results in error about double assignment
    ret@algorithm <- paste0(algorithm(obj), "_ionized")
    return(ret)
})
