#' @include main.R
#' @include feature_groups-set.R
#' @include feature_groups-screening.R
NULL

# merges screening info from screenInfo slots or as.data.table() tables
mergeScreeningSetInfos <- function(setObjects, sInfos = lapply(setObjects, screenInfo),
                                   rmSetCols = TRUE, markSets = TRUE)
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
                # take first non NA value
                scrInfo[, (col) := {
                    ret <- .SD[[1]] # set to first by default: in case all are NA and to ensure correct type
                    for (v in .SD)
                    {
                        if (!is.na(v))
                        {
                            ret <- v
                            break
                        }
                    }
                    ret
                }, by = seq_len(nrow(scrInfo)), .SDcols = allCols]
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
    
    if (nrow(scrInfo) > 0 && markSets)
    {
        # add set presence
        scrInfo[, sets := mapply(name, group, FUN = function(n, g) {
            ret <- names(setObjects)
            return(paste0(ret[sapply(setObjects, function(so) screenInfo(so)[name == n & group == g, .N] > 0)],
                          collapse = ","))
        })]
    }
    
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

setMethod("[", c("featureGroupsScreeningSet", "ANY", "ANY", "missing"), function(x, i, j, ..., rGroups,
                                                                                 suspects = NULL, sets = NULL,
                                                                                 drop = TRUE)
{
    checkmate::assertCharacter(suspects, null.ok = TRUE)
    
    x <- callNextMethod(x, i, j, ..., rGroups = rGroups, sets = sets, drop = drop)
    x <- syncScreeningSetObjects(x)
    
    if (!is.null(suspects))
    {
        x@setObjects <- lapply(x@setObjects, "[", suspects = suspects)
        # --> groups may have been removed
        x <- x[, unique(unlist(sapply(x@setObjects, groupNames)))]
        x <- syncScreeningSetObjects(x)
    }    
    
    return(x)
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
                                                          onlyHits = onlyHits),
                                     markSets = is.null(collapseSuspects))
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
    
    unsetOrNULL <- function(o) if (is.null(o)) rep(list(NULL), length(sets(fGroups))) else lapply(sets(fGroups), unset, obj = o)
    
    unsetMSPeakLists <- unsetOrNULL(MSPeakLists)
    unsetFormulas <- unsetOrNULL(formulas)
    unsetCompounds <- unsetOrNULL(compounds)
    
    fGroups@setObjects <- mapply(setObjects(fGroups), unsetMSPeakLists, unsetFormulas, unsetCompounds,
                                 FUN = annotateSuspects, MoreArgs = list(...), SIMPLIFY = FALSE)
    
    return(syncScreeningSetObjects(fGroups))
})

setMethod("filter", "featureGroupsScreeningSet", function(obj, ..., onlyHits = NULL,
                                                          selectHitsBy = NULL, selectBestFGroups = FALSE,
                                                          maxLevel = NULL, maxFormRank = NULL, maxCompRank = NULL,
                                                          minAnnSimForm = NULL, minAnnSimComp = NULL, minAnnSimBoth = NULL,
                                                          absMinFragMatches = NULL, relMinFragMatches = NULL,
                                                          sets = NULL, negate = FALSE)
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
    
    # filter functionality from fGroupsSet
    obj <- callNextMethod(obj, ..., negate = negate)
    obj <- syncScreeningSetObjects(obj)
    
    # filter functionality from screening (no need to pass ...)
    obj@setObjects <- lapply(obj@setObjects, filter, onlyHits = onlyHits, selectHitsBy = selectHitsBy,
                             selectBestFGroups = selectBestFGroups, maxLevel = maxLevel, maxFormRank = maxFormRank,
                             maxCompRank = maxCompRank, minAnnSimForm = minAnnSimForm, minAnnSimComp = minAnnSimComp,
                             minAnnSimBoth = minAnnSimBoth, absMinFragMatches = absMinFragMatches,
                             relMinFragMatches = relMinFragMatches, negate = negate)
    # --> groups may have been removed
    obj <- obj[, unique(unlist(sapply(obj@setObjects, groupNames)))]
    obj <- syncScreeningSetObjects(obj)
    
    return(obj)
})

setMethod("screenSuspects", "featureGroupsSet", function(fGroups, suspects, rtWindow, mzWindow,
                                                         adduct, skipInvalid, onlyHits)
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
    
    unsetFGroupsList <- sapply(sets(fGroups), unset, obj = fGroups, simplify = FALSE)
    setObjects <- mapply(unsetFGroupsList, suspects, adducts(fGroups), SIMPLIFY = FALSE,
                         FUN = function(fg, s, a) screenSuspects(fg, s, rtWindow, mzWindow, a,
                                                                 skipInvalid, onlyHits))
    
    return(featureGroupsScreeningSet(screenInfo = mergeScreeningSetInfos(setObjects), setObjects = setObjects,
                                     groups = copy(groupTable(fGroups)), analysisInfo = analysisInfo(fGroups),
                                     groupInfo = groupInfo(fGroups), features = getFeatures(fGroups),
                                     ftindex = copy(groupFeatIndex(fGroups))))
})


featureGroupsSetScreeningUnset <- setClass("featureGroupsSetScreeningUnset",
                                           contains = "featureGroupsScreening")
setMethod("unset", "featureGroupsScreeningSet", function(obj, set)
{
    assertSets(obj, set, FALSE)
    
    uobj <- callNextMethod()
    
    obj <- obj[, sets = set]
    sInfo <- mergeScreeningSetInfos(setObjects(obj), rmSetCols = FALSE, markSets = FALSE)
    
    ret <- featureGroupsSetScreeningUnset(screenInfo = sInfo, groups = groupTable(uobj),
                                          groupInfo = groupInfo(uobj), analysisInfo = analysisInfo(uobj),
                                          features = getFeatures(uobj), ftindex = groupFeatIndex(uobj))
    # override after constructing: parent constructor already sets algorithm,
    # which results in error about double assignment
    ret@algorithm <- paste0(algorithm(obj), "_unset")
    return(ret)
})
