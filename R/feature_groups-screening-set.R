#' @include main.R
#' @include feature_groups-set.R
#' @include feature_groups-screening.R
NULL

# merges screening info from screenInfo slots or as.data.table() tables
mergeScreeningSetInfos <- function(setObjects, setThreshold, sInfos = lapply(setObjects, screenInfo), rmSetCols = TRUE)
{
    rmCols <- c("mz", "fragments_mz")
    unCols <- c("rt", "formula", "SMILES", "InChI", "InChIKey", "neutralMass",  "d_rt", "d_mz", "fragments_formula")
    
    if (length(setObjects) > 1)
    {
        sets <- names(setObjects)
        
        getAllCols <- function(cols)
        {
            cols <- unlist(lapply(cols, paste0, "-", sets))
            return(cols[sapply(cols, function(x) !is.null(scrInfo[[x]]))])
        }
        
        renameDupCols <- function(si, suf)
        {
            cols <- setdiff(names(si), c("name", "group", unCols))
            if (length(cols) > 0)
            {
                si <- copy(si)
                setnames(si, cols, paste0(cols, suf))
            }
            return(si)
        }
        
        scrInfo <- ReduceWithArgs(x = sInfos, paste0("-", names(setObjects)), f = function(l, r, sl, sr)
        {
            # suffix non-unique columns columns
            l <- copy(l); r <- copy(r)
            
            merge(renameDupCols(l, sl), renameDupCols(r, sr), suffixes = c(sl, sr), by = c("name", "group"),
                  all = TRUE)
        })
        
        if (nrow(scrInfo) > 0)
        {
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
        }
        
        if (rmSetCols)
        {
            rmc <- getAllCols(rmCols)
            if (length(rmc) > 0)
                scrInfo[, (rmc) := NULL]
        }
    }
    else if (length(sInfos) == 1)
    {
        scrInfo <- copy(sInfos[[1]])
        if (rmSetCols)
        {
            rmc <- intersect(rmCols, names(scrInfo))
            if (length(rmc) > 0)
                scrInfo[, (rmc) := NULL]
        }
    }
    else
        scrInfo <- data.table()
    
    if (nrow(scrInfo) > 0)
    {
        # add set presence
        scrInfo[, sets := mapply(name, group, FUN = function(n, g) {
            ret <- names(setObjects)
            return(paste0(ret[sapply(setObjects, function(so) screenInfo(so)[name == n & group == g, .N] > 0)],
                          collapse = ","))
        })]
        scrInfo[, setCoverage := (sapply(sets, countCharInStr, ch = ",") + 1) / length(setObjects)]
        if (setThreshold > 0)
            scrInfo <- scrInfo[setCoverage >= setThreshold]
    }
    
    return(scrInfo[])
}

syncScreeningSetObjects <- function(obj)
{
    # BUG? can't call "[" directly here to subset??
    # obj@setObjects <- lapply(obj@setObjects, "[", i = analyses(obj), j = groupNames(obj))
    obj@setObjects <- lapply(obj@setObjects, function(x) x[analyses(obj), groupNames(obj)])
    obj@setObjects <- pruneList(obj@setObjects, checkEmptyElements = TRUE)
    obj@screenInfo <- mergeScreeningSetInfos(obj@setObjects, obj@setThreshold)
    return(obj)
}

featureGroupsScreeningSet <- setClass("featureGroupsScreeningSet",
                                      slots = c(screenInfo = "data.table", setThreshold = "numeric"),
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
    if (nrow(ret) > 0)
        ret <- mergeScreenInfoWithDT(ret, screenInfo(x), collapseSuspects, onlyHits)
    return(ret)    
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
                                                         adduct, skipInvalid, onlyHits,
                                                         setThreshold = 0)
{
    if (!is.null(adduct))
        stop("adduct argument not supported for sets!")
    
    if (checkmate::testDataFrame(suspects))
    {
        assertSuspectList(suspects, TRUE, skipInvalid)
        suspects <- sapply(sets(fGroups), function(s) suspects, simplify = FALSE) # same for all set
    }
    else
    {
        checkmate::assertList(suspects, "data.frame", any.missing = FALSE, all.missing = FALSE,
                              len = length(sets(fGroups)), names = "unique")
        checkmate::assertSubset(names(suspects), sets(fGroups), empty.ok = FALSE)
    }
    checkmate::assertNumber(setThreshold, lower = 0, upper = 1)
    
    # sync order
    suspects <- suspects[sets(fGroups)]
    
    unsetFGroupsList <- sapply(sets(fGroups), unset, obj = fGroups, simplify = FALSE)
    setObjects <- Map(unsetFGroupsList, suspects,
                      f = function(fg, s) screenSuspects(fg, s, rtWindow = rtWindow, mzWindow = mzWindow,
                                                         adduct = NULL, skipInvalid = skipInvalid, onlyHits = onlyHits))
    
    scr <- mergeScreeningSetInfos(setObjects, setThreshold)
    if (onlyHits)
        fGroups <- fGroups[, scr$group]
    
    return(featureGroupsScreeningSet(screenInfo = scr, setThreshold = setThreshold, setObjects = setObjects,
                                     groupAlgo = fGroups@groupAlgo, groupArgs = fGroups@groupArgs,
                                     groupVerbose = fGroups@groupVerbose, groups = copy(groupTable(fGroups)),
                                     analysisInfo = analysisInfo(fGroups), groupInfo = groupInfo(fGroups),
                                     features = getFeatures(fGroups), ftindex = copy(groupFeatIndex(fGroups)),
                                     annotations = copy(annotations(fGroups))))
})


featureGroupsSetScreeningUnset <- setClass("featureGroupsSetScreeningUnset",
                                           contains = "featureGroupsScreening")
setMethod("unset", "featureGroupsScreeningSet", function(obj, set)
{
    assertSets(obj, set, FALSE)
    
    uobj <- callNextMethod()
    
    obj <- obj[, sets = set]
    sInfo <- mergeScreeningSetInfos(setObjects(obj), obj@setThreshold, rmSetCols = FALSE)
    sInfo[, c("sets", "setCoverage") := NULL]
    
    ret <- featureGroupsSetScreeningUnset(screenInfo = sInfo, groups = groupTable(uobj),
                                          groupInfo = groupInfo(uobj), analysisInfo = analysisInfo(uobj),
                                          features = getFeatures(uobj), ftindex = groupFeatIndex(uobj),
                                          annotations = annotations(uobj))
    # override after constructing: parent constructor already sets algorithm,
    # which results in error about double assignment
    ret@algorithm <- paste0(algorithm(obj), "_unset")
    return(ret)
})
