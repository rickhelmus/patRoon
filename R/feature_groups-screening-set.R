#' @include main.R
#' @include feature_groups-set.R
#' @include feature_groups-screening.R
NULL

# merges screening info from screenInfo slots
mergeScreeningSetInfos <- function(setObjects, sInfos = lapply(setObjects, screenInfo), rmSetCols = TRUE)
{
    rmCols <- c("mz", "fragments_mz")
    unCols <- c("rt", "formula", "SMILES", "InChI", "InChIKey", "neutralMass",  "d_rt", "d_mz", "fragments_formula")
    
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
    
    if (length(setObjects) > 1)
    {
        sets <- names(setObjects)
        
        getAllCols <- function(cols)
        {
            cols <- unlist(lapply(cols, paste0, "-", sets))
            return(cols[sapply(cols, function(x) !is.null(scrInfo[[x]]))])
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
        scrInfo <- renameDupCols(scrInfo, paste0("-", names(setObjects)[1]))
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
    }
    else
        scrInfo[, sets := character()]
    
    return(scrInfo[])
}

syncScreeningSetObjects <- function(obj)
{
    # BUG? can't call "[" directly here to subset??
    # obj@setObjects <- lapply(obj@setObjects, "[", i = analyses(obj), j = groupNames(obj))
    obj@setObjects <- lapply(obj@setObjects, function(x) x[analyses(obj), groupNames(obj)])
    newsi <- mergeScreeningSetInfos(obj@setObjects)

    # retain form/comp ranks    
    oldsi <- screenInfo(obj)
    for (col in c("formRank", "compRank"))
    {
        if (!is.null(oldsi[[col]]))
            newsi[, (col) := oldsi[group %in% newsi$group][[col]]]
    }

    obj@screenInfo <- newsi[]
    return(obj)
}

featureGroupsScreeningSet <- setClass("featureGroupsScreeningSet",
                                      slots = c(screenInfo = "data.table"),
                                      contains = c("featureGroupsSet", "workflowStepSet"))

setMethod("initialize", "featureGroupsScreeningSet",
          function(.Object, ...) callNextMethod(.Object, algorithm = "screening-set", ...))

setMethod("screenInfo", "featureGroupsScreeningSet", function(obj) obj@screenInfo)

#' @export
setMethod("show", "featureGroupsScreeningSet", function(object)
{
    callNextMethod(object)
    doScreeningShow(object)
})

setMethod("[", c("featureGroupsScreeningSet", "ANY", "ANY", "missing"), function(x, i, j, ..., rGroups,
                                                                                 suspects = NULL, sets = NULL,
                                                                                 drop = TRUE)
{
    checkmate::assertCharacter(suspects, null.ok = TRUE)
    assertSets(x, sets, TRUE)
    
    x <- callNextMethod(x, i, j, ..., rGroups = rGroups, sets = sets, drop = drop)
    
    if (!is.null(suspects))
    {
        x@setObjects <- lapply(x@setObjects, "[", suspects = suspects)
        # --> groups may have been removed
        x <- x[, unique(unlist(lapply(x@setObjects, groupNames)))]
    }    

    if (!is.null(sets))
    {
        x@setObjects <- x@setObjects[sets]
        x <- syncScreeningSetObjects(x)
    }
    
    return(x)
})

#' @export
setMethod("delete", "featureGroupsScreeningSet", function(obj, i = NULL, j = NULL, ...)
{
    obj <- callNextMethod()
    obj <- syncScreeningSetObjects(obj)
    return(obj)
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

setMethod("annotateSuspects", "featureGroupsScreeningSet", function(fGroups, MSPeakLists, formulas, compounds, ...)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertClass, . ~ MSPeakLists + formulas + compounds,
           c("MSPeakListsSet", "formulasSet", "compoundsSet"), null.ok = TRUE, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    unsetMSPeakLists <- checkAndUnSetOther(sets(fGroups), MSPeakLists, "MSPeakLists", TRUE)
    unsetFormulas <- checkAndUnSetOther(sets(fGroups), formulas, "formulas", TRUE)
    unsetCompounds <- checkAndUnSetOther(sets(fGroups), compounds, "compounds", TRUE)
    
    fGroups@setObjects <- Map(setObjects(fGroups), unsetMSPeakLists, unsetFormulas, unsetCompounds,
                              f = annotateSuspects, MoreArgs = list(...))
    
    # clear old rank cols if present
    rankCols <- c("formRank", "compRank")
    if (any(rankCols %in% names(screenInfo(fGroups))))
        fGroups@screenInfo[, intersect(rankCols, names(screenInfo(fGroups))) := NULL]
    
    fGroups <- syncScreeningSetObjects(fGroups)
    
    # add non set specific ranks
    allRankCols <- getAllSuspSetCols(c("formRank", "compRank"), names(screenInfo(fGroups)), sets(fGroups))
    if (any(grepl("^formRank", allRankCols)))
    {
        fGroups@screenInfo[!is.na(formula) & group %in% groupNames(formulas), formRank := mapply(group, formula, FUN = function(g, f)
        {
            unFTable <- unique(formulas[[g]], by = "neutral_formula")
            r <- which(f == unFTable$neutral_formula)
            return(if (length(r) > 0) r[1] else NA_integer_)
        })][]
    }

    if (any(grepl("^compRank", allRankCols)))
    {
        fGroups@screenInfo[!is.na(InChIKey) & group %in% groupNames(compounds), compRank := mapply(group, InChIKey, FUN = function(g, ik)
        {
            r <- which(getIKBlock1(ik) == compounds[[g]]$InChIKey1)
            return(if (length(r) > 0) r[1] else NA_integer_)
        })][]
    }

    return(fGroups)
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
    if (!is.null(adduct))
        stop("adduct argument not supported for sets!")
    
    if (checkmate::testDataFrame(suspects))
    {
        assertSuspectList(suspects, FALSE, skipInvalid)
        suspects <- sapply(sets(fGroups), function(s) suspects, simplify = FALSE) # same for all set
    }
    else
    {
        checkmate::assertList(suspects, "data.frame", any.missing = FALSE, all.missing = FALSE,
                              len = length(sets(fGroups)))
        checkmate::assert(
            checkmate::checkNames(names(suspects), "unnamed"),
            checkmate::checkNames(names(suspects), "unique", must.include = sets(fGroups)),
            .var.name = "suspects"
        )
        if (checkmate::testNames(names(suspects), "unnamed"))
            names(suspects) <- sets(fGroups)
    }
    
    # sync order
    suspects <- suspects[sets(fGroups)]
    
    unsetFGroupsList <- sapply(sets(fGroups), unset, obj = fGroups, simplify = FALSE)
    setObjects <- Map(unsetFGroupsList, suspects,
                      f = function(fg, s) screenSuspects(fg, s, rtWindow = rtWindow, mzWindow = mzWindow,
                                                         adduct = NULL, skipInvalid = skipInvalid, onlyHits = onlyHits))
    
    scr <- mergeScreeningSetInfos(setObjects)
    if (onlyHits)
        fGroups <- fGroups[, scr$group]
    
    return(featureGroupsScreeningSet(screenInfo = scr, setObjects = setObjects,
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
    sInfo <- mergeScreeningSetInfos(setObjects(obj), rmSetCols = FALSE)
    if (length(sInfo) > 0)
    {
        sInfo[, sets := NULL]
        # restore set specific columns
        setnames(sInfo, sub(paste0("\\-", set, "$"), "", names(sInfo)))
    }
    
    ret <- featureGroupsSetScreeningUnset(screenInfo = sInfo, groups = groupTable(uobj),
                                          groupInfo = groupInfo(uobj), analysisInfo = analysisInfo(uobj),
                                          features = getFeatures(uobj), ftindex = groupFeatIndex(uobj),
                                          annotations = annotations(uobj))
    # override after constructing: parent constructor already sets algorithm,
    # which results in error about double assignment
    ret@algorithm <- paste0(algorithm(obj), "_unset")
    return(ret)
})
