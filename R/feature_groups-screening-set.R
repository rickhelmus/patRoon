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
        
        for (col in unCols)
        {
            allCols <- getAllCols(col)
            if (length(allCols) > 0)
            {
                # take first non NA value
                scrInfo[, (col) := {
                    ret <- .SD[[1]] # set to first by default: in case all are NA and to ensure correct type
                    if (nrow(scrInfo) > 0)
                    {    
                        for (v in .SD)
                        {
                            if (!is.na(v))
                            {
                                ret <- v
                                break
                            }
                        }
                    }
                    ret
                }, by = seq_len(nrow(scrInfo)), .SDcols = allCols]
                scrInfo[, (allCols) := NULL]
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

    # retain sets form/comp ranks and estIDLevel    
    oldsi <- screenInfo(obj)
    for (col in c("formRank", "compRank", "estIDLevel"))
    {
        if (!is.null(oldsi[[col]]))
            newsi[oldsi, (col) := get(paste0("i.", col)), on = c("group", "name")]
    }

    obj@screenInfo <- newsi[]
    return(obj)
}

#' @param set \setsWF The name of the set.
#' @param sets \setsWF A \code{character} with name(s) of the sets to keep (or remove if \code{negate=TRUE}).
#'
#' @section Sets workflows: \setsWFClass{featureGroupsScreeningSet}{featureGroupsScreening}
#'
#'   \setsWFNewMethodsSO{featureGroupsScreeningUnset}{Only the screening results present in the specified set are kept.}
#'
#'   \setsWFChangedMethods{
#'
#'   \item \code{annotateSuspects} Suspect annotation is performed per set. Thus, formula/compound ranks, estimated
#'   identification levels etc are calculated for each set. Subsequently, these results are merged in the final
#'   \code{screenInfo}. In addition, an overall \code{formRank} and \code{compRank} column is created based on the
#'   rankings of the suspect candidate in the set consensus data. Furthermore, an overall \code{estIDLevel} is generated
#'   that is based on the 'best' estimated identification level among the sets data (\emph{i.e.} the lowest). In case
#'   there is a tie between sub-levels (\emph{e.g.} \samp{3a} and \samp{3b}), then the sub-level is stripped
#'   (\emph{e.g.} \samp{3}).
#'
#'   \item \code{filter} All filters releated to estimated identification levels and formula/compound rankings  are
#'   applied to the overall set data (see above). All others are applied to set specific data: in this case candidates
#'   are only removed if none of the set data confirms to the filter.
#'
#'   }
#'
#'   This class derives also from \code{\link{featureGroupsSet}}. Please see its documentation for more relevant details
#'   with sets workflows.
#'
#'   Note that the \code{formRank} and \code{compRank} columns are \emph{not} updated when the data is subset.
#'
#' @rdname featureGroupsScreening-class
#' @export
featureGroupsScreeningSet <- setClass("featureGroupsScreeningSet",
                                      slots = c(screenInfo = "data.table"),
                                      contains = c("featureGroupsSet", "workflowStepSet"))

setMethod("initialize", "featureGroupsScreeningSet",
          function(.Object, ...) callNextMethod(.Object, algorithm = "screening-set", ...))

#' @rdname featureGroupsScreening-class
#' @export
setMethod("screenInfo", "featureGroupsScreeningSet", function(obj) obj@screenInfo)

setMethod("mergedConsensusNames", "featureGroupsScreeningSet", function(obj) sets(obj))

#' @rdname featureGroupsScreening-class
#' @export
setMethod("show", "featureGroupsScreeningSet", function(object)
{
    callNextMethod(object)
    doScreeningShow(object)
})

#' @rdname featureGroupsScreening-class
#' @export
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

#' @rdname featureGroupsScreening-class
#' @export
setMethod("delete", "featureGroupsScreeningSet", function(obj, i = NULL, j = NULL, ...)
{
    oldn <- length(obj)
    obj <- callNextMethod()
    if (length(obj) != oldn)
        obj <- syncScreeningSetObjects(obj)
    return(obj)
})

#' @rdname featureGroupsScreening-class
#' @export
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

#' @rdname featureGroupsScreening-class
#' @export
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
    
    # clear old set cols if present
    cols <- c("formRank", "compRank", "estIDLevel")
    if (any(cols %in% names(screenInfo(fGroups))))
    {
        fGroups@screenInfo <- copy(fGroups@screenInfo)
        fGroups@screenInfo[, intersect(cols, names(screenInfo(fGroups))) := NULL]
    }
    
    fGroups <- syncScreeningSetObjects(fGroups)
    
    # add non set specific columns

    cols <- getAllSuspCols("formRank", names(screenInfo(fGroups)), mergedConsensusNames(fGroups))
    if (length(cols) > 0)
    {
        fGroups@screenInfo[, formRank := {
            if (is.na(formula) || !group %in% groupNames(formulas) || all(is.na(unlist(mget(cols)))))
                NA_integer_
            else
            {
                r <- which(formula == formulas[[group]]$neutral_formula)
                if (length(r) > 0)
                    r[1]
                else
                    NA_integer_
            }
        }, by = seq_len(nrow(fGroups@screenInfo))][]
    }

    cols <- getAllSuspCols("compRank", names(screenInfo(fGroups)), mergedConsensusNames(fGroups))
    if (length(cols) > 0)
    {
        fGroups@screenInfo[, compRank := {
            if (is.na(InChIKey) || !group %in% groupNames(compounds) || all(is.na(unlist(mget(cols)))))
                NA_integer_
            else
            {
                r <- which(getIKBlock1(InChIKey) == compounds[[group]]$InChIKey1)
                if (length(r) > 0)
                    r[1]
                else
                    NA_integer_
            }
        }, by = seq_len(nrow(fGroups@screenInfo))][]
    }
    
    # add best ID Levels
    cols <- getAllSuspCols("estIDLevel", names(screenInfo(fGroups)), mergedConsensusNames(fGroups))
    if (length(cols) > 0)
    {
        fGroups@screenInfo[, estIDLevel := {
            allIDs <- unlist(mget(cols))
            allIDs <- allIDs[!is.na(allIDs)]
            if (length(allIDs) == 0)
                NA_character_
            else
            {
                numIDs <- numericIDLevel(allIDs)
                whMin <- which(numIDs == min(numIDs))
                if (!allSame(allIDs[whMin]))
                    as.character(numIDs[whMin[1]]) # strip sublevel if not all the same
                else
                    allIDs[whMin[1]]
            }
        }, by = seq_len(nrow(fGroups@screenInfo))][]
    }

    return(fGroups)
})

#' @rdname featureGroupsScreening-class
#' @export
setMethod("filter", "featureGroupsScreeningSet", function(obj, ..., onlyHits = NULL,
                                                          selectHitsBy = NULL, selectBestFGroups = FALSE,
                                                          maxLevel = NULL, maxFormRank = NULL, maxCompRank = NULL,
                                                          minAnnSimForm = NULL, minAnnSimComp = NULL, minAnnSimBoth = NULL,
                                                          absMinFragMatches = NULL, relMinFragMatches = NULL,
                                                          negate = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertFlag, . ~ onlyHits + selectBestFGroups + negate, null.ok = c(TRUE, FALSE, FALSE), fixed = list(add = ac))
    checkmate::assertChoice(selectHitsBy, choices = c("intensity", "level"), null.ok = TRUE, add = ac)
    aapply(checkmate::assertCount, . ~ maxLevel + maxFormRank + maxCompRank + absMinFragMatches + relMinFragMatches,
           null.ok = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ minAnnSimForm + minAnnSimComp + minAnnSimBoth, null.ok = TRUE, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    oldsi <- screenInfo(obj)
    # NOTE: we do onlyHits later, as otherwise doSuspectFilter() might cause set synchronization (via delete()) whereas
    # the setObjects are not updated yet
    obj <- doSuspectFilter(obj, onlyHits = NULL, selectHitsBy, selectBestFGroups, maxLevel, maxFormRank, maxCompRank,
                           minAnnSimForm, minAnnSimComp, minAnnSimBoth, absMinFragMatches, relMinFragMatches, negate)
    newsi <- screenInfo(obj)
    suspFiltered <- !isTRUE(all.equal(oldsi, screenInfo(obj)))
    
    if (suspFiltered)
    {
        # update setObjects
        obj@setObjects <- lapply(obj@setObjects, function(so)
        {
            sosi <- copy(screenInfo(so))
            sosi[, keep := FALSE]
            sosi[newsi, keep := TRUE, on = c("group", "name")] # mark overlap
            so@screenInfo <- sosi[keep == TRUE, -"keep"]
            return(so)
        })
    }
    
    if (!is.null(onlyHits))
        obj <- doSuspectFilter(obj, onlyHits = onlyHits, selectHitsBy = NULL, selectBestFGroups = FALSE, maxLevel = NULL,
                               maxFormRank = NULL, maxCompRank = NULL, minAnnSimForm = NULL, minAnnSimComp = NULL,
                               minAnnSimBoth = NULL, absMinFragMatches = NULL, relMinFragMatches = NULL,
                               negate = negate)
    
    # filter functionality from fGroupsSet
    if (...length() > 0)
        obj <- callNextMethod(obj, ..., negate = negate)
    
    if (...length() > 0 || suspFiltered || !is.null(onlyHits))
        obj <- syncScreeningSetObjects(obj)

    return(obj)
})

#' @rdname suspect-screening
#' @export
setMethod("screenSuspects", "featureGroupsSet", function(fGroups, suspects, rtWindow, mzWindow,
                                                         adduct, skipInvalid, onlyHits)
{
    verifyNoAdductIonizationArg(adduct)
    
    if (checkmate::testDataFrame(suspects))
    {
        assertSuspectList(suspects, FALSE, skipInvalid)
        if (any(c("mz", "adduct", "fragments_mz") %in% names(suspects)) && length(sets(fGroups)) > 1)
            warning("The suspect list seems to contain an mz, adduct or fragments_mz column, ",
                    "which are generally specific to the ionization mode used. ",
                    "These columns most likely need to removed since the same suspect list will be used for all sets.",
                    call. = FALSE)
        suspects <- sapply(sets(fGroups), function(s) suspects, simplify = FALSE) # same for all sets
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
                                     groupQualities = copy(groupQualities(fGroups)),
                                     groupScores = copy(groupScores(fGroups)),
                                     annotations = copy(annotations(fGroups))))
})

#' @rdname suspect-screening
#' @export
setMethod("screenSuspects", "featureGroupsScreeningSet", function(fGroups, suspects, rtWindow, mzWindow,
                                                                  adduct, skipInvalid, onlyHits, amend = FALSE)
{
    aapply(checkmate::assertFlag, . ~ onlyHits + amend)
    
    fGroupsScreened <- callNextMethod(fGroups, suspects, rtWindow, mzWindow, adduct, skipInvalid, onlyHits)
    if (!amend)
        return(fGroupsScreened)
    
    # amend screening results
    
    fGroups@setObjects <- Map(fGroups@setObjects, fGroupsScreened@setObjects, f = function(so, sos)
    {
        so@screenInfo <- rbind(so@screenInfo, sos@screenInfo, fill = TRUE)
        so@screenInfo <- unique(so@screenInfo, by = c("name", "group"))
        return(so)
    })
    fGroups <- syncScreeningSetObjects(fGroups)
    
    if (onlyHits)
        fGroups <- fGroups[, fGroups@screenInfo$group]
    
    return(fGroups)
})


#' @rdname featureGroupsScreening-class
#' @export
featureGroupsSetScreeningUnset <- setClass("featureGroupsSetScreeningUnset",
                                           contains = "featureGroupsScreening")

#' @rdname featureGroupsScreening-class
#' @export
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
