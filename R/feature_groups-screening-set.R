# SPDX-FileCopyrightText: 2016 - 2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
#' @include feature_groups-set.R
#' @include feature_groups-screening.R
NULL

# merges screening info from screenInfo slots
mergeScreeningSetInfos <- function(setObjects, sInfos = lapply(setObjects, screenInfo), rmSetCols = TRUE)
{
    rmCols <- c("mz", "fragments_mz")
    unCols <- c("rt", "formula", "SMILES", "InChI", "InChIKey", "neutralMass", "d_rt", "d_mz", "LC50_SMILES")
    
    renameDupCols <- function(si, suf, all)
    {
        except <- c("name", "group")
        if (!all)
            except <- c(except, unCols)
        cols <- setdiff(names(si), except)
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
        
        scrInfo <- ReduceWithArgs(x = sInfos, names(setObjects), f = function(l, r, sl, sr)
        {
            # suffix non-unique columns columns
            l <- copy(l); r <- copy(r)
            ssl <- paste0("-", sl); ssr <- paste0("-", sr)
            
            if (sr == names(setObjects)[2])
                l <- renameDupCols(l, ssl, TRUE) # rename left only once (ie when right is the second set)
            
            merge(l, renameDupCols(r, ssr, TRUE), suffixes = c(ssl, ssr), by = c("name", "group"), all = TRUE)
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
        scrInfo <- renameDupCols(scrInfo, paste0("-", names(setObjects)[1]), FALSE)
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
            return(paste0(ret[sapply(sInfos, function(si) si[name == n & group == g, .N] > 0)],
                          collapse = ","))
        })]
    }
    else
        scrInfo[, sets := character()]
    
    return(scrInfo[])
}

syncScreeningSetObjects <- function(obj, unsetFGroups)
{
    newsi <- mergeScreeningSetInfos(unsetFGroups)

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
#'   \item \code{filter} All filters related to estimated identification levels and formula/compound rankings  are
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
                                      contains = "featureGroupsSet")

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
setMethod("[", c("featureGroupsScreeningSet", "ANY", "ANY", "missing"), function(x, i, j, ..., ni, rGroups,
                                                                                 suspects = NULL, sets = NULL,
                                                                                 reorder = FALSE, drop = TRUE)
{
    checkmate::assertCharacter(suspects, null.ok = TRUE)
    assertSets(x, sets, TRUE)
    
    curSets <- get("sets", pos = 2)(x)
    
    x <- callNextMethod(x, i, j, ..., ni = ni, rGroups = rGroups, sets = sets, reorder = reorder, drop = drop)
    
    if (!is.null(suspects))
        x <- x[, x@screenInfo[name %in% suspects]$group]

    if (!is.null(sets))
    {
        # get rid of set specific columns for removed sets
        rmSets <- setdiff(curSets, sets)
        rmCols <- getAllMergedConsCols(names(x@screenInfo), rmSets)
        if (length(rmCols) > 0)
            x@screenInfo[, (rmCols) := NULL]
        
        # update sets assignments and get rid of set specific rows
        newSets <- sets
        x@screenInfo[, sets := {
            sv <- unlist(strsplit(sets, ",", fixed = TRUE))
            paste0(intersect(newSets, sv), collapse = ",")
        }, by = seq_len(nrow(x@screenInfo))][]
        x@screenInfo <- x@screenInfo[nzchar(sets) == TRUE]
    }
    
    return(x)
})

#' @rdname featureGroupsScreening-class
#' @export
setMethod("delete", "featureGroupsScreeningSet", doSFGroupsScreeningDelete)

#' @rdname featureGroupsScreening-class
#' @export
setMethod("annotateSuspects", "featureGroupsScreeningSet", function(fGroups, MSPeakLists, formulas, compounds,
                                                                    absMzDev = 0.005,
                                                                    specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
                                                                    checkFragments = c("mz", "formula", "compound"),
                                                                    formulasNormalizeScores = "max",
                                                                    compoundsNormalizeScores = "max",
                                                                    IDFile = system.file("misc", "IDLevelRules.yml",
                                                                                         package = "patRoon"),
                                                                    logPath = file.path("log", "ident"))
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertClass, . ~ MSPeakLists + formulas + compounds,
           c("MSPeakListsSet", "formulasSet", "compoundsSet"), null.ok = TRUE, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    unsetFGroups <- sapply(sets(fGroups), unset, obj = fGroups, simplify = FALSE)
    unsetMSPeakLists <- checkAndUnSetOther(sets(fGroups), MSPeakLists, "MSPeakLists", TRUE)
    unsetFormulas <- checkAndUnSetOther(sets(fGroups), formulas, "formulas", TRUE)
    unsetCompounds <- checkAndUnSetOther(sets(fGroups), compounds, "compounds", TRUE)
    
    logPath <- if (is.null(logPath)) rep(list(NULL), length(sets(fGroups))) else file.path(logPath, sets(fGroups))
    
    unsetFGroups <- Map(unsetFGroups, unsetMSPeakLists, unsetFormulas, unsetCompounds, logPath = logPath,
                        f = annotateSuspects, MoreArgs = list(absMzDev = absMzDev,
                                                              specSimParams = specSimParams,
                                                              checkFragments = checkFragments,
                                                              formulasNormalizeScores = formulasNormalizeScores,
                                                              compoundsNormalizeScores = compoundsNormalizeScores,
                                                              IDFile = IDFile))
    
    # re-generate so that annotation specific columns are added
    # NOTE: this will also clearout any previously added annotation columns
    fGroups@screenInfo <- mergeScreeningSetInfos(unsetFGroups)
    
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

    fGroups@screenInfo <- assignSetsIDLs(fGroups@screenInfo, mergedConsensusNames(fGroups))

    return(fGroups)
})

#' @rdname featureGroupsScreening-class
#' @export
setMethod("filter", "featureGroupsScreeningSet", function(obj, ..., onlyHits = NULL,
                                                          selectHitsBy = NULL, selectBestFGroups = FALSE,
                                                          maxLevel = NULL, maxFormRank = NULL, maxCompRank = NULL,
                                                          minAnnSimForm = NULL, minAnnSimComp = NULL, minAnnSimBoth = NULL,
                                                          absMinFragMatches = NULL, relMinFragMatches = NULL,
                                                          minRF = NULL, maxLC50 = NULL, negate = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertFlag, . ~ onlyHits + selectBestFGroups + negate, null.ok = c(TRUE, FALSE, FALSE), fixed = list(add = ac))
    checkmate::assertChoice(selectHitsBy, choices = c("intensity", "level"), null.ok = TRUE, add = ac)
    aapply(checkmate::assertCount, . ~ maxLevel + maxFormRank + maxCompRank + absMinFragMatches + relMinFragMatches,
           null.ok = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ minAnnSimForm + minAnnSimComp + minAnnSimBoth + minRF + maxLC50, null.ok = TRUE,
           fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    obj <- doSuspectFilter(obj, onlyHits = onlyHits, selectHitsBy, selectBestFGroups, maxLevel, maxFormRank, maxCompRank,
                           minAnnSimForm, minAnnSimComp, minAnnSimBoth, absMinFragMatches, relMinFragMatches, minRF,
                           maxLC50, negate)
    
    # filter functionality from fGroupsSet
    if (...length() > 0)
        obj <- callNextMethod(obj, ..., negate = negate)
    
    return(obj)
})

#' @rdname pred-quant
#' @export
setMethod("predictRespFactors", "featureGroupsScreeningSet", function(obj, calibrants, ...)
{
    if (length(obj) == 0)
        return(obj)
    
    checkmate::assertList(calibrants, types = "data.frame", any.missing = FALSE, len = length(sets(obj)))
    
    unsetFGroups <- sapply(sets(fGroups), unset, obj = fGroups, simplify = FALSE)
    unsetFGroups <- Map(unsetFGroups, calibrants, f = predictRespFactors, MoreArgs = list(...))
    obj <- syncScreeningSetObjects(obj, unsetFGroups)
    
    return(obj)
    
})

#' @rdname pred-tox
#' @export
setMethod("predictTox", "featureGroupsScreeningSet", doPredictToxSuspects)

#' @rdname pred-quant
#' @export
setMethod("calculateConcs", "featureGroupsScreeningSet", function(fGroups, featureAnn = NULL, areas = FALSE)
{
    # dummy method so that featureAnn can default to NULL
    callNextMethod(fGroups, featureAnn, areas)
})

#' @rdname pred-tox
#' @export
setMethod("calculateTox", "featureGroupsScreeningSet", function(fGroups, featureAnn = NULL)
{
    # dummy method so that featureAnn can default to NULL
    callNextMethod(fGroups, featureAnn)
})

setMethod("findMobilities", "featureGroupsScreeningSet", function(fGroups, mobPeaksParam, mzWindow = 0.005,
                                                                  IMSWindow = 0.01, clusterMethod = "distance",
                                                                  minIntensityIMS = 0, maxMSRTWindow = 2,
                                                                  chromPeaksParam = NULL, EICRTWindow = 20,
                                                                  peakRTWindow = 5, calcArea = "integrate",
                                                                  fallbackEIC = TRUE, parallel = TRUE,
                                                                  fromSuspects = FALSE, minMobilityMatches = 0)
{
    ac <- checkmate::makeAssertCollection()
    assertFindMobilitiesArgs(mobPeaksParam, mzWindow, IMSWindow, clusterMethod, minIntensityIMS, maxMSRTWindow,
                             chromPeaksParam, EICRTWindow, peakRTWindow, calcArea, fallbackEIC, parallel, ac)
    checkmate::assertFlag(fromSuspects, add = ac)
    checkmate::assertCount(minMobilityMatches, add = ac)
    checkmate::reportAssertions(ac)
    
    if (length(fGroups) == 0)
        return(fGroups) # nothing to do...

    anaInfo <- analysisInfo(fGroups)
    for (s in sets(fGroups))
    {
        anasSet <- anaInfo[set == s]$analysis
        if (fromSuspects)
        {
            fGroups@features <- assignFeatureMobilitiesSuspects(fGroups@features, screenInfo(unset(fGroups, s)),
                                                                IMSWindow,
                                                                \(ft, a) if (!a %chin% anasSet) ft[0] else ft)
        }
    }
    fGroups@features <- assignFeatureMobilitiesPeaks(fGroups@features, mobPeaksParam, mzWindow, IMSWindow, clusterMethod,
                                                     minIntensityIMS, maxMSRTWindow)
    fGroups@features <- reintegrateMobilityFeatures(fGroups@features, EICRTWindow, peakRTWindow, calcArea,
                                                    chromPeaksParam, fallbackEIC, parallel)
    fGroups <- clusterFGroupMobilities(fGroups, IMSWindow, TRUE)
    
    gInfo <- groupInfo(fGroups)
    scr <- expandAndUpdateScreenInfoForIMS(screenInfo(fGroups), gInfo)
    
    mySets <- sets(fGroups)
    fgSetNames <- sapply(mySets, function(s) names(fGroups[, sets = s]), simplify = FALSE)

    # Prune hits that are not actually in set: the data for mobility feature groups is copied from their parents, while
    # the actual mobility feature group may not be in the same sets.
    for (i in seq_len(nrow(scr)))
    {
        curSets <- unlist(strsplit(scr$sets[i], ",", fixed = TRUE))
        groupSets <- mySets[sapply(fgSetNames, function(n) scr$group[i] %chin% n)] # actual set presence of fGroup
        
        # update set assignment
        newSets <- intersect(curSets, groupSets)
        set(scr, i = i, j = "sets", paste0(newSets, collapse = ","))
        
        # set copied values from to NA from cleared out sets
        rmSets <- setdiff(curSets, newSets)
        NACols <- getAllMergedConsCols(names(scr), rmSets)
        set(scr, i = i, j = NACols, NA)
    }

    scr <- finalizeScreenInfoForIMS(scr, gInfo, minMobilityMatches, IMSWindow)
    fGroups@screenInfo <- scr
    
    return(fGroups)
})

#' @rdname suspect-screening
#' @export
setMethod("screenSuspects", "featureGroupsSet", function(fGroups, suspects, rtWindow, mzWindow, IMSWindow,
                                                         adduct, skipInvalid, prefCalcChemProps, neutralChemProps,
                                                         minMobilityMatches, onlyHits)
{
    verifyNoAdductIonizationArg(adduct)
    
    suspects <- assertAndPrepareSuspectsSets(suspects, sets(fGroups), skipInvalid)

    unsetFGroupsList <- sapply(sets(fGroups), unset, obj = fGroups, simplify = FALSE)
    setObjects <- Map(unsetFGroupsList, suspects,
                      f = function(fg, s) screenSuspects(fg, s, rtWindow = rtWindow, mzWindow = mzWindow,
                                                         IMSWindow = IMSWindow, adduct = NULL, skipInvalid = skipInvalid,
                                                         prefCalcChemProps = prefCalcChemProps,
                                                         neutralChemProps = neutralChemProps,
                                                         minMobilityMatches = minMobilityMatches, onlyHits = onlyHits))
    
    scr <- mergeScreeningSetInfos(setObjects)
    if (onlyHits)
        fGroups <- fGroups[, scr$group]
    
    return(featureGroupsScreeningSet(screenInfo = scr, groupAlgo = fGroups@groupAlgo, groupArgs = fGroups@groupArgs,
                                     groupVerbose = fGroups@groupVerbose, groups = copy(groupTable(fGroups)),
                                     groupInfo = copy(groupInfo(fGroups)), features = getFeatures(fGroups),
                                     ftindex = copy(groupFeatIndex(fGroups)),
                                     groupQualities = copy(groupQualities(fGroups)),
                                     groupScores = copy(groupScores(fGroups)), ISTDs = copy(internalStandards(fGroups)),
                                     ISTDAssignments = internalStandardAssignments(fGroups),
                                     annotations = copy(annotations(fGroups)),
                                     concentrations = copy(concentrations(fGroups)),
                                     toxicities = copy(toxicities(fGroups))))
})

#' @rdname suspect-screening
#' @export
setMethod("screenSuspects", "featureGroupsScreeningSet", doScreenSuspectsAmend)


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
    sInfo <- copy(screenInfo(obj))
    if (length(sInfo) > 0)
    {
        sInfo <- removeDTColumnsIfPresent(sInfo, c("formRank", "compRank", "estIDLevel", "sets"))
        # restore set specific columns
        setnames(sInfo, sub(paste0("\\-", set, "$"), "", names(sInfo)))
    }
    
    ret <- featureGroupsSetScreeningUnset(screenInfo = sInfo, groups = groupTable(uobj),
                                          groupInfo = groupInfo(uobj), features = getFeatures(uobj),
                                          ftindex = groupFeatIndex(uobj), annotations = annotations(uobj),
                                          groupQualities = groupQualities(uobj), groupScores = groupScores(uobj),
                                          ISTDs = internalStandards(uobj),
                                          ISTDAssignments = internalStandardAssignments(uobj),
                                          concentrations = concentrations(uobj), toxicities = toxicities(uobj))
    # override after constructing: parent constructor already sets algorithm,
    # which results in error about double assignment
    ret@algorithm <- paste0(algorithm(obj), "_unset")
    return(ret)
})
