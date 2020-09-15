#' @include main.R
#' @include feature_groups.R
#' @include utils-screening.R
NULL

#' @rdname suspect-screening
featureGroupsScreening <- setClass("featureGroupsScreening",
                                   slots = c(screenInfo = "data.table"),
                                   contains = "featureGroups")

setMethod("initialize", "featureGroupsScreening",
          function(.Object, ...) callNextMethod(.Object, algorithm = "screening", ...))


setMethod("screenInfo", "featureGroupsScreening", function(obj) obj@screenInfo)

setMethod("[", c("featureGroupsScreening", "ANY", "ANY", "missing"), function(x, i, j, ..., drop = TRUE)
{
    x <- callNextMethod()
    x@screenInfo <- x@screenInfo[group %in% names(x)]
    return(x)
})

setMethod("as.data.table", "featureGroupsScreening",
          function(x, ..., collapseSuspects = ",", onlyHits = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertString(collapseSuspects, null.ok = TRUE, add = ac)
    checkmate::assertFlag(onlyHits, add = ac)
    checkmate::reportAssertions(ac)
    
    ret <- callNextMethod(x, ...)
    if (nrow(ret) > 0)
    {
        si <- copy(screenInfo(x))
        setnames(si, c("rt", "mz"), c("susp_rt", "susp_mz"))
        
        if (!is.null(collapseSuspects))
        {
            si[, name := paste0(name, collapse = collapseSuspects), by = "group"]
            si <- unique(si, by = "group")
        }
        
        ret <- merge(ret, si, by = "group", all.x = !onlyHits, sort = FALSE)
    }
    return(ret)
})

#' @templateVar normParam compoundsNormalizeScores,formulasNormalizeScores
#' @templateVar noNone TRUE
#' @template norm-args
setMethod("annotateSuspects", "featureGroupsScreening", function(fGroups, MSPeakLists, formulas, compounds,
                                                                 absMzDev = 0.005, relMinMSMSIntensity = 0.05,
                                                                 checkFragments = c("mz", "formula", "compound"),
                                                                 formulasNormalizeScores = "max",
                                                                 compoundsNormalizeScores = "max",
                                                                 IDLevelRules = defaultIDLevelRules())
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertClass, . ~ MSPeakLists + formulas + compounds,
           c("MSPeakLists", "formulas", "compounds"), null.ok = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ absMzDev + relMinMSMSIntensity, lower = 0,
           finite = TRUE, fixed = list(add = ac))
    checkmate::assertSubset(checkFragments, c("mz", "formula", "compound"), add = ac)
    aapply(assertNormalizationMethod, . ~ formulasNormalizeScores + compoundsNormalizeScores, withNone = FALSE,
           fixed = list(add = ac))
    checkmate::assertDataFrame(IDLevelRules, types = c("numeric", "character", "logical"),
                               all.missing = TRUE, min.rows = 1, add = ac)
    assertHasNames(IDLevelRules,
                   c("level", "subLevel", "type", "score", "relative", "value", "higherThanNext", "mustExist"),
                   add = ac)
    checkmate::reportAssertions(ac)
    
    hash <- makeHash(fGroups, MSPeakLists, formulas, compounds, absMzDev,
                     relMinMSMSIntensity, checkFragments, formulasNormalizeScores,
                     compoundsNormalizeScores, IDLevelRules)
    cd <- loadCacheData("annotateSuspects", hash)
    if (!is.null(cd))
        return(cd)
    
    si <- copy(screenInfo(fGroups))
    
    for (i in seq_len(nrow(si)))
    {
        gName <- si$group[i]
        MSMSList <- if (!is.null(MSPeakLists)) MSPeakLists[[gName]][["MSMS"]] else NULL
        fTable <- if (!is.null(formulas)) formulas[[gName]] else NULL
        fScRanges <- if (!is.null(formulas)) formulas@scoreRanges[[gName]] else NULL
        cTable <- if (!is.null(compounds)) compounds[[gName]] else NULL
        cScRanges <- if (!is.null(compounds)) compounds@scoreRanges[[gName]] else NULL
        
        suspFormRank <- NA_integer_
        if (!is.null(fTable) && !is.null(si[["formula"]]) && !is.na(si$formula[i]))
        {
            unFTable <- unique(fTable, by = "formula")
            suspFormRank <- which(si$formula[i] == unFTable$neutral_formula)
            suspFormRank <- if (length(suspFormRank) > 0) suspFormRank[1] else NA_integer_
        }
        
        suspIK1 <- if (!is.null(si[["InChIKey"]]) && !is.na(si$InChIKey[i])) getIKBlock1(si$InChIKey[i]) else NULL
        annSim <- 0; suspCompRank <- NA_integer_
        if (!is.null(MSMSList) && !is.null(cTable) && !is.null(suspIK1))
        {
            suspCompRank <- which(suspIK1 == cTable$InChIKey1)
            suspCompRank <- if (length(suspCompRank) > 0) suspCompRank[1] else NA_integer_
            
            if (!is.na(suspCompRank) && !is.null(cTable[["fragInfo"]][[suspCompRank]]))
                annSim <- annotatedMSMSSimilarity(cTable[["fragInfo"]][[suspCompRank]],
                                                  MSMSList, absMzDev, relMinMSMSIntensity)
        }
        
        set(si, i, c("suspFormRank", "suspCompRank", "annotatedMSMSSimilarity"), list(suspFormRank, suspCompRank, annSim))
        set(si, i, "estIDLevel",
            estimateIdentificationLevel(si$d_rt[i], suspIK1, si$formula[i], annSim,
                                        if (!is.null(si[["fragments_mz"]])) si$fragments_mz[i] else NULL,
                                        if (!is.null(si[["fragments_formula"]])) si$fragments_formula[i] else NULL,
                                        checkFragments, MSMSList, fTable, fScRanges,
                                        formulasNormalizeScores, cTable,
                                        mCompNames = if (!is.null(compounds)) mergedCompoundNames(compounds) else NULL,
                                        cScRanges, compoundsNormalizeScores, absMzDev, IDLevelRules))
    }
    
    fGroups@screenInfo <- si
    
    saveCacheData("annotateSuspects", fGroups, hash)
    
    return(fGroups)
})

setMethod("filter", "featureGroupsScreening", function(obj, ..., onlyHits = FALSE,
                                                       selectBy = NULL, maxLevel = NULL, negate = FALSE)
{
    # UNDONE: doc that selectBy only applies to hits, in case of ties: first hit
    # UNDONE: mention that filter with onlyHits may need to be repeated
    
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertFlag, . ~ onlyHits + negate, fixed = list(add = ac))
    checkmate::assertChoice(selectBy, c("intensity", "level"), null.ok = TRUE, add = ac)
    checkmate::assertInt(maxLevel, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (!is.null(selectBy))
    {
        gTab <- as.data.table(obj, onlyHits = TRUE)
        doKeep <- function(v, d) is.na(v) | length(v) == 1 | order(v, decreasing = d) == 1
        if (selectBy == "intensity")
        {
            gTab[, avgInts := rowMeans(.SD), .SDcol = analyses(obj)]
            gTab <- gTab[, keep := doKeep(avgInts, !negate), by = "name"]
        }
        else # keep best ID level
        {
            if (is.null(gTab[["estIDLevel"]]))
                stop("Cannot select by identification level: no annotation data available (did you run annotateSuspects()?). ")
            gTab <- gTab[, keep := doKeep(estIDLevel, negate), by = "name"]
        }
        
        if (any(!gTab$keep))
        {
            # merge-in keep column so we can subset screenInfo
            si <- copy(screenInfo(obj))
            si[gTab, keep := i.keep, on = c("group", "name")]
            setorderv(si, "name")
            obj@screenInfo <- si[keep == TRUE, -"keep"]
        }
    }
    
    if (!is.null(maxLevel) && !is.null(screenInfo(obj)[["estIDLevel"]]))
    {
        pred <- if (negate) function(x) x > maxLevel else function(x) x <= maxLevel
        obj@screenInfo <- screenInfo(obj)[nzchar(estIDLevel) & pred(numericIDLevel(estIDLevel))]
    }

    # NOTE: do last in case previous steps removed hits 
    if (onlyHits)
    {
        sGroups <- unique(screenInfo(obj)$group)
        if (negate)
            obj <- obj[, setdiff(names(obj), sGroups)]
        else
            obj <- obj[, sGroups]
    }
    
    if (...length() > 0)
        obj <- callNextMethod(obj, ..., negate)
    
    return(obj)
})

#' @details \code{groupFeaturesScreening} uses results from \code{screenSuspects}
#'   to transform an existing \code{\link{featureGroups}} object by (1) renaming
#'   any matched feature groups by the respective name of the suspect and (2)
#'   filtering out any feature groups that were not matched. A common workflow
#'   is to first obtain and group features (with \emph{e.g.}
#'   \code{\link{findFeatures}} and \code{\link{groupFeatures}}), screen them
#'   with \code{screenSuspects}, convert the \code{featureGroups} object that was
#'   used for screening with this method function and continue any further
#'   workflow steps such as compound identification as with 'regular'
#'   \code{featureGroups}.
#'
#' @param fGroups The \code{\link{featureGroups}} object that should be
#'   transformed (and was used to obtain the screening results).
#' @param scr The screening results table returned by \code{screenSuspects}.
#'
#' @return \code{groupFeaturesScreening} returns a modified \code{featureGroups}
#'   object in which those feature groups that were not matched by any suspects
#'   are removed and others renamed by the respective suspect name. In case of
#'   duplicate suspect results, feature group names are made unique with
#'   \code{\link{make.unique}}.
#'
#' @note Please note that \code{groupFeaturesScreening} method can only
#'   transform the \code{featureGroups} object that was used to obtain the given
#'   screening results.
#'
#' @rdname suspect-screening
#' @aliases groupFeaturesScreening
#' @export
setMethod("groupFeaturesScreening", "featureGroups", function(fGroups, suspects, rtWindow, mzWindow,
                                                              adduct, skipInvalid)
{
    if (!is.null(adduct))
        adduct <- checkAndToAdduct(adduct)
    
    checkmate::assertFlag(skipInvalid) # not in assert collection, should fail before assertSuspectList
    
    ac <- checkmate::makeAssertCollection()
    assertSuspectList(suspects, adduct, skipInvalid, add = ac)
    aapply(checkmate::assertNumber, . ~ rtWindow + mzWindow, lower = 0, finite = TRUE, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    # do this before checking cache to ensure proper errors/warnings are thrown!
    suspects <- prepareSuspectList(suspects, adduct, skipInvalid)
    
    hash <- makeHash(fGroups, suspects, rtWindow, mzWindow, adduct, skipInvalid)
    cd <- loadCacheData("screenSuspects", hash)
    if (!is.null(cd))
        return(cd)

    scr <- doScreenSuspects(fGroups, suspects, rtWindow, mzWindow, adduct, skipInvalid)    
    
    ret <- featureGroupsScreening(screenInfo = scr, groups = copy(groupTable(fGroups)),
                                  analysisInfo = analysisInfo(fGroups), groupInfo = groupInfo(fGroups),
                                  features = getFeatures(fGroups), ftindex = copy(groupFeatIndex(fGroups)))
    
    saveCacheData("screenSuspects", ret, hash)
    
    return(ret)
})
