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

setMethod("[", c("featureGroupsScreening", "ANY", "ANY", "missing"), function(x, i, j, ..., rGroups,
                                                                              suspects = NULL, drop = TRUE)
{
    checkmate::assertCharacter(suspects, null.ok = TRUE)
    
    x <- callNextMethod(x, i, j, ..., rGroups = rGroups, drop = drop)
    
    if (!is.null(suspects))
        x <- x[, x@screenInfo[name %in% suspects]$group]
    
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
            # only keep unique and remove suspect specific columns
            # UNDONE: keep specific columns if only one suspect?
            si <- unique(si[, c("group", "name"), with = FALSE], by = "group")
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
                                                                 simMSMSMethod = "cosine",
                                                                 checkFragments = c("mz", "formula", "compound"),
                                                                 formulasNormalizeScores = "max",
                                                                 compoundsNormalizeScores = "max",
                                                                 IDFile = system.file("inst", "misc", "IDLevelRules.yml",
                                                                                      package = "patRoon"))
{
    # UNDONE: prog bar
    # UNDONE/document: annSimBoth falls back to annSimForm/annSimComp if no formulas available
    
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertClass, . ~ MSPeakLists + formulas + compounds,
           c("MSPeakLists", "formulas", "compounds"), null.ok = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ absMzDev + relMinMSMSIntensity, lower = 0,
           finite = TRUE, fixed = list(add = ac))
    checkmate::assertChoice(simMSMSMethod, c("cosine", "jaccard"))
    checkmate::assertSubset(checkFragments, c("mz", "formula", "compound"), add = ac)
    aapply(assertNormalizationMethod, . ~ formulasNormalizeScores + compoundsNormalizeScores, withNone = FALSE,
           fixed = list(add = ac))
    checkmate::assertFileExists(IDFile, "r", add = ac)
    checkmate::reportAssertions(ac)
    
    hash <- makeHash(fGroups, MSPeakLists, formulas, compounds, absMzDev,
                     relMinMSMSIntensity, simMSMSMethod, checkFragments, formulasNormalizeScores,
                     compoundsNormalizeScores, makeFileHash(IDFile))
    cd <- loadCacheData("annotateSuspects", hash)
    if (!is.null(cd))
        return(cd)
    
    IDLevelRules <- yaml::yaml.load_file(IDFile, eval.expr = FALSE)
    
    if (!checkmate::test_named(IDLevelRules))
        stop("No valid rules could be loaded")
    if (!all(grepl("^[[:digit:]]+[[:alpha:]]?$", names(IDLevelRules))))
        stop("Levels should be defined as a number and may optionally followed by one character (e.g. 3, 2b etc)")
    
    IDLevelRules <- IDLevelRules[order(names(IDLevelRules))] # sort to ensure lowest levels will be tested first
    
    
    mzWithin <- function(mz1, mz2) abs(mz1 - mz2) <= absMzDev
    
    si <- copy(screenInfo(fGroups))
    
    for (i in seq_len(nrow(si)))
    {
        gName <- si$group[i]
        MSMSList <- if (!is.null(MSPeakLists)) MSPeakLists[[gName]][["MSMS"]] else NULL
        fTable <- if (!is.null(formulas)) formulas[[gName]] else NULL
        fScRanges <- if (!is.null(formulas)) formulas@scoreRanges[[gName]] else NULL
        cTable <- if (!is.null(compounds)) compounds[[gName]] else NULL
        cScRanges <- if (!is.null(compounds)) compounds@scoreRanges[[gName]] else NULL
        
        suspFormRank <- NA_integer_; annSimForm <- annSimBoth <- NA_real_
        if (!is.null(fTable) && !is.null(si[["formula"]]) && !is.na(si$formula[i]))
        {
            unFTable <- unique(fTable, by = "formula")
            suspFormRank <- which(si$formula[i] == unFTable$neutral_formula)
            suspFormRank <- if (length(suspFormRank) > 0) suspFormRank[1] else NA_integer_
            if (!is.na(suspFormRank))
                annSimForm <- annSimBoth <- annotatedMSMSSimilarity(annotatedPeakList(formulas,
                                                                                      precursor = unFTable$formula[suspFormRank],
                                                                                      groupName = gName, MSPeakLists = MSPeakLists),
                                                                    absMzDev, relMinMSMSIntensity, simMSMSMethod)
        }
        
        suspIK1 <- if (!is.null(si[["InChIKey"]]) && !is.na(si$InChIKey[i])) getIKBlock1(si$InChIKey[i]) else NULL
        annSimComp <- NA_real_; suspCompRank <- NA_integer_
        if (!is.null(MSMSList) && !is.null(cTable) && !is.null(suspIK1))
        {
            suspCompRank <- which(suspIK1 == cTable$InChIKey1)
            suspCompRank <- if (length(suspCompRank) > 0) suspCompRank[1] else NA_integer_
            
            if (!is.na(suspCompRank) && !is.null(cTable[["fragInfo"]][[suspCompRank]]))
            {
                annSimComp <- annotatedMSMSSimilarity(annotatedPeakList(compounds, index = suspCompRank,
                                                                        groupName = gName, MSPeakLists = MSPeakLists),
                                                      absMzDev, relMinMSMSIntensity, simMSMSMethod)
                
                if (!is.na(suspFormRank))
                    annSimBoth <- annotatedMSMSSimilarity(annotatedPeakList(compounds, index = suspCompRank,
                                                                            groupName = gName, MSPeakLists = MSPeakLists,
                                                                            formulas = formulas),
                                                          absMzDev, relMinMSMSIntensity, simMSMSMethod)
                else
                    annSimBoth <- annSimComp
            }
        }
        
        fragMZMatches <- fragFormMatches <- fragFormCompMatches <- NA_integer_
        fragMZs <- fragForms <- NULL
        maxSuspFrags <- maxFragMatches <- NA_integer_
        if (!is.null(MSMSList) && !is.null(si[["fragments_mz"]]) &&
            !is.na(si[["fragments_mz"]][i]) && nzchar(si[["fragments_mz"]][i]) &&
            "mz" %in% checkFragments)
        {
            fragMZs <- as.numeric(unlist(strsplit(si[["fragments_mz"]][i], ";")))
            maxSuspFrags <- length(fragMZs)
            maxFragMatches <- sum(sapply(MSMSList$mz, function(mz1) any(sapply(fragMZs, mzWithin, mz1 = mz1))))
        }
        if (!is.null(si[["fragments_formula"]]) && !is.na(si[["fragments_formula"]][i]) &&
            nzchar(si[["fragments_formula"]][i]))
        {
            fragForms <- unlist(strsplit(si[["fragments_formula"]][i], ";"))
            maxSuspFrags <- max(NAToZero(maxSuspFrags), length(fragForms))
            
            if (!is.null(fTable) && "formula" %in% checkFragments)
            {
                frTable <- fTable[byMSMS == TRUE & si$formula[i] == neutral_formula]
                if (nrow(frTable) > 0)
                {
                    fi <- getFragmentInfoFromForms(MSMSList, frTable)
                    maxFragMatches <- max(NAToZero(maxFragMatches), sum(fragForms %in% fi$formula))
                }
            }
            
            if (!is.null(cTable) && "compound" %in% checkFragments && !is.na(suspCompRank) &&
                !is.null(cTable[["fragInfo"]][[suspCompRank]]))
                maxFragMatches <- max(NAToZero(maxFragMatches), sum(fragForms %in% cTable[["fragInfo"]][[suspCompRank]]$formula))
        }

        estIDLevel <- estimateIdentificationLevel(si$name[i], si$group[i], si$d_rt[i], suspIK1, si$formula[i],
                                                  annSimForm, annSimComp, annSimBoth,
                                                  maxSuspFrags, maxFragMatches, fTable, suspFormRank, fScRanges,
                                                  formulasNormalizeScores, cTable, suspCompRank,
                                                  mCompNames = if (!is.null(compounds)) mergedCompoundNames(compounds) else NULL,
                                                  cScRanges, compoundsNormalizeScores, absMzDev, IDLevelRules)
        
        set(si, i,
            c("suspFormRank", "suspCompRank", "annSimForm", "annSimComp", "annSimBoth",
              "maxFrags", "maxFragMatches", "estIDLevel"),
            list(suspFormRank, suspCompRank, annSimForm, annSimComp, annSimBoth, maxSuspFrags, maxFragMatches, estIDLevel))
    }
    
    rmCols <- c("suspFormRank", "suspCompRank", "annSimForm", "annSimComp", "annSimBoth", "maxFrags", "maxFragMatches")
    rmCols <- rmCols[sapply(rmCols, function(col) !is.null(si[[col]]) && all(is.na(si[[col]])))]
    if (length(rmCols) > 0)
        si <- si[, setdiff(names(si), rmCols), with = FALSE]
    
    fGroups@screenInfo <- si
    
    saveCacheData("annotateSuspects", fGroups, hash)
    
    return(fGroups)
})

setMethod("filter", "featureGroupsScreening", function(obj, ..., onlyHits = NULL,
                                                       selectHitsBy = NULL, selectBestFGroups = FALSE,
                                                       maxLevel = NULL, maxFormRank = NULL, maxCompRank = NULL,
                                                       minAnnSimForm = NULL, minAnnSimComp = NULL, minAnnSimBoth = NULL,
                                                       minFragMatches = NULL, negate = FALSE)
{
    # UNDONE: doc that selectHitsBy/selectFGroupsBy only applies to hits, in case of ties: first hit
    # UNDONE: mention that filter with onlyHits may need to be repeated
    # UNDONE: cache?
    # UNDONE: minFragMatches --> split in abs/rel thresholds?
    # UNDONE: keep or remove NA values with colFilter()? document what happens
    # UNDONE: properly document negate:
    #   - selectHitsBy: select worst hit
    #   - onlyHits: if TRUE only select non-hits, if FALSE select hits, if NULL nothing
    
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertFlag, . ~ onlyHits + selectBestFGroups + negate, null.ok = c(TRUE, FALSE, FALSE), fixed = list(add = ac))
    checkmate::assertChoice(selectHitsBy, choices = c("intensity", "level"), null.ok = TRUE, add = ac)
    aapply(checkmate::assertCount, . ~ maxLevel + maxFormRank + maxCompRank + minFragMatches, null.ok = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ minAnnSimForm + minAnnSimComp + minAnnSimBoth, null.ok = TRUE, fixed = list(add = ac))
    checkmate::reportAssertions(ac)

    if (nrow(screenInfo(obj)) > 0)
    {
        colFilter <- function(pred, what, col)
        {
            val <- get(what)
            if (!is.null(val))
            {
                if (is.null(screenInfo(obj)[[col]]))
                    warning(sprintf("Cannot apply %s filter: no annotation data available (did you run annotateSuspects()?).", what))
                else
                {
                    if (negate)
                        doPred <- function(x, v) is.na(x) | !nzchar(x) | !pred(x, v)
                    else
                        doPred <- function(x, v) !is.na(x) & nzchar(x) & pred(x, v)
                    obj@screenInfo <- screenInfo(obj)[doPred(get(col), val)]
                }
            }
            return(obj)
        }
        minPred <- function(x, v) x >= v
        maxPred <- function(x, v) x <= v
        levPred <- function(x, v) maxPred(numericIDLevel(x), v)
        
        obj <- colFilter(levPred, "maxLevel", "estIDLevel")
        obj <- colFilter(maxPred, "maxFormRank", "suspFormRank")
        obj <- colFilter(maxPred, "maxCompRank", "suspCompRank")
        obj <- colFilter(minPred, "minAnnSimForm", "annSimForm")
        obj <- colFilter(minPred, "minAnnSimComp", "annSimComp")
        obj <- colFilter(minPred, "minAnnSimBoth", "annSimBoth")
        obj <- colFilter(minPred, "minFragMatches", "maxFragMatches")
        
        # do here so that only duplicates not yet filtered out in previous steps are considered
        if (!is.null(selectHitsBy) || selectBestFGroups)
        {
            doKeep <- function(v, d) is.na(v) | length(v) == 1 | seq_along(v) == order(v, decreasing = d)[1]
            doSelectFilter <- function(si, by, byCol)
            {
                if (by == "level" && is.null(si[["estIDLevel"]]))
                    warning("Cannot select by identification level: no annotation data available (did you run annotateSuspects()?).")
                else
                {
                    gTab <- as.data.table(obj, collapseSuspects = NULL, onlyHits = TRUE)
                    
                    if (by == "intensity")
                    {
                        gTab[, avgInts := rowMeans(.SD), .SDcol = analyses(obj)]
                        gTab <- gTab[, keep := doKeep(avgInts, !negate), by = byCol]
                    }
                    else # select by best hit
                        gTab <- gTab[, keep := doKeep(estIDLevel, negate), by = byCol]
                    
                    if (any(!gTab$keep))
                    {
                        # merge-in keep column so we can subset screenInfo
                        si <- copy(si)
                        si[gTab, keep := i.keep, on = c("group", "name")]
                        setorderv(si, "name")
                        obj@screenInfo <- si[keep == TRUE, -"keep"]
                    }
                }
                return(obj@screenInfo)
            }
            
            if (!is.null(selectHitsBy))
                obj@screenInfo <- doSelectFilter(obj@screenInfo, selectHitsBy, "name")
            if (selectBestFGroups)
                obj@screenInfo <- doSelectFilter(obj@screenInfo, "level", "group")
        }
    }
    
    # NOTE: do last in case previous steps removed hits 
    if (!is.null(onlyHits))
    {
        sGroups <- unique(screenInfo(obj)$group)
        if (negate && onlyHits)
            obj <- obj[, setdiff(names(obj), sGroups)]
        else
            obj <- obj[, sGroups]
    }
    
    if (...length() > 0)
        obj <- callNextMethod(obj, ..., negate = negate)
    
    return(obj)
})

#' @details \code{screenSuspects} uses results from \code{screenSuspects}
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
#' @return \code{screenSuspects} returns a modified \code{featureGroups}
#'   object in which those feature groups that were not matched by any suspects
#'   are removed and others renamed by the respective suspect name. In case of
#'   duplicate suspect results, feature group names are made unique with
#'   \code{\link{make.unique}}.
#'
#' @note Please note that \code{screenSuspects} method can only
#'   transform the \code{featureGroups} object that was used to obtain the given
#'   screening results.
#'
#' @rdname suspect-screening
#' @export
setMethod("screenSuspects", "featureGroups", function(fGroups, suspects, rtWindow, mzWindow,
                                                      adduct, skipInvalid, onlyHits)
{
    if (!is.null(adduct))
        adduct <- checkAndToAdduct(adduct)
    
    checkmate::assertFlag(skipInvalid) # not in assert collection, should fail before assertSuspectList
    
    ac <- checkmate::makeAssertCollection()
    assertSuspectList(suspects, adduct, skipInvalid, add = ac)
    aapply(checkmate::assertNumber, . ~ rtWindow + mzWindow, lower = 0, finite = TRUE, fixed = list(add = ac))
    checkmate::assertFlag(onlyHits, add = ac)
    checkmate::reportAssertions(ac)
    
    # do this before checking cache to ensure proper errors/warnings are thrown!
    suspects <- prepareSuspectList(suspects, adduct, skipInvalid)
    
    hash <- makeHash(fGroups, suspects, rtWindow, mzWindow, adduct, skipInvalid)
    cd <- loadCacheData("screenSuspects", hash)
    if (!is.null(cd))
        return(cd)

    scr <- doScreenSuspects(fGroups, suspects, rtWindow, mzWindow, adduct, skipInvalid)

    if (onlyHits)
        fGroups <- fGroups[, scr$group]
    
    ret <- featureGroupsScreening(screenInfo = scr, groups = copy(groupTable(fGroups)),
                                  analysisInfo = analysisInfo(fGroups), groupInfo = groupInfo(fGroups),
                                  features = getFeatures(fGroups), ftindex = copy(groupFeatIndex(fGroups)))
    
    saveCacheData("screenSuspects", ret, hash)
    
    return(ret)
})
