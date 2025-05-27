# SPDX-FileCopyrightText: 2016 - 2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
#' @include feature_groups.R
#' @include feature_groups-set.R
#' @include utils-screening.R
NULL

#' Class for suspect screened feature groups.
#'
#' This class derives from \code{\link{featureGroups}} and adds suspect screening information.
#'
#' @param obj,object,x The \code{featureGroupsScreening} object.
#' @param i,j,reorder See \code{\link{featureGroups}}.
#' @param \dots Further arguments passed to the base method.
#' @param onlyHits If \itemize{
#'
#'   \item \code{negate=FALSE} and \code{onlyHits=TRUE} then all feature groups without suspect hits will be removed.
#'   Otherwise nothing will be done.
#'
#'   \item \code{negate=TRUE} then \code{onlyHits=TRUE} will select feature groups without suspect hits,
#'   \code{onlyHits=FALSE} will only retain feature groups with suspect matches and this filter is ignored if
#'   \code{onlyHits=NULL}.
#'
#'   }
#'
#' @slot screenInfo A (\code{\link{data.table}}) with results from suspect screening. This table will be amended with
#'   ID confidence data when \code{\link{estimateIDConfidence}} is run.
#' 
#' @template ms2q-slot
#'
#' @templateVar class featureGroupsScreening
#' @template class-hierarchy
#'
#' @seealso \code{\link{featureGroups}}
#'
#' @export
featureGroupsScreening <- setClass("featureGroupsScreening",
                                   slots = c(screenInfo = "data.table", MS2QuantMeta = "list"),
                                   contains = "featureGroups")

setMethod("initialize", "featureGroupsScreening",
          function(.Object, ...) callNextMethod(.Object, algorithm = "screening", ...))

setMethod("clearMobilities", "featureGroupsScreening", doFGroupScreeningClearMobilities)

#' @describeIn featureGroupsScreening Returns a table with screening information
#'   (see \code{screenInfo} slot).
#' @aliases screenInfo
#' @export
setMethod("screenInfo", "featureGroupsScreening", function(obj) obj@screenInfo)

#' @describeIn featureGroupsScreening Shows summary information for this object.
#' @export
setMethod("show", "featureGroupsScreening", function(object)
{
    callNextMethod(object)
    doScreeningShow(object)
})

#' @describeIn featureGroupsScreening Subset on analyses, feature groups and/or
#'   suspects.
#'   
#' @param suspects An optional \code{character} vector with suspect names. If
#'   specified, only \code{featureGroups} will be kept that are assigned to
#'   these suspects.
#' @param drop Ignored.
#'
#' @export
setMethod("[", c("featureGroupsScreening", "ANY", "ANY", "missing"), function(x, i, j, ..., suspects = NULL,
                                                                              reorder = FALSE, drop = TRUE)
{
    checkmate::assertCharacter(suspects, null.ok = TRUE)
    
    x <- callNextMethod(x, i, j, ..., reorder = reorder, drop = drop)
    
    if (!is.null(suspects))
        x <- x[, x@screenInfo[name %in% suspects]$group]
    
    return(x)
})

#' @templateVar where featureGroupsScreening
#' @templateVar what feature groups or screening results
#' @template delete
#' @param k The \code{k} argument is used to delete screening results (instead of features) and should be: \itemize{
#'
#'   \item a \code{character} vector with suspect names that should be removed
#'
#'   \item a \code{function} that is called with the screening info table and should return a \code{logical} vector for
#'   each suspect row to be removed
#'
#'   \item \code{NA} to remove all screening results, which is especially useful when paired with the \code{j} argument,
#'   \emph{i.e.} to remove all screening results for a particular set of feature groups.
#'
#'   \item \code{NULL} to not touch screening results and only perform deletion as the \code{\link{featureGroups}}
#'   method.
#'
#'   }
#'
#'   Setting both \code{i} and \code{k} is currently not supported.
#'
#' @export
setMethod("delete", "featureGroupsScreening", doSFGroupsScreeningDelete)

#' @rdname pred-quant
#' @export
setMethod("predictRespFactors", "featureGroupsScreening", function(obj, calibrants, eluent, organicModifier, pHAq,
                                                                   concUnit = "ugL", calibConcUnit = concUnit)
{
    checkPackage("MS2Quant", "kruvelab/MS2Quant")
    
    ac <- checkmate::makeAssertCollection()
    assertQuantEluent(eluent, obj, add = ac)
    checkmate::assertChoice(organicModifier, c("MeOH", "MeCN"), add = ac)
    checkmate::assertNumber(pHAq, finite = TRUE, add = ac)
    aapply(assertConcUnit, . ~ concUnit + calibConcUnit, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    calibrants <- assertAndPrepareQuantCalib(calibrants, calibConcUnit)

    if (length(obj) == 0)
        return(obj)
    
    scr <- screenInfo(obj)
    if (is.null(scr[["SMILES"]]) || all(is.na(scr$SMILES)))
        stop("Suspects lack necessary SMILES information to perform calculations, aborting...", call. = FALSE)
    if (any(is.na(scr$SMILES)))
        warning("Some suspect SMILES are NA and will be ignored", call. = FALSE)
    
    inp <- screenInfo(obj)[, c("group", "SMILES"), with = FALSE]
    inp <- inp[!is.na(SMILES)]
    # avoid duplicate calculations if there happen to be suspects with the same SMILES
    inp <- unique(inp, by = c("group", "SMILES"))
    
    printf("Predicting response factors from SMILES with MS2Quant for %d suspects...\n", nrow(inp))
    res <- predictRespFactorsSMILES(inp, groupInfo(obj), calibrants, eluent, organicModifier, pHAq, concUnit)
    
    if (!is.null(scr[["RF_SMILES"]]))
        scr[, RF_SMILES := NULL] # clearout for merge below
    scr <- merge(scr, res$RFs[, c("group", "SMILES", "RF_SMILES"), with = FALSE], by = c("group", "SMILES"),
                 sort = FALSE, all.x = TRUE)
    
    obj@screenInfo <- scr
    obj@MS2QuantMeta <- res$MD
    
    return(obj)
})

#' @rdname pred-tox
#' @export
setMethod("predictTox", "featureGroupsScreening", doPredictToxSuspects)

#' @rdname pred-quant
#' @export
setMethod("calculateConcs", "featureGroupsScreening", function(fGroups, featureAnn = NULL, areas = FALSE)
{
    checkmate::assertFlag(areas)
    
    if (!is.null(featureAnn) && length(featureAnn) > 0)
        fGroups <- callNextMethod()
    else
        fGroups@concentrations <- data.table()

    scr <- screenInfo(fGroups)
    
    if (is.null(scr[["RF_SMILES"]]))
    {
        cat("Screening results lacks predicted response factors and will not be used for quantitation.",
            "You can use predictRespFactors() to add suspect response factors.\n")
        return(fGroups)
    }
    
    resp <- scr[, c("group", "name", "SMILES", "RF_SMILES"), with = FALSE]
    resp[, type := "suspect"]
    setnames(resp, c("SMILES", "name", "RF_SMILES"), c("candidate", "candidate_name", "RF"))
    
    concs <- calcFeatureConcs(fGroups, resp, areas)

    fGroups@concentrations <- finalizeFeaturePredTab(rbind(fGroups@concentrations, concs, fill = TRUE))
    
    return(fGroups)
})

#' @rdname pred-tox
#' @export
setMethod("calculateTox", "featureGroupsScreening", function(fGroups, featureAnn = NULL)
{
    if (!is.null(featureAnn) && length(featureAnn) > 0)
        fGroups <- callNextMethod()
    else
        fGroups@toxicities <- data.table()
    
    scr <- screenInfo(fGroups)
    
    if (is.null(scr[["LC50_SMILES"]]))
    {
        cat("Screening results lacks predicted toxicity values and will not be used.",
            "You can use predictTox() to add suspect toxicity values.\n")
        return(fGroups)
    }
    
    LC50s <- scr[, c("group", "name", "SMILES", "LC50_SMILES"), with = FALSE]
    LC50s[, type := "suspect"]
    setnames(LC50s, c("SMILES", "name", "LC50_SMILES"), c("candidate", "candidate_name", "LC50"))
    
    fGroups@toxicities <- finalizeFeaturePredTab(rbind(fGroups@toxicities, LC50s, fill = TRUE))
    
    return(fGroups)
})

#' @param checkFragments Which type(s) of MS/MS fragments from workflow data should be checked to evaluate the number of
#'   suspect fragment matches (\emph{i.e.} from the \code{fragments_mz}/\code{fragments_formula} columns in the suspect
#'   list). Valid values are: \code{"mz"}, \code{"formula"}, \code{"compounds"}. The former uses \emph{m/z} values in
#'   the specified \code{MSPeakLists} object, whereas the others use the formulae that were annotated to MS/MS peaks in
#'   the given \code{formulas} or \code{compounds} objects. Multiple values are possible: in this case the maximum
#'   number of fragment matches will be reported.
#' @rdname id-conf
#' @export
setMethod("estimateIDConfidence", "featureGroupsScreening", function(obj, MSPeakLists = NULL, formulas = NULL,
                                                                     compounds = NULL,
                                                                     absMzDev = defaultLim("mz", "medium"),
                                                                     checkFragments = c("mz", "formula", "compound"),
                                                                     formulasNormalizeScores = "max",
                                                                     compoundsNormalizeScores = "max",
                                                                     IDFile = system.file("misc", "IDLevelRules.yml",
                                                                                          package = "patRoon"),
                                                                     logPath = file.path("log", "ident"))
{
    # NOTE: keep args in sync with sets method
    
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertClass, . ~ MSPeakLists + formulas + compounds,
           c("MSPeakLists", "formulas", "compounds"), null.ok = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ absMzDev, lower = 0,
           finite = TRUE, fixed = list(add = ac))
    checkmate::assertSubset(checkFragments, c("mz", "formula", "compound"), add = ac)
    aapply(assertNormalizationMethod, . ~ formulasNormalizeScores + compoundsNormalizeScores, withNone = FALSE,
           fixed = list(add = ac))
    checkmate::assertFileExists(IDFile, "r", add = ac)
    if (!is.null(logPath))
        assertCanCreateDir(logPath, add = ac)
    checkmate::reportAssertions(ac)

    hash <- makeHash(obj, MSPeakLists, formulas, compounds, absMzDev, checkFragments, formulasNormalizeScores,
                     compoundsNormalizeScores, makeFileHash(IDFile))
    cd <- loadCacheData("estimateIDConfidence", hash)
    if (!is.null(cd))
        return(cd)
    
    IDLevelRules <- readIDLRules(IDFile)
    
    if (nrow(screenInfo(obj)) == 0)
    {
        cat("No suspect hits, nothing to annotate")
        return(obj)
    }
        
        
    mzWithin <- function(mz1, mz2) abs(mz1 - mz2) <= absMzDev

    si <- copy(screenInfo(obj))
    annCols <- suspAnnCols()
    si <- si[, setdiff(names(si), annCols), with = FALSE] # remove any previous results
    
    printf("Annotating %d suspects...\n", nrow(si))
    prog <- openProgBar(0, nrow(si))
    
    for (i in seq_len(nrow(si)))
    {
        gName <- si$group[i]
        MSMSList <- if (!is.null(MSPeakLists)) MSPeakLists[[gName]][["MSMS"]] else NULL
        fTable <- if (!is.null(formulas)) formulas[[gName]] else NULL
        fScRanges <- if (!is.null(formulas)) formulas@scoreRanges[[gName]] else NULL
        cTable <- if (!is.null(compounds)) compounds[[gName]] else NULL
        cScRanges <- if (!is.null(compounds)) compounds@scoreRanges[[gName]] else NULL
        
        formRank <- NA_integer_; annSimForm <- NA_real_
        if (!is.null(fTable) && !is.null(si[["formula"]]) && !is.na(si$formula[i]))
        {
            formRank <- which(si$formula[i] == fTable$neutral_formula)
            formRank <- if (length(formRank) > 0) formRank[1] else NA_integer_
            if (!is.na(formRank))
                annSimForm <- formulas[[gName]]$annSim[formRank]
        }
        
        suspIK1 <- if (!is.null(si[["InChIKey"]]) && !is.na(si$InChIKey[i])) getIKBlock1(si$InChIKey[i]) else NULL
        annSimComp <- annSimBoth <- NA_real_; compRank <- NA_integer_
        if (!is.null(MSMSList) && !is.null(cTable) && !is.null(suspIK1))
        {
            compRank <- which(suspIK1 == cTable$InChIKey1)
            compRank <- if (length(compRank) > 0) compRank[1] else NA_integer_
            
            if (!is.na(compRank) && !is.null(cTable[["fragInfo"]][[compRank]]))
            {
                annSimComp <- compounds[[gName]]$annSim[compRank]

                if (!is.null(compounds[[gName]][["annSimBoth"]]))
                    annSimBoth <- compounds[[gName]]$annSimBoth[compRank]
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

            if (!is.null(fTable) && "formula" %in% checkFragments && !is.na(formRank) &&
                !is.null(fTable[["fragInfo"]][[formRank]]))
                maxFragMatches <- max(NAToZero(maxFragMatches), sum(fragForms %in% fTable[["fragInfo"]][[formRank]]$ion_formula))

            if (!is.null(cTable) && "compound" %in% checkFragments && !is.na(compRank) &&
                !is.null(cTable[["fragInfo"]][[compRank]]))
                maxFragMatches <- max(NAToZero(maxFragMatches), sum(fragForms %in% cTable[["fragInfo"]][[compRank]]$ion_formula))
        }

        maxFragMatchesRel <- NA_real_
        if (!is.na(maxFragMatches))
            maxFragMatchesRel <- maxFragMatches / maxSuspFrags
        
        fTableNorm <- if (!is.null(fTable))
        {
            normalizeAnnScores(fTable, formScoreNames(TRUE), fScRanges, mergedConsensusNames(formulas),
                               formulasNormalizeScores == "minmax")
        }
        else
            NULL
        cTableNorm <- if (!is.null(cTable))
        {
            normalizeAnnScores(cTable, compScoreNames(TRUE), cScRanges, mergedConsensusNames(compounds),
                               compoundsNormalizeScores == "minmax")
        }
        else
            NULL
        
        estIDLevel <- estimateIdentificationLevel(si$name[i], si$group[i], si$d_rt[i], suspIK1, si$formula[i],
                                                  annSimForm, annSimComp, annSimBoth,
                                                  maxSuspFrags, maxFragMatches, fTable, fTableNorm, formRank,
                                                  mFormNames = if (!is.null(formulas)) mergedConsensusNames(formulas) else character(),
                                                  cTable, cTableNorm, compRank,
                                                  mCompNames = if (!is.null(compounds)) mergedConsensusNames(compounds) else character(),
                                                  absMzDev, IDLevelRules, logPath)
        
        set(si, i,
            c("formRank", "compRank", "annSimForm", "annSimComp", "annSimBoth",
              "maxFrags", "maxFragMatches", "maxFragMatchesRel", "estIDLevel"),
            list(formRank, compRank, annSimForm, annSimComp, annSimBoth, maxSuspFrags, maxFragMatches,
                 maxFragMatchesRel, estIDLevel))
        
        setTxtProgressBar(prog, i)
    }
    
    rmCols <- annCols[sapply(annCols, function(col) !is.null(si[[col]]) && all(is.na(si[[col]])))]
    if (length(rmCols) > 0)
        si <- si[, setdiff(names(si), rmCols), with = FALSE]
    
    obj@screenInfo <- si
    
    close(prog)
    
    saveCacheData("estimateIDConfidence", obj, hash)
    
    return(obj)
})

#' @describeIn featureGroupsScreening Performs rule based filtering. This method builds on the comprehensive filter
#'   functionality from the base \code{\link{filter,featureGroups-method}}. It adds several filters to select
#'   \emph{e.g.} the best ranked suspects or those with a minimum estimated identification level. \strong{NOTE}: most
#'   filters \emph{only} affect suspect hits, not feature groups. Set \code{onlyHits=TRUE} to subsequently remove any
#'   feature groups that lost any suspect matches due to these filter steps.
#'
#' @param selectHitsBy Should be \code{"intensity"} or \code{"level"}. For cases where the same suspect is matched to
#'   multiple feature groups, only the suspect to the feature group with highest mean intensity
#'   (\code{selectHitsBy="intensity"}) or best identification level (\code{selectHitsBy="level"}) is kept. In case of
#'   ties only the first hit is kept. Set to \code{NULL} to ignore this filter. If \code{negate=TRUE} then only those
#'   hits with lowest mean intensity/poorest identification level are kept.
#' @param selectBestFGroups If \code{TRUE} then for any cases where a single feature group is matched to several
#'   suspects only the suspect assigned to the feature group with best identification score is kept. In case of ties
#'   only the first is kept.
#' @param maxLevel,maxFormRank,maxCompRank,minAnnSimForm,minAnnSimComp,minAnnSimBoth Filter suspects by maximum
#'   identification level (\emph{e.g.} \code{"3a"}), formula/compound rank or with minimum formula/compound/combined
#'   annotation similarity. Set to \code{NULL} to ignore.
#' @param absMinFragMatches,relMinFragMatches Only retain suspects with this minimum number MS/MS matches with the
#'   fragments specified in the suspect list (\emph{i.e.} \code{fragments_mz}/\code{fragments_formula}).
#'   \code{relMinFragMatches} sets the minimum that is relative (\samp{0-1}) to the maximum number of MS/MS fragments
#'   specified in the \code{fragments_*} columns of the suspect list. Set to \code{NULL} to ignore.
#' @param minRF Filter suspect hits by the given minimum predicted response factor (as calculated by
#'   \code{\link[=predictRespFactors]{predictRespFactors}}). Set to \code{NULL} to ignore.
#' @param maxLC50 Filter suspect hits by the given maximum toxicity (LC50) (as calculated by
#'   \code{\link[=predictTox]{predictTox}}). Set to \code{NULL} to ignore.
#' @param negate If set to \code{TRUE} then filtering operations are performed in opposite manner.
#' @param applyIMS For IMS workflows: whether the filters are only applied to IMS parents (\code{applyIMS=FALSE}), only
#'   to mobility features (\code{applyIMS=TRUE}) or to both (\code{applyIMS="both"}). Other feature groups will always
#'   be kept. The \code{negate} option does not affect \code{applyIMS}.
#'
#' @return \code{filter} returns a filtered \code{featureGroupsScreening} object.
#'
#' @note \code{filter} removes suspect hits with \code{NA} values when any of the filters related to minimum or maximum
#'   values are applied (unless \code{negate=TRUE}).
#'
#' @export
setMethod("filter", "featureGroupsScreening", function(obj, ..., onlyHits = NULL,
                                                       selectHitsBy = NULL, selectBestFGroups = FALSE,
                                                       maxLevel = NULL, maxFormRank = NULL, maxCompRank = NULL,
                                                       minAnnSimForm = NULL, minAnnSimComp = NULL, minAnnSimBoth = NULL,
                                                       absMinFragMatches = NULL, relMinFragMatches = NULL,
                                                       minRF = NULL, maxLC50 = NULL, negate = FALSE, applyIMS = "both")
{
    # NOTE: keep args and method in sync with featureGroupsScreeningSet method
    
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertFlag, . ~ onlyHits + selectBestFGroups + negate, null.ok = c(TRUE, FALSE, FALSE), fixed = list(add = ac))
    checkmate::assertChoice(selectHitsBy, choices = c("intensity", "level"), null.ok = TRUE, add = ac)
    aapply(checkmate::assertCount, . ~ maxLevel + maxFormRank + maxCompRank + absMinFragMatches + relMinFragMatches,
           null.ok = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ minAnnSimForm + minAnnSimComp + minAnnSimBoth + minRF + maxLC50, null.ok = TRUE,
           fixed = list(add = ac))
    assertApplyIMSArg(applyIMS, add = ac)
    checkmate::reportAssertions(ac)
    
    obj <- doSuspectFilter(obj, onlyHits, selectHitsBy, selectBestFGroups, maxLevel, maxFormRank, maxCompRank,
                           minAnnSimForm, minAnnSimComp, minAnnSimBoth, absMinFragMatches, relMinFragMatches, minRF,
                           maxLC50, negate, applyIMS)

    if (...length() > 0)
        obj <- callNextMethod(obj, ..., negate = negate, applyIMS = applyIMS)
    
    return(obj)
})

#' @param fromSuspects If \code{TRUE} then mobilities values are directly copied from suspect list data (if available).
#'   See the \verb{Post mobility assignment} section.
#'
#' @details In suspect screening workflows \code{assignMobilities} also assigns reference mobility and \acronym{CCS}
#'   values to suspect hits, and can filter hits if \code{IMSMatchParams} is set. This is similarly performed as
#'   \code{\link{screenSuspects}}, please see its documentation for more details.
#'
#' @section Post mobility assignment: \subsection{Suspect screening workflows}{In suspect screening workflows the
#'   \code{fromSuspects} arguments can be set to alternatively perform mobility assignment directly from the suspect
#'   list data (replacing Steps 1-2). The feature mobility is simply assigned from the suspect data and the mobility
#'   range is derived from the \code{IMSWindow} argument. Relationships with IMS parents (Step 2) are similarly formed.
#'   An advantage of this approach is that no mobility peak detection is needed, which may useful for low intensity
#'   features where this could be difficult. However, it is strongly recommended to set \code{fallbackEICs=FALSE} to
#'   still have a means of verification that the mobility feature is actually present.
#'
#'   Setting \code{fromSuspects=TRUE} is primarily intended for workflows where (1) the mobility of a suspect is
#'   accurately known upfront or (2) IMS data should only be used as a rough filtering step for feature data. In the
#'   latter case accurate feature mobility assignment is not of interest and the suspect IMS data is typically not
#'   accurately known (\emph{e.g.} predicted), hence, for these workflows the tolerance specified by \code{IMSWindow}
#'   should be increased.
#'
#'   If both \code{fromSuspects} and \code{mobPeakParams} are set, regular mobility assignment (Steps 1-2) is performed
#'   for features without suspect hit.}
#'
#' @template IMSMatchParams-arg
#'
#' @rdname assignMobilities_feat
#' @aliases assignMobilities,featureGroupsScreening-method
#' @export
setMethod("assignMobilities", "featureGroupsScreening", function(obj, mobPeakParams = NULL, chromPeakParams = NULL,
                                                                 EIMParams = getDefEIMParams(),
                                                                 EICParams = getDefEICParams(),
                                                                 peakRTWindow = defaultLim("retention", "narrow"),
                                                                 fallbackEIC = TRUE, calcArea = "integrate",
                                                                 IMSWindow = defaultLim("mobility", "medium"),
                                                                 CCSParams = NULL, parallel = TRUE,
                                                                 fromSuspects = FALSE, IMSMatchParams = NULL)
{
    ac <- checkmate::makeAssertCollection()
    assertFindMobilitiesArgs(mobPeakParams, chromPeakParams, EIMParams, EICParams, peakRTWindow, fallbackEIC,
                             calcArea, IMSWindow, CCSParams, parallel, ac)
    checkmate::assertFlag(fromSuspects, add = ac)
    assertIMSMatchParams(IMSMatchParams, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    if (length(obj) == 0 || (is.null(mobPeakParams) && !fromSuspects && is.null(CCSParams)))
        return(obj) # nothing to do...

    if (!is.null(mobPeakParams) || fromSuspects)
    {
        obj <- checkAssignedMobilityFGroups(obj)
        if (fromSuspects)
            obj@features <- assignFeatureMobilitiesSuspects(obj@features, IMSWindow, screenInfo(obj))
        if (!is.null(mobPeakParams))
            obj@features <- assignFeatureMobilitiesPeaks(obj@features, mobPeakParams, EIMParams)
        obj@features <- reintegrateMobilityFeatures(obj@features, chromPeakParams, EICParams, peakRTWindow, fallbackEIC,
                                                    calcArea, parallel)
        obj <- clusterFGroupMobilities(obj, IMSWindow, FALSE)
    }

    if (!is.null(CCSParams))
        obj <- assignFGroupsCCS(obj, CCSParams)
    
    gInfo <- groupInfo(obj)
    scr <- expandAndUpdateScreenInfoForIMS(screenInfo(obj), gInfo)
    scr <- finalizeScreenInfoForIMS(scr, gInfo, IMSMatchParams)
    obj@screenInfo <- scr
    
    return(obj)
})

#' Target and suspect screening
#'
#' Utilities to screen for analytes with known or suspected identity.
#'
#' Besides 'full non-target analysis', where compounds may be identified with little to no prior knowledge, a common
#' strategy is to screen for compounds with known or suspected identity. This may be a generally favorable approach if
#' possible, as it can significantly reduce the load on data interpretation.
#'
#' \code{screenSuspects} is used to perform suspect screening. The input \code{\link{featureGroups}} object will be
#' screened for suspects by \emph{m/z} values and optionally retention times. Afterwards, any feature groups not matched
#' may be kept or removed, depending whether a full non-target analysis is desired.
#'
#' @param fGroups The \code{\link{featureGroups}} object that should be screened.
#' @param suspects A \code{data.frame} with suspect information. See the \verb{Suspect list format} section below.
#'
#'   \setsWF Can also be a \code{list} with suspect lists to be used for each set (otherwise the same suspect lists is
#'   used for all sets). The \code{list} can be named with the sets names to mark which suspect list is to be used with
#'   which set (\emph{e.g.} \code{suspects=list(positive=suspsPos, negative=suspsNeg)}).
#' @param rtWindow,mzWindow The retention time window (in seconds) and \emph{m/z} window that will be used for matching
#'   a suspect (+/- feature data).
#' @param adduct An \code{\link{adduct}} object (or something that can be converted to it with \code{\link{as.adduct}}).
#'   Examples: \code{"[M-H]-"}, \code{"[M+Na]+"}. May be \code{NULL}, see \verb{Suspect list format} and \verb{Matching
#'   of suspect masses} sections below.
#' @param skipInvalid If set to \code{TRUE} then suspects with invalid data (\emph{e.g.} missing names or other missing
#'   data) will be ignored with a warning. Similarly, any suspects for which mass calculation failed (when no \code{mz}
#'   column is present in the suspect list), for instance, due to invalid \code{SMILES}, will be ignored with a warning.
#' @param onlyHits If \code{TRUE} then all feature groups not matched by any of the suspects will be removed.
#'
#' @template IMSMatchParams-arg
#'
#' @section Suspect list format: the \code{suspects} argument for \code{screenSuspects} should be a \code{data.frame}
#'   with the following mandatory and optional columns:
#'
#'   \itemize{
#'
#'   \item \code{name} The suspect name. Must be file-compatible. (\strong{mandatory})
#'
#'   \item \code{rt} The retention time (in seconds) for the suspect. If specified the suspect will only be matched if
#'   its retention matches the experimental value (tolerance defined by the \code{rtWindow} argument).
#'   (\strong{optional})
#'
#'   \item \code{neutralMass},\code{formula},\code{SMILES},\code{InChI} The neutral monoisotopic mass, chemical formula,
#'   SMILES or InChI for the suspect. (data from one of these columns are \strong{mandatory} in case no value from the
#'   \code{mz} column is available for a suspect)
#'
#'   \item \code{mz} The ionized \emph{m/z} of the suspect. (\strong{mandatory} unless it can be calculated from one of
#'   the aforementioned columns)
#'
#'   \item \code{adduct} A \code{character} that can be converted with \code{\link{as.adduct}}. Can be used to
#'   automatically calculate values for the \code{mz} column. (\strong{mandatory} unless data from the \code{mz} column
#'   is available, the \code{adduct} argument is set or \code{fGroups} has adduct annotations)
#'
#'   \item \code{fragments_mz},\code{fragments_formula} One or more MS/MS fragments (specified as \emph{m/z} or
#'   formulae, respectively). Multiple values can be specified by separating them with a semicolon (\verb{;}). This data
#'   is used by \code{\link{estimateIDConfidence}} to report detected MS/MS fragments and calculate identification levels.
#'   (\strong{optional})
#'
#'   \item \code{mobility},\code{CCS} The mobility or \acronym{CCS} value of the suspect. These values may be used to
#'   filter out suspects, see the \code{IMSMatchParams} argument. Multiple values for a single suspect can be specified
#'   by separating them with a semicolon(\verb{;}). Adduct specific columns may be added by suffixing the adduct to the
#'   column name, \emph{e.g.} \code{mobility_[M+H]+} and \code{CCS_[M-H]-}. (\strong{optional})
#'
#'   }
#'
#' @section Matching of suspect masses: How the mass of a suspect is matched with the mass of a feature depends on the
#'   available data: \itemize{
#'
#'   \item If the suspect has data from the \code{mz} column of the suspect list, then this data is matched with the
#'   detected feature \emph{m/z}.
#'
#'   \item Otherwise, if the suspect has data in the \code{adduct} column of the suspect list, this data is used to
#'   calculate its \code{mz} value, which is then used like above.
#'
#'   \item In the last case, the neutral mass of the suspect is matched with the neutral mass of the feature. Hence,
#'   either the \code{adduct} argument needs to be specified, or the \code{featureGroups} input object must have adduct
#'   annotations.
#'
#'   }
#'
#' @section IMS reference assignment: If both adduct specific and non-adduct specific reference values are available,
#'   then non-adduct specific data is chosen (unless \code{NA}) as reference for the suspect hit. Otherwise, data is
#'   taken from the adduct specific data corresponding to the adduct assigned to the feature group (or \code{adduct}
#'   argument). If multiple mobility or \acronym{CCS} values for a suspect are specified in the suspect list, then the
#'   reference value is chosen which is the closest to that of the feature.
#'
#' @templateVar whatCP suspect list
#' @template chemPropCalc
#'
#' @return \code{screenSuspects} returns a \code{\link{featureGroupsScreening}} object, which is a copy of the input
#'   \code{fGroups} object amended with additional screening information.
#'
#' @note \code{screenSuspects} may use the suspect names to base file names used for reporting, logging etc. Therefore,
#'   it is important that these are file-compatible names. For this purpose, \code{screenSuspects} will automatically
#'   try to convert long, non-unique and/or otherwise incompatible suspect names.
#'
#' @seealso \code{featureGroupsScreening}
#'
#' @name suspect-screening
#' @aliases screenSuspects
#' @aliases screenSuspects,featureGroups-method
#' @export
setMethod("screenSuspects", "featureGroups", function(fGroups, suspects, rtWindow, mzWindow, IMSMatchParams,
                                                      adduct, skipInvalid, prefCalcChemProps, neutralChemProps,
                                                      onlyHits)
{
    checkmate::assertFlag(skipInvalid) # not in assert collection, should fail before assertSuspectList

    # NOTE: skip adduct check if fGroups is empty    
    needsAdduct <- is.null(adduct) && nrow(annotations(fGroups)) == 0 && length(fGroups) > 0
    
    ac <- checkmate::makeAssertCollection()
    assertSuspectList(suspects, needsAdduct = needsAdduct, skipInvalid, add = ac)
    aapply(checkmate::assertNumber, . ~ rtWindow + mzWindow, lower = 0, finite = TRUE,
           fixed = list(add = ac))
    assertIMSMatchParams(IMSMatchParams, null.ok = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ skipInvalid + prefCalcChemProps + neutralChemProps + onlyHits,
           fixed = list(add = ac))
    checkmate::reportAssertions(ac)

    if (!is.null(adduct))
        adduct <- checkAndToAdduct(adduct, fGroups)
    
    # do this before checking cache to ensure proper errors/warnings are thrown!
    suspects <- prepareSuspectList(suspects, adduct, skipInvalid, checkDesc = TRUE,
                                   prefCalcChemProps = prefCalcChemProps, neutralChemProps = neutralChemProps)
    
    hash <- makeHash(fGroups, suspects, rtWindow, mzWindow, IMSMatchParams, adduct, skipInvalid, prefCalcChemProps,
                     neutralChemProps, onlyHits)
    cd <- loadCacheData("screenSuspects", hash)
    if (!is.null(cd))
        return(cd)

    scr <- doScreenSuspects(fGroups, suspects, rtWindow, mzWindow, IMSMatchParams, adduct, skipInvalid)

    if (onlyHits)
        fGroups <- fGroups[, scr$group]
    
    ret <- featureGroupsScreening(screenInfo = scr, groups = copy(groupTable(fGroups)),
                                  groupInfo = copy(groupInfo(fGroups)), features = getFeatures(fGroups),
                                  ftindex = copy(groupFeatIndex(fGroups)),
                                  groupQualities = copy(groupQualities(fGroups)),
                                  groupScores = copy(groupScores(fGroups)), ISTDs = copy(internalStandards(fGroups)),
                                  ISTDAssignments = internalStandardAssignments(fGroups),
                                  annotations = copy(annotations(fGroups)),
                                  concentrations = copy(concentrations(fGroups)),
                                  toxicities = copy(toxicities(fGroups)))
    
    saveCacheData("screenSuspects", ret, hash)
    
    return(ret)
})

#' @param amend If \code{TRUE} then screening results will be \emph{amended} to the original object.
#' @rdname suspect-screening
#' @export
setMethod("screenSuspects", "featureGroupsScreening", doScreenSuspectsAmend)
