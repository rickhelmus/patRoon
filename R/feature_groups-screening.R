#' @include main.R
#' @include feature_groups.R
#' @include feature_groups-set.R
#' @include utils-screening.R
NULL

screeningSlots <- c(screenInfo = "data.table")

#' Class for suspect screened feature groups.
#'
#' This class derives from \code{\link{featureGroups}} and adds suspect screening information.
#'
#' @param obj,object,x,fGroups The \code{featureGroupsScreening} object.
#' @param \dots Further arguments passed to the base method.
#' @param onlyHits For \code{as.data.table}: if \code{TRUE} then only feature groups with suspect hits are reported.
#'
#'   For \code{filter} \itemize{
#'
#'   \item if \code{negate=FALSE} and \code{onlyHits=TRUE} then all feature groups without suspect hits will be removed.
#'   Otherwise nothing will be done.
#'
#'   \item if \code{negate=TRUE} then \code{onlyHits=TRUE} will select feature groups without suspect hits,
#'   \code{onlyHits=FALSE} will only retain feature groups with suspect matches and this filter is ignored if
#'   \code{onlyHits=NULL}.
#'
#'   }
#'
#' @slot screenInfo A (\code{\link{data.table}}) with results from suspect screening. This table will be amended with
#'   annotation data when \code{annotateSuspects} is run.
#'
#' @section Suspect annotation: The \code{annotateSuspects} method is used to annotate suspects after
#'   \code{\link{screenSuspects}} was used to collect suspect screening results and other workflow steps such as formula
#'   and compound annotation steps have been completed. The annotation results, which can be acquired with the
#'   \code{as.data.table} and \code{screenInfo} methods, amends the current screening data with the following columns:
#'
#'   \itemize{
#'
#'   \item \code{formRank},\code{compRank} The rank of the suspect within the formula/compound annotation results.
#'
#'   \item \code{annSimForm},\code{annSimComp},\code{annSimBoth} A similarity measure between measured and annotated
#'   MS/MS peaks from annotation of formulae, compounds or both. The similarity is calculated as the spectral similarity
#'   between a peaklist with (a) all MS/MS peaks and (b) only annotated peaks. Thus, a value of one means that all MS/MS
#'   peaks were annotated. If both formula and compound annotations are available then \code{annSimBoth} is calculated
#'   after combining all the annotated peaks, otherwise \code{annSimBoth} equals the available value for
#'   \code{annSimForm} or \code{annSimComp}. The similarity calculation can be configured with the
#'   \code{relMinMSMSIntensity} and \code{simMSMSMethod} arguments to \code{annotateSuspects}.
#'
#'   \item \code{maxFrags} The maximum number of MS/MS fragments that can be matched for this suspect (based on the
#'   \code{fragments_*} columns from the suspect list).
#'
#'   \item \code{maxFragMatches},\code{maxFragMatchesRel} The absolute and relative amount of experimental MS/MS peaks
#'   that were matched from the fragments specified in the suspect list. The value for \code{maxFragMatchesRel} is
#'   relative to the value for \code{maxFrags}. The calculation of this column is influenced by the
#'   \code{checkFragments} argument to \code{annotateSuspects}.
#'
#'   \item \code{estIDLevel} Provides an \emph{estimation} of the identification level, roughly following that of
#'   \insertCite{Schymanski2014}{patRoon}. However, please note that this value is only an estimation, and manual
#'   interpretation is still necessary to assign final identification levels. The estimation is done through a set of
#'   rules, see the \verb{Identification level rules} section below.
#'
#'   }
#'
#'   Note that only columns are present is sufficient data is available for their calculation.
#'
#' @section Identification level rules: The estimation of identification levels is configured through a YAML file which
#'   specifies the rules for each level. The default file is shown below.
#'
#' @eval paste0("@@section Identification level rules: \\preformatted{", patRoon:::readAllFile(system.file("misc",
#'   "IDLevelRules.yml", package = "patRoon")), "}")
#'
#' @section Identification level rules: Most of the file should be self-explanatory. Some notes:
#'
#'   \itemize{
#'
#'   \item Each rule is either a field of \code{suspectFragments} (minimum number of MS/MS fragments matched from
#'   suspect list), \code{retention} (maximum retention deviation from suspect list), \code{rank} (the maximum
#'   annotation rank from formula or compound annotations), \code{all} (this level is always matched) or any of the
#'   scorings available from the formula or compound annotations.
#'
#'   \item In case any of the rules could be applied to either formula or compound annotations, the annotation type must
#'   be specified with the \code{type} field (\code{formula} or \code{compound}).
#'
#'   \item Identification levels should start with a number and may optionally be followed by a alphabetic character.
#'   The lowest levels are checked first.
#'
#'   \item If \code{relative=yes} then the relative scoring will be used for testing.
#'
#'   \item For \code{suspectFragments}: if the number of fragments from the suspect list (\code{maxFrags} column) is
#'   less then the minimum rule value, the minimum is adjusted to the number of available fragments.
#'
#'   }
#'
#'   A template rules file can be generated with the \code{\link{genIDLevelRulesFile}} function, and this file can
#'   subsequently passed to \code{annotateSuspects}. The file format is highly flexible and (sub)levels can be added or
#'   removed if desired. Note that the default file is currently only suitable when annotation is performed with GenForm
#'   and MetFrag, for other algorithms it is crucial to modify the rules.
#'
#' @templateVar class featureGroupsScreening
#' @template class-hierarchy
#'
#' @references \insertAllCited{} \cr \cr \insertRef{Stein1994}{patRoon}
#'
#' @seealso \code{\link{featureGroups}}
#'
#' @export
featureGroupsScreening <- setClass("featureGroupsScreening",
                                   slots = c(screenInfo = "data.table"),
                                   contains = "featureGroups")

setMethod("initialize", "featureGroupsScreening",
          function(.Object, ...) callNextMethod(.Object, algorithm = "screening", ...))

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
#' @param i,j,rGroups Used for subsetting data analyses, feature groups and
#'   replicate groups, see \code{\link{featureGroups}}.
#' @param suspects An optional \code{character} vector with suspect names. If
#'   specified, only \code{featureGroups} will be kept that are assigned to
#'   these suspects.
#' @param drop Ignored.
#'
#' @export
setMethod("[", c("featureGroupsScreening", "ANY", "ANY", "missing"), function(x, i, j, ..., rGroups,
                                                                              suspects = NULL, drop = TRUE)
{
    checkmate::assertCharacter(suspects, null.ok = TRUE)
    
    x <- callNextMethod(x, i, j, ..., rGroups = rGroups, drop = drop)
    
    if (!is.null(suspects))
        x <- x[, x@screenInfo[name %in% suspects]$group]
    
    return(x)
})

#' @rdname featureGroupsScreening-class
#' @export
setMethod("delete", "featureGroupsScreening", function(obj, i = NULL, j = NULL, ...)
{
    obj <- callNextMethod()
    obj@screenInfo <- obj@screenInfo[group %in% names(obj)]
    return(obj)
})

#' @describeIn featureGroupsScreening Obtain a summary table (a
#'   \code{\link{data.table}}) with retention, \emph{m/z}, intensity and
#'   optionally other feature data. Furthermore, the output table will be merged
#'   with information from \code{screenInfo}, such as suspect names and other
#'   properties and annotation data.
#'
#' @param collapseSuspects If a \code{character} then any suspects that were
#'   matched to the same feature group are collapsed to a single row and suspect
#'   names are separated by the value of \code{collapseSuspects}. If \code{NULL}
#'   then no collapsing occurs, and each suspect match is reported on a single
#'   row. Note that some columns will not be reported when collapsing is
#'   enabled.
#'
#' @export
setMethod("as.data.table", "featureGroupsScreening", function(x, ..., collapseSuspects = ",", onlyHits = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertString(collapseSuspects, null.ok = TRUE, add = ac)
    checkmate::assertFlag(onlyHits, add = ac)
    checkmate::reportAssertions(ac)
    
    ret <- callNextMethod(x, ...)
    if (nrow(ret) > 0)
        ret <- mergeScreenInfoWithDT(ret, screenInfo(x), collapseSuspects, onlyHits)

    return(ret)
})

#' @describeIn featureGroupsScreening Incorporates annotation data obtained during the workflow to annotate suspects
#'   with matched known MS/MS fragments, formula/candidate ranks and automatic estimation of identification levels. See
#'   the \verb{Suspect annotation} section for more details. The estimation of identification levels for each suspect is
#'   logged in the \code{log/ident} directory.
#'
#' @templateVar normParam compoundsNormalizeScores,formulasNormalizeScores
#' @templateVar noNone TRUE
#' @template norm-args
#'
#' @param MSPeakLists,formulas,compounds Annotation data (\code{\link{MSPeakLists}}, \code{\link{formulas}} and
#'   \code{\link{compounds}}) obtained for this \code{featureGroupsScreening} object. All arguments can be \code{NULL}
#'   to exclude it from the annotation.
#' @param absMzDev Maximum absolute \emph{m/z} deviation.
#' @param relMinMSMSIntensity Minimum relative intensity (\samp{0-1}) threshold applied when calculating annotation
#'   similarities.
#' @param simMSMSMethod Either \code{"cosine"} or \code{"jaccard"}: used to compare MS/MS peak lists for annotation
#'   similarity calculation.
#' @param checkFragments Which type(s) of MS/MS fragments from workflow data should be checked to evaluate the number of
#'   suspect fragment matches (\emph{i.e.} from the \code{fragments_mz}/\code{fragments_formula} columns in the suspect
#'   list). Valid values are: \code{"mz"}, \code{"formula"}, \code{"compounds"}. The former uses \emph{m/z} values in
#'   the specified \code{MSPeakLists} object, whereas the others use the formulae that were annotated to MS/MS peaks in
#'   the given \code{formulas} or \code{compounds} objects. Multiple values are possible: in this case the maximum
#'   number of fragment matches will be reported.
#' @param IDFile A file path to a YAML file with rules used for estimation of identification levels. See the
#'   \verb{Suspect annotation} section for more details. If not specified then a default rules file will be used.
#' @param logPath A directory path to store logging information. If \code{NULL} then logging is disabled.
#'
#' @return \code{annotateSuspects} returns a \code{featureGroupsScreening} object, which is a
#'   \code{\link{featureGroups}} object amended with annotation data.
#'
#' @note The \code{relMinMSMSIntensity} filter argument to \code{annotateSuspects} is applied \emph{after} removing the
#'   precursor ion from the peak lists (if present). Thus, intensity scales may be different when this filter is applied
#'   when the most abundant peak resulted from the precursor ion.
#'
#' @author Rick Helmus <\email{r.helmus@@uva.nl}>, Emma Schymanski <\email{emma.schymanski@@uni.lu}> (contributions to
#'   identification level rules), Bas van de Velde (contributions to spectral similarity calculation).
#'
#' @section Source: Cosine spectral similarity calculation was based on the code from \code{SpectrumSimilarity()}
#'   function of \pkg{OrgMassSpecR}.
#'
#' @aliases annotateSuspects
#' @export
setMethod("annotateSuspects", "featureGroupsScreening", function(fGroups, MSPeakLists, formulas, compounds,
                                                                 absMzDev = 0.005, relMinMSMSIntensity = 0.05,
                                                                 simMSMSMethod = "cosine",
                                                                 checkFragments = c("mz", "formula", "compound"),
                                                                 formulasNormalizeScores = "max",
                                                                 compoundsNormalizeScores = "max",
                                                                 IDFile = system.file("misc", "IDLevelRules.yml",
                                                                                      package = "patRoon"),
                                                                 logPath = file.path("log", "ident"))
{
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
    if (!is.null(logPath))
        assertCanCreateDir(logPath, add = ac)
    checkmate::reportAssertions(ac)

    hash <- makeHash(fGroups, MSPeakLists, formulas, compounds, absMzDev,
                     relMinMSMSIntensity, simMSMSMethod, checkFragments, formulasNormalizeScores,
                     compoundsNormalizeScores, makeFileHash(IDFile))
    cd <- loadCacheData("annotateSuspects", hash)
    if (!is.null(cd))
        return(cd)
    
    IDLevelRules <- readYAML(IDFile)
    
    if (!checkmate::test_named(IDLevelRules))
        stop("No valid rules could be loaded")
    if (!all(grepl("^[[:digit:]]+[[:alpha:]]?$", names(IDLevelRules))))
        stop("Levels should be defined as a number and may optionally followed by one character (e.g. 3, 2b etc)")
    
    IDLevelRules <- IDLevelRules[order(names(IDLevelRules))] # sort to ensure lowest levels will be tested first

    if (nrow(screenInfo(fGroups)) == 0)
    {
        cat("No suspect hits, nothing to annotate")
        return(fGroups)
    }
        
        
    mzWithin <- function(mz1, mz2) abs(mz1 - mz2) <= absMzDev

    si <- copy(screenInfo(fGroups))
    annCols <- c("formRank", "compRank", "annSimForm", "annSimComp", "annSimBoth", "maxFrags", "maxFragMatches",
                 "maxFragMatchesRel", "estIDLevel")
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
        
        formRank <- NA_integer_; annSimForm <- annSimBoth <- NA_real_
        if (!is.null(fTable) && !is.null(si[["formula"]]) && !is.na(si$formula[i]))
        {
            formRank <- which(si$formula[i] == fTable$neutral_formula)
            formRank <- if (length(formRank) > 0) formRank[1] else NA_integer_
            if (!is.na(formRank))
                annSimForm <- annSimBoth <- annotatedMSMSSimilarity(annotatedPeakList(formulas,
                                                                                      index = formRank,
                                                                                      groupName = gName,
                                                                                      MSPeakLists = MSPeakLists),
                                                                    absMzDev, relMinMSMSIntensity, simMSMSMethod)
        }
        
        suspIK1 <- if (!is.null(si[["InChIKey"]]) && !is.na(si$InChIKey[i])) getIKBlock1(si$InChIKey[i]) else NULL
        annSimComp <- NA_real_; compRank <- NA_integer_
        if (!is.null(MSMSList) && !is.null(cTable) && !is.null(suspIK1))
        {
            compRank <- which(suspIK1 == cTable$InChIKey1)
            compRank <- if (length(compRank) > 0) compRank[1] else NA_integer_
            
            if (!is.na(compRank) && !is.null(cTable[["fragInfo"]][[compRank]]))
            {
                annSimComp <- annotatedMSMSSimilarity(annotatedPeakList(compounds, index = compRank,
                                                                        groupName = gName, MSPeakLists = MSPeakLists),
                                                      absMzDev, relMinMSMSIntensity, simMSMSMethod)
                
                if (!is.na(formRank))
                    annSimBoth <- annotatedMSMSSimilarity(annotatedPeakList(compounds, index = compRank,
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
        
        estIDLevel <- estimateIdentificationLevel(si$name[i], si$group[i], si$d_rt[i], suspIK1, si$formula[i],
                                                  annSimForm, annSimComp, annSimBoth,
                                                  maxSuspFrags, maxFragMatches, fTable, formRank,
                                                  mFormNames = if (!is.null(formulas)) mergedConsensusNames(formulas) else character(),
                                                  fScRanges, formulasNormalizeScores, cTable, compRank,
                                                  mCompNames = if (!is.null(compounds)) mergedConsensusNames(compounds) else character(),
                                                  cScRanges, compoundsNormalizeScores, absMzDev, IDLevelRules, logPath)
        
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
    
    fGroups@screenInfo <- si
    
    close(prog)
    
    saveCacheData("annotateSuspects", fGroups, hash)
    
    return(fGroups)
})

#' @describeIn featureGroupsScreening Performs rule based filtering. This method
#'   builds on the comprehensive filter functionality from the base
#'   \code{\link{filter,featureGroups-method}}. It adds several filters to
#'   select \emph{e.g.} the best ranked suspects or those with a minimum
#'   estimated identification level. \strong{NOTE}: most filters \emph{only}
#'   affect suspect hits, not feature groups. Set \code{onlyHits=TRUE} to
#'   subsequently remove any feature groups that lost any suspect matches due to
#'   other filter steps.
#'
#' @param selectHitsBy Should be \code{"intensity"} or \code{"level"}. For cases
#'   where the same suspect is matched to multiple feature groups, only the
#'   suspect to the feature group with highest mean intensity
#'   (\code{selectHitsBy="intensity"}) or best identification level
#'   (\code{selectHitsBy="level"}) is kept. In case of ties only the first hit
#'   is kept. Set to \code{NULL} to ignore this filter. If \code{negate=TRUE}
#'   then only those hits with lowest mean intensity/poorest identification
#'   level are kept.
#' @param selectBestFGroups If \code{TRUE} then for any cases where a single
#'   feature group is matched to several suspects only the suspect assigned to
#'   the feature group with best identification score is kept. In case of ties
#'   only the first is kept.
#' @param
#' maxLevel,maxFormRank,maxCompRank,minAnnSimForm,minAnnSimComp,minAnnSimBoth
#' Filter suspects by maximum identification level (\emph{e.g.} \code{"3a"}),
#' formula/compound rank or with minimum formula/compound/combined annotation
#' similarity. Set to \code{NULL} to ignore.
#' @param absMinFragMatches,relMinFragMatches Only retain suspects with this
#'   minimum number MS/MS matches with the fragments specified in the suspect
#'   list (\emph{i.e.} \code{fragments_mz}/\code{fragments_formula}).
#'   \code{relMinFragMatches} sets the minimum that is relative (\samp{0-1}) to
#'   the maximum number of MS/MS fragments specified in the \code{fragments_*}
#'   columns of the suspect list. Set to \code{NULL} to ignore.
#' @param negate If set to \code{TRUE} then filtering operations are performed
#'   in opposite manner.
#'
#' @return \code{filter} returns a filtered \code{featureGroupsScreening}
#'   object.
#'
#' @note \code{filter} removes suspect hits with \code{NA} values when any of
#'   the filters related to minimum or maximum values are applied (unless
#'   \code{negate=TRUE}).
#'
#' @export
setMethod("filter", "featureGroupsScreening", function(obj, ..., onlyHits = NULL,
                                                       selectHitsBy = NULL, selectBestFGroups = FALSE,
                                                       maxLevel = NULL, maxFormRank = NULL, maxCompRank = NULL,
                                                       minAnnSimForm = NULL, minAnnSimComp = NULL, minAnnSimBoth = NULL,
                                                       absMinFragMatches = NULL, relMinFragMatches = NULL, negate = FALSE)
{
    # NOTE: keep args and method in sync with featureGroupsScreeningSet method
    
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertFlag, . ~ onlyHits + selectBestFGroups + negate, null.ok = c(TRUE, FALSE, FALSE), fixed = list(add = ac))
    checkmate::assertChoice(selectHitsBy, choices = c("intensity", "level"), null.ok = TRUE, add = ac)
    aapply(checkmate::assertCount, . ~ maxLevel + maxFormRank + maxCompRank + absMinFragMatches + relMinFragMatches,
           null.ok = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ minAnnSimForm + minAnnSimComp + minAnnSimBoth, null.ok = TRUE, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    obj <- doSuspectFilter(obj, onlyHits, selectHitsBy, selectBestFGroups, maxLevel, maxFormRank, maxCompRank,
                           minAnnSimForm, minAnnSimComp, minAnnSimBoth, absMinFragMatches, relMinFragMatches, negate)

    if (...length() > 0)
        obj <- callNextMethod(obj, ..., negate = negate)
    
    return(obj)
})


#' @details \code{screenSuspects} is used to perform suspect screening. The input \code{\link{featureGroups}} object
#'   will be screened for suspects by \emph{m/z} values and optionally retention times. Afterwards, any feature groups
#'   not matched may be kept or removed, depending whether a full non-target analysis is desired.
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
#'   is used by \code{\link{annotateSuspects}} to report detected MS/MS fragments and calculate identification levels.
#'   (\strong{optional})
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
#'
#' @return \code{screenSuspects} returns a \code{\link{featureGroupsScreening}} object, which is a copy of the input
#'   \code{fGroups} object amended with additional screening information.
#'
#' @note For \code{screenSuspects} in some cases you may need to install
#'   \href{http://openbabel.org/wiki/Main_Page}{OpenBabel} (\emph{e.g.} when only InChI data is available for mass
#'   calculation).
#'
#' @seealso \code{featureGroupsScreening}
#'
#' @references \insertRef{OBoyle2011}{patRoon}
#'
#' @rdname suspect-screening
#' @aliases screenSuspects
#' @export
setMethod("screenSuspects", "featureGroups", function(fGroups, suspects, rtWindow, mzWindow,
                                                      adduct, skipInvalid, onlyHits)
{
    checkmate::assertFlag(skipInvalid) # not in assert collection, should fail before assertSuspectList

    # NOTE: skip adduct check if fGroups is empty    
    needsAdduct <- is.null(adduct) && nrow(annotations(fGroups)) == 0 && length(fGroups) > 0
    
    ac <- checkmate::makeAssertCollection()
    assertSuspectList(suspects, needsAdduct = needsAdduct, skipInvalid, add = ac)
    aapply(checkmate::assertNumber, . ~ rtWindow + mzWindow, lower = 0, finite = TRUE, fixed = list(add = ac))
    checkmate::assertFlag(onlyHits, add = ac)
    checkmate::reportAssertions(ac)

    if (!is.null(adduct))
        adduct <- checkAndToAdduct(adduct, fGroups)
    
    # do this before checking cache to ensure proper errors/warnings are thrown!
    suspects <- prepareSuspectList(suspects, adduct, skipInvalid)
    
    hash <- makeHash(fGroups, suspects, rtWindow, mzWindow, adduct, skipInvalid, onlyHits)
    cd <- loadCacheData("screenSuspects", hash)
    if (!is.null(cd))
        return(cd)

    scr <- doScreenSuspects(fGroups, suspects, rtWindow, mzWindow, skipInvalid)

    if (onlyHits)
        fGroups <- fGroups[, scr$group]
    
    ret <- featureGroupsScreening(screenInfo = scr, groups = copy(groupTable(fGroups)),
                                  analysisInfo = analysisInfo(fGroups), groupInfo = groupInfo(fGroups),
                                  features = getFeatures(fGroups), ftindex = copy(groupFeatIndex(fGroups)),
                                  groupQualities = copy(groupQualities(fGroups)),
                                  groupScores = copy(groupScores(fGroups)), iSTDs = copy(internalStandards(fGroups)),
                                  annotations = copy(annotations(fGroups)))
    
    saveCacheData("screenSuspects", ret, hash)
    
    return(ret)
})

#' @param amend If \code{TRUE} then screening results will be \emph{amended} to the original object.
#' @rdname suspect-screening
#' @export
setMethod("screenSuspects", "featureGroupsScreening", function(fGroups, suspects, rtWindow, mzWindow,
                                                               adduct, skipInvalid, onlyHits, amend = FALSE)
{
    aapply(checkmate::assertFlag, . ~ onlyHits + amend)
    
    fGroupsScreened <- callNextMethod(fGroups, suspects, rtWindow, mzWindow, adduct, skipInvalid, onlyHits)
    if (!amend)
        return(fGroupsScreened)
    
    # amend screening results
    
    fGroups@screenInfo <- rbind(fGroups@screenInfo, fGroupsScreened@screenInfo, fill = TRUE)
    fGroups@screenInfo <- unique(fGroups@screenInfo, by = c("name", "group"))
    
    if (onlyHits)
        fGroups <- fGroups[, fGroups@screenInfo$group]
    
    return(fGroups)
})

#' @details \code{numericIDLevel} Extracts the numeric part of a given
#'   identification level (\emph{e.g.} \code{"3a"} becomes \samp{3}).
#' @param level The identification level to be converted.
#' @rdname suspect-screening
#' @export
numericIDLevel <- function(level)
{
    checkmate::assertCharacter(level, any.missing = TRUE, min.chars = 1)
    ret <- integer(length(level))
    ret[is.na(level)] <- NA_integer_
    ret[!is.na(level)] <- as.integer(gsub("[[:alpha:]]*", "", level[!is.na(level)]))
    return(ret)
}

#' @details \code{genIDLevelRulesFile} Generates a template YAML file that is
#'   used to configure the rules for automatic estimation of identification
#'   levels. This file can then be used as input for
#'   \code{\link{annotateSuspects}}.
#' @param out The file path to the target file.
#' @param inLevels,exLevels A \link[=regex]{regular expression} for the
#'   identification levels to include or exclude, respectively. For instance,
#'   \code{exLevels="4|5"} would exclude level 4 and 5 from the output file. Set
#'   to \code{NULL} to ignore.
#' @rdname suspect-screening
#' @export
genIDLevelRulesFile <- function(out, inLevels = NULL, exLevels = NULL)
{
    aapply(checkmate::assertCharacter, . ~ inLevels + exLevels, null.ok = TRUE)
    checkmate::assertPathForOutput(basename(out), overwrite = TRUE)
    
    defFile <- system.file("misc", "IDLevelRules.yml", package = "patRoon")
    
    if (is.null(inLevels) && is.null(exLevels))
        file.copy(defFile, out, overwrite = TRUE)
    else
    {
        rules <- readYAML(defFile)
        if (!is.null(inLevels))
            rules <- rules[grepl(inLevels, names(rules))]
        if (!is.null(exLevels))
            rules <- rules[!grepl(exLevels, names(rules))]
        # UNDONE: this quotes ID levels without sub-level, fix?
        yaml::write_yaml(rules, out, indent = 4)
    }
    invisible(NULL)
}
