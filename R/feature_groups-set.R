#' @include main.R
#' @include feature_groups.R
#' @include features-set.R
NULL

#' @param set \setsWF The name of the set.
#' @param sets \setsWF For \code{[}: a \code{character} with name(s) of the sets to keep.
#'
#'   For \code{overlap} and \code{unique}: If \code{TRUE} then the \code{which} argument changes its meaning and is used
#'   to specify the names of the sets to be compared.
#'
#' @slot groupAlgo,groupArgs,groupVerbose \setsWF Grouping parameters that were used when this object was created. Used
#'   by \code{adducts<-} and \code{selectIons} when these methods perform a re-grouping of features.
#' @slot annotations,ISTDAssignments \setsWF As the \code{featureGroups} slots, but contains the data per set.
#' @slot annotationsChanged Set internally by \code{adducts()<-} and applied as soon as \code{reGroup=TRUE}. 
#'
#' @section Sets workflows: \setsWFClass{featureGroupsSet}{featureGroups}
#'
#'   \setsWFNewMethodsFeat{featureGroupsUnset}{The adduct annotations for the selected set are used to convert all
#'   feature (group) masses to ionic \emph{m/z} values. The annotations persist in the converted object. }
#'
#'   \setsWFChangedMethods{
#'
#'   \item \code{adducts}, \code{adducts<-} require the \code{set} argument. The order of the data that is
#'   returned/changed follows that of the \code{annotations} slot. Furthermore, \code{adducts<-} will perform a
#'   re-grouping of features when its \code{reGroup} parameter is set to \code{TRUE}. The implications for this are
#'   discussed below. Note that no adducts are changed \emph{until} \code{reGroup=TRUE}.
#'
#'   \item the subset operator (\code{[}) has specific arguments to choose (feature presence in) sets. See the argument
#'   descriptions.
#'
#'   \item \code{as.data.table}: normalization of intensities is performed per set.
#'
#'   \item \code{export} Only allows to export data from one set. The \code{unset} method is used prior to exporting the
#'   data.
#'
#'   \item \code{overlap} and \code{unique} allow to handle data per set. See the \code{sets} argument description.
#'
#'   \item \code{selectIons} Will perform a re-grouping of features. The implications of this are discussed below.
#'   
#'   \item \code{normInts} Performs normalization for each set \emph{independently}.
#'
#'   }
#'
#'   A re-grouping of features occurs if \code{selectIons} is called or \code{adducts<-} is used with
#'   \code{reGroup=TRUE}. Afterwards, it is very likely that feature group names are changed. Since data generated later
#'   in the workflow (\emph{e.g.} annotation steps) rely on feature group names, these objects are \strong{not valid}
#'   anymore, and \strong{must} be re-generated.
#'
#' @rdname featureGroups-class
#' @export
featureGroupsSet <- setClass("featureGroupsSet",
                             slots = c(groupAlgo = "character", groupArgs = "list", groupVerbose = "logical",
                                       annotationsChanged = "list"),
                             contains = "featureGroups")

#' @rdname featureGroups-class
#' @export
setMethod("sets", "featureGroupsSet", function(obj) sets(getFeatures(obj)))

#' @rdname featureGroups-class
#' @export
setMethod("internalStandardAssignments", "featureGroupsSet", function(fGroups, set = NULL)
{
    if (is.null(set))
        return(fGroups@ISTDAssignments)
    assertSets(fGroups, set, FALSE)
    return(fGroups@ISTDAssignments[[set]])
})

#' @rdname featureGroups-class
#' @export
setMethod("adducts", "featureGroupsSet", function(obj, set, ...)
{
    assertSets(obj, set, FALSE)
    s <- set # workaround for below to distinguish between DT and parent set
    ann <- annotations(obj)[set == s]
    return(setNames(ann$adduct, ann$group))
})

#' @rdname featureGroups-class
#' @param reGroup \setsWF Set to \code{TRUE} to re-group the features after the adduct annotations are changed. See the
#'   \verb{Sets workflow} section for more details.
#' @export
setMethod("adducts<-", "featureGroupsSet", function(obj, value, set, reGroup = TRUE)
{
    # UNDONE: this function definition gives a warning with R checking as value must be the last argument, however,
    # putting it as last will not work with the other arguments as their names must then be specified explicitly and
    # defaults don't work.
    
    assertSets(obj, set, FALSE)
    checkmate::assertFlag(reGroup)
    
    s <- set
    ann <- annotations(obj)[set == s]
    
    checkmate::assertCharacter(value, min.chars = 1, any.missing = FALSE, len = nrow(ann))
    
    if (checkmate::testNamed(value))
    {
        checkmate::assertNames(names(value), permutation.of = ann$group, .var.name = "value")
        value <- value[ann$group] # ensure correct order
    }
    else
        names(value) <- ann$group

    updatedAnn <- updateAnnAdducts(ann, groupInfo(obj), value)
    
    if (!isTRUE(all.equal(annotations(obj), updatedAnn)))
    {
        if (length(obj@annotationsChanged) == 0)
            obj@annotationsChanged <- copy(obj@annotations)
        obj@annotationsChanged <- rbind(obj@annotationsChanged[set != s], updatedAnn)
        
        if (reGroup)
        {
            usFGroups <- sapply(sets(obj), unset, obj = obj, simplify = FALSE)
            usFGroups <- Map(usFGroups, sets(obj), f = function(ufg, us)
            {
                ufg@annotations <- obj@annotationsChanged[set == us]
                return(ufg)
            })
            obj <- do.call(makeSet, c(unname(usFGroups), list(groupAlgo = obj@groupAlgo, groupArgs = obj@groupArgs,
                                                              verbose = obj@groupVerbose, labels = names(usFGroups),
                                                              adducts = NULL)))
            obj@annotationsChanged <- data.table()
        }
    }
    
    return(obj)
})

#' @rdname featureGroups-class
#' @export
setMethod("delete", "featureGroupsSet", function(obj, i = NULL, j = NULL, ...)
{
    # HACK: subset annotations/ISTD assignments here as format with sets is different
    ann <- annotations(obj)
    if (nrow(ann) > 0)
        obj@annotations <- data.table() # disable subsetting in parent method
    ISTDAssign <- internalStandardAssignments(obj)
    if (length(ISTDAssign) > 0)
        obj@ISTDAssignments <- list()
    
    obj <- callNextMethod(obj, i, j, ...)
    
    if (nrow(ann) > 0)
        obj@annotations <- ann[set %in% sets(obj) & group %in% names(obj)]
    if (length(ISTDAssign) > 0)
        obj@ISTDAssignments <- lapply(ISTDAssign, function(ia) ia[names(ia) %chin% names(obj)])

    mySets <- sets(obj)
    if (nrow(internalStandards(obj)) > 0)
        obj@ISTDs <- obj@ISTDs[sapply(sets, function(s) any(mySets %chin% unlist(strsplit(s, ",", fixed = TRUE))))]
    if (nrow(toxicities(obj)) > 0)
        obj@toxicities <- obj@toxicities[sapply(sets, function(s) any(mySets %chin% unlist(strsplit(s, ",", fixed = TRUE))))]
    
    return(obj)
})

#' @rdname featureGroups-class
#' @export
setMethod("show", "featureGroupsSet", function(object)
{
    callNextMethod(object)
    printf("Sets: %s\n", paste0(sets(object), collapse = ", "))
})

#' @rdname featureGroups-class
#' @export
setMethod("featureTable", "featureGroupsSet", function(obj) featureTable(obj@features))

#' @rdname featureGroups-class
#' @export
setMethod("[", c("featureGroupsSet", "ANY", "ANY", "missing"), function(x, i, j, ..., ni, rGroups, sets = NULL,
                                                                        reorder = FALSE, drop = TRUE)
{
    assertSets(x, sets, TRUE)
    if (!is.null(sets))
        i <- mergeAnaSubsetArgWithSets(i, sets, analysisInfo(x), reorder)
    
    return(callNextMethod(x, i, j, ..., ni = ni, reorder = reorder, rGroups = rGroups))
})

# UNDONE: mention that object will be unset
#' @rdname featureGroups-class
#' @export
setMethod("export", "featureGroupsSet", function(obj, type, out, set) export(unset(obj, set), type, out))

#' @rdname featureGroups-class
#' @export
setMethod("unique", "featureGroupsSet", function(x, which, ..., sets = FALSE)
{
    checkmate::assertFlag(sets)
    if (sets)
    {
        mySets <- get("sets", pos = 2)(x)
        checkmate::assertSubset(which, mySets, empty.ok = FALSE, add = ac)
        ai <- analysisInfo(x)
        which <- unique(ai[set %in% which]$group)
    }
    callNextMethod(x, which = which, ...)
})

#' @rdname featureGroups-class
#' @export
setMethod("overlap", "featureGroupsSet", function(fGroups, which, exclusive, sets = FALSE)
{
    mySets <- get("sets", pos = 2)(fGroups)
    
    checkmate::assertFlag(sets)
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assert(checkmate::checkSubset(which, if (sets) mySets else replicateGroups(fGroups), empty.ok = FALSE),
                      checkmate::checkList(which, "character", any.missing = FALSE),
                      .var.name = "which", add = ac)
    checkmate::assertFlag(exclusive, add = ac)
    checkmate::reportAssertions(ac)
    
    if (sets)
    {
        if (is.list(which))
            stop("which cannot be a list when sets=TRUE", call. = FALSE)
        
        if (length(which) < 2 || length(fGroups) == 0)
            return(fGroups) # nothing to do...
        
        if (exclusive)
            ret <- unique(fGroups, which = which, sets = TRUE)
        else
            ret <- fGroups[, sets = which]
        
        ret <- minSetsFGroupsFilter(ret, relThreshold = 1, verbose = FALSE)
    }
    else
        ret <- callNextMethod(fGroups, which = which, exclusive = exclusive)
    
    return(ret)
})

#' @rdname featureGroups-class
#' @export
setMethod("selectIons", "featureGroupsSet", function(fGroups, components, prefAdduct, ...)
{
    setLen <- length(sets(fGroups))
    
    checkmate::assertClass(components, "componentsSet")
    checkmate::assert(checkmate::checkCharacter(prefAdduct, any.missing = FALSE, min.len = 1, max.len = setLen),
                      checkmate::checkList(prefAdduct, types = c("adduct", "character"), any.missing = FALSE,
                                           min.len = 1, max.len = setLen),
                      .var.name = "prefAdduct")
    prefAdduct <- rep(prefAdduct, length.out = setLen)
    
    if (length(components) == 0 || is.null(componentInfo(components)[["neutral_mass"]]))
    {
        cat("Components are empty or lack adduct/isotope annotations, skipping...\n")
        return(fGroups)
    }
    
    # annotate and merge ions for all set objects
    usFGroups <- sapply(sets(fGroups), unset, obj = fGroups, simplify = FALSE)
    usComponents <- checkAndUnSetOther(sets(fGroups), components, "components")
    usFGroups <- Map(usFGroups, usComponents, prefAdduct, f = selectIons, MoreArgs = list(...))
    
    # and re-group with new adduct information
    return(do.call(makeSet, c(unname(usFGroups), list(groupAlgo = fGroups@groupAlgo, groupArgs = fGroups@groupArgs,
                   verbose = fGroups@groupVerbose, labels = names(usFGroups), adducts = NULL))))
})

#' @rdname featureGroups-class
#' @export
setMethod("normInts", "featureGroupsSet", function(fGroups, featNorm, groupNorm, normFunc, standards, ISTDRTWindow,
                                                   ISTDMZWindow, minISTDs, ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertSubset(featNorm, c("tic", "istd", "conc", "none"))
    checkmate::assertFlag(groupNorm, add = ac)
    checkmate::reportAssertions(ac)
    
    if (length(fGroups) == 0)
        return(fGroups)
    
    baseArgs <- c(list(featNorm = featNorm, groupNorm = groupNorm, normFunc = normFunc, ISTDRTWindow = ISTDRTWindow,
                       ISTDMZWindow = ISTDMZWindow, minISTDs = minISTDs), list(...))
    if (featNorm != "istd" && !groupNorm)
        return(do.call(callNextMethod, c(list(fGroups, standards = NULL), baseArgs))) # no need to do per set

    # reset
    fGroups@ISTDs <- data.table()
    fGroups@ISTDAssignments <- list()
    
    usFGroups <- sapply(sets(fGroups), unset, obj = fGroups, simplify = FALSE)
    
    if (featNorm == "istd")
    {
        # NOTE: skipInvalid is set to FALSE. We don't have readily access to it (it's in ...), but if the user set it to
        # TRUE, it will be picked up later by screenSuspects() in the non sets method anyway.
        standards <- assertAndPrepareSuspectsSets(standards, sets(fGroups), skipInvalid = FALSE)
        usFGroups <- Map(usFGroups, standards = standards, f = normInts, MoreArgs = baseArgs)
        
        # merge ISTD slots
        fGroups@ISTDs <- mergeScreeningSetInfos(usFGroups, lapply(usFGroups, internalStandards))
        fGroups@ISTDAssignments <- lapply(usFGroups, slot, "ISTDAssignments")
    }
    else
        usFGroups <- do.call(lapply, c(list(usFGroups, normInts), baseArgs))
    
    
    allNormFeats <- Reduce(modifyList, lapply(usFGroups, featureTable))
    fGroups@features@features <- Map(featureTable(fGroups), allNormFeats[analyses(fGroups)], f = function(ft, ftN)
    {
        ft <- copy(ft)
        if (nrow(ft) == 0)
            ft[, c("intensity_rel", "area_rel") := numeric()]
        else
            ft[match(ftN$group, group), c("intensity_rel", "area_rel") := .(ftN$intensity_rel, ftN$area_rel)]
        return(ft)
    })
    
    return(fGroups)
})

#' @rdname pred-quant
#' @export
setMethod("calculateConcs", "featureGroupsSet", function(fGroups, featureAnn, areas = FALSE)
{
    # set null.ok to TRUE here to allow calculations from screening results. The non-sets methods called by
    # doCalcConcSets will assert !NULL if needed.
    checkmate::assertClass(featureAnn, c("featureAnnotations", "workflowStepSet"), null.ok = TRUE)
    
    return(doCalcConcSets(fGroups, featureAnn, areas))
})

#' @rdname pred-tox
#' @export
setMethod("calculateTox", "featureGroupsSet", function(fGroups, featureAnn)
{
    # set null.ok to TRUE here to allow calculations from screening results. The non-sets methods called by
    # doCalcToxSets will assert !NULL if needed.
    checkmate::assertClass(featureAnn, c("featureAnnotations", "workflowStepSet"), null.ok = TRUE)
    
    return(doCalcToxSets(fGroups, featureAnn))
})

#' @return The \code{featuresSet} method (for \link[=sets-workflow]{sets workflows}) returns a
#'   \code{\link{featureGroupsSet}} object.
#' @rdname groupFeatures
#' @export
setMethod("groupFeatures", "featuresSet", function(obj, algorithm, ..., verbose = TRUE)
{
    checkmate::assertChoice(algorithm, c("openms", "xcms", "xcms3", "kpic2"))
    
    otherArgs <- list(...)
    
    fGroups <- do.call(callNextMethod, c(list(obj = obj, algorithm = algorithm, verbose = verbose), otherArgs))
    obj <- getFeatures(fGroups) # may have been changed (eg in initialize())
    
    ret <- featureGroupsSet(groupAlgo = algorithm, groupArgs = otherArgs, groupVerbose = verbose,
                            groups = groupTable(fGroups), groupInfo = groupInfo(fGroups), features = obj,
                            ftindex = groupFeatIndex(fGroups), algorithm = makeSetAlgorithm(list(fGroups)))
    ret@annotations <- getAnnotationsFromSetFeatures(ret)
    
    return(ret)
})

#' @param groupAlgo groupAlgo The name of the feature grouping algorithm. See the \code{algorithm} argument of
#'   \code{\link{groupFeatures}} for details.
#' @param groupArgs A \code{list} with arguments directly passed to \code{groupFeatures} (can be named). Example:
#'   \code{groupArgs=list(maxAlignMZ=0.002)}.
#' @param verbose If set to \code{FALSE} then no text output is shown.
#'
#' @rdname makeSet
#' @export
setMethod("makeSet", "featureGroups", function(obj, ..., groupAlgo, groupArgs = NULL, verbose = TRUE,
                                               adducts = NULL, labels = NULL)
{
    if (is.null(labels) && is.null(adducts))
        stop("The labels and adducts arguments are not set (NULL). ",
             "Please set the labels argument, as automatic labelling requires adducts.")
    
    fGroupsList <- list(obj, ...)
    ac <- checkmate::makeAssertCollection()
    assertMakeSetArgs(fGroupsList, "featureGroups", adducts, TRUE, labels, ac)
    checkmate::assertList(groupArgs, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (is.null(groupArgs))
        groupArgs <- list()
    
    if (all(lengths(fGroupsList) == 0))
        stop("Cannot make set if all feature groups objects are empty")
    
    allAnas <- unlist(lapply(fGroupsList, analyses))
    if (anyDuplicated(allAnas))
        stop("Some objects have non-unique analyses: ", paste0(unique(allAnas[duplicated(allAnas)]), collapse = ","))
    
    if (!is.null(adducts))
    {
        adducts <- prepareMakeSetAdducts(fGroupsList, adducts, labels)
        adductsChr <- lapply(adducts, as.character)
        names(fGroupsList) <- names(adducts)
    }
    else
    {
        names(fGroupsList) <- labels
        for (i in seq_along(fGroupsList))
        {
            if (length(fGroupsList[[i]]) > 0 && nrow(annotations(fGroupsList[[i]])) == 0)
                stop("Missing feature ion annotations. Either set the adducts argument or run selectIons()")
        }
        adducts <- adductsChr <- setNames(rep(list(NULL), length(fGroupsList)), names(fGroupsList))
    }
    
    assignAdductsToFeatTab <- function(ft, mcol, add, ann)
    {
        ft <- copy(ft)
        if (nrow(ft) == 0)
            ft[, adduct := character()]
        else if (!is.null(add))
            ft[, adduct := add]
        else
            ft[, adduct := ann[match(ft[[mcol]], group)]$adduct]
        return(ft)
    }
    
    # prepare features: add adducts needed for neutralization later
    fGroupsList <- Map(fGroupsList, adductsChr, f = function(fGroups, add)
    {
        ann <- annotations(fGroups)
        fGroups@features@features <- lapply(featureTable(fGroups), assignAdductsToFeatTab, "group", add, ann)
        return(fGroups)
    })
    
    allFeats <- sapply(fGroupsList, getFeatures, simplify = FALSE)
    featSet <- doMakeFeaturesSet(allFeats, adducts)
    
    # convert all unset fGroups to pseudo features
    # neutralize pseudo features & group
    # assign new group names to featSet
    # construct final fGroups
    
    fgFeat <- convertFGroupsToPseudoFeatures(fGroupsList)
    fgFeat@features <- Map(featureTable(fgFeat), adductsChr, lapply(fGroupsList, annotations), f = assignAdductsToFeatTab,
                           MoreArgs = list(mcol = "ID"))
    fgFeat <- neutralizeFeatures(fgFeat, adduct = NULL)
    setGroups <- groupPseudoFeatures(fgFeat, groupAlgo, c(groupArgs, list(verbose = verbose)))
    featSet@features <- Map(featureTable(featSet), analysisInfo(featSet)$set, f = function(ft, s)
    {
        ft <- copy(ft)
        ft[, groupOld := group] # UNDONE: keep this?
        ft[, group := featureTable(setGroups)[[s]][match(ft$group, ID)]$group]
        return(ft)
    })
    
    # HACK: add a dummy column for empty feature tables to ensure a row is still added
    grpInts <- rbindlist(lapply(featureTable(featSet), function(ft)
    {
        if (nrow(ft) == 0)
            return(data.table(dummy = NA))
        setnames(transpose(ft[, "intensity", with = FALSE]), ft$group)
    }), fill = TRUE)
    if (!is.null(grpInts[["dummy"]]))
        grpInts[, dummy := NULL]
    setnafill(grpInts, fill = 0)
    setcolorder(grpInts, rownames(groupInfo(setGroups)))
    
    grpFeatInds <- rbindlist(lapply(featureTable(featSet), function(ft)
    {
        if (nrow(ft) == 0)
            return(data.table(dummy = NA))
        setnames(as.data.table(transpose(list(seq_len(nrow(ft))))), ft$group)
    }), fill = TRUE)
    if (!is.null(grpFeatInds[["dummy"]]))
        grpFeatInds[, dummy := NULL]
    setnafill(grpFeatInds, fill = 0)
    setcolorder(grpFeatInds, rownames(groupInfo(setGroups)))
    
    ret <- featureGroupsSet(groupAlgo = groupAlgo, groupArgs = groupArgs, groupVerbose = verbose,
                            groups = grpInts, groupInfo = groupInfo(setGroups), features = featSet,
                            ftindex = grpFeatInds, algorithm = paste0(groupAlgo, "-set"))
    ret@annotations <- getAnnotationsFromSetFeatures(ret)
    
    return(ret)
})

#' @rdname makeSet
#' @export
setMethod("makeSet", "featureGroupsSet", function(obj, ...)
{
    stop("Making a set from set objects is not supported", call. = FALSE)
})

#' @rdname featureGroups-class
#' @export
featureGroupsUnset <- setClass("featureGroupsUnset", contains = "featureGroups")

#' @rdname featureGroups-class
#' @export
setMethod("unset", "featureGroupsSet", function(obj, set)
{
    # UNDONE: mention that group names remain the same and thus represent neutral masses
    # UNDONE: or rename?
    
    assertSets(obj, set, FALSE)
    obj <- obj[, sets = set]
    
    gInfo <- groupInfo(obj)
    ann <- copy(annotations(obj))
    if (nrow(gInfo) > 0)
    {
        adducts <- sapply(unique(ann$adduct), as.adduct)
        gInfo$mzs <- calculateMasses(gInfo$mzs, adducts[ann$adduct], type = "mz")
        ann <- ann[, -"set"]
    }
    ISTDs <- copy(internalStandards(obj))
    if (nrow(ISTDs) > 0)
        ISTDs <- ISTDs[, -"sets"]
    ISTDAssign <- if (length(internalStandardAssignments(obj)) > 0) internalStandardAssignments(obj, set) else list()

    tox <- copy(toxicities(obj))
    if (nrow(tox) > 0)
        tox <- tox[, -"sets"]
    
    
    return(featureGroupsUnset(groups = copy(groupTable(obj)), groupInfo = gInfo,
                              features = unset(getFeatures(obj), set), ftindex = copy(groupFeatIndex(obj)),
                              groupQualities = copy(groupQualities(obj)), groupScores = copy(groupScores(obj)),
                              annotations = ann, ISTDs = ISTDs, ISTDAssignments = ISTDAssign,
                              concentrations = copy(concentrations(obj)), toxicities = tox,
                              algorithm = paste0(algorithm(obj), "_unset")))
})
