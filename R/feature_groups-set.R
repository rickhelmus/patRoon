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
#' @slot annotations \setsWF As the \code{featureGroups} slot, but contains the annotation data per set.
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
#'   discussed below.
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
                             slots = c(groupAlgo = "character", groupArgs = "list", groupVerbose = "logical"),
                             contains = "featureGroups")

#' @rdname featureGroups-class
#' @export
setMethod("sets", "featureGroupsSet", function(obj) sets(getFeatures(obj)))

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
        obj@annotations <- rbind(obj@annotations[set != s], updatedAnn)
        
        if (reGroup)
        {
            usFGroups <- sapply(sets(obj), unset, obj = obj, simplify = FALSE)
            obj <- do.call(makeSet, c(unname(usFGroups), list(groupAlgo = obj@groupAlgo, groupArgs = obj@groupArgs,
                                                              verbose = obj@groupVerbose, labels = names(usFGroups),
                                                              adducts = NULL)))
        }
    }
    
    return(obj)
})

#' @rdname featureGroups-class
#' @export
setMethod("delete", "featureGroupsSet", function(obj, i = NULL, j = NULL, ...)
{
    # HACK: subset annotations here as format with sets is different
    ann <- annotations(obj)
    if (nrow(ann) > 0)
        obj@annotations <- data.table() # disable subsetting in parent method
    
    obj <- callNextMethod(obj, i, j, ...)
    
    if (nrow(ann) > 0)
        obj@annotations <- ann[set %in% sets(obj) & group %in% names(obj)]
    
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
setMethod("[", c("featureGroupsSet", "ANY", "ANY", "missing"), function(x, i, j, ..., rGroups, sets = NULL, drop = TRUE)
{
    assertSets(x, sets, TRUE)
    if (!is.null(sets))
        i <- mergeAnaSubsetArgWithSets(i, sets, analysisInfo(x))
    
    return(callNextMethod(x, i, j, ..., rGroups = rGroups))
})

# UNDONE: mention that object will be unset
#' @rdname featureGroups-class
#' @export
setMethod("export", "featureGroupsSet", function(obj, type, out, set) export(unset(obj, set), type, out))

#' @rdname featureGroups-class
#' @export
setMethod("as.data.table", "featureGroupsSet", function(x, average = FALSE, areas = FALSE, features = FALSE,
                                                        qualities = FALSE, regression = FALSE, averageFunc = mean,
                                                        normFunc = NULL, FCParams = NULL)
{
    # NOTE keep args in sync with featureGroupsScreeningSet
    
    anaInfo <- analysisInfo(x)
    
    # HACK: add annotations later as format with sets is different
    ann <- x@annotations
    if (nrow(ann) > 0)
        x@annotations <- data.table()
    
    # NOTE: we normalize hereafter per set afterwards
    ret <- callNextMethod(x, average = average, areas = areas, features = features, qualities = qualities,
                          regression = regression, averageFunc = averageFunc, normFunc = NULL, FCParams = FCParams)
    
    if (!is.null(ret[["analysis"]])) # add set column if feature data is present
    {
        ret[, set := anaInfo[match(analysis, anaInfo$analysis), "set"]]
        setcolorder(ret, c("group", "group_ret", "group_mz", "set", "analysis"))
    }
    
    if (!is.null(normFunc))
    {
        # do normalization here to do so per set
        
        if (features)
            ret[, c("area", "intensity") := .(if (all(area == 0)) area else area / normFunc(area),
                                              if (all(intensity == 0)) intensity else intensity / normFunc(intensity)),
                by = c("set", "group")]
        else
        {
            rowSeq <- seq_len(nrow(ret)) # BUG? can't put in expression below directly...
            for (s in sets(x))
            {
                anaInfo <- analysisInfo(x)
                anaInfo <- anaInfo[anaInfo$set == s, ]
                intCols <- if (average) unique(anaInfo$group) else anaInfo$analysis
                ret[, (intCols) := {
                    v <- unlist(.SD)
                    if (all(v == 0)) .SD else as.list(v / normFunc(v))
                }, by = rowSeq, .SDcols = intCols]
            }
        }
    }
    
    if (nrow(ann) > 0)
    {
        if (features && !average)
            ret <- merge(ret, ann, by = c("group", "set"))
        else
        {
            # collapse annotation info for each group
            ann <- copy(ann)
            ann[, adduct := paste0(adduct, collapse = ","), by = "group"]
            ann <- unique(ann, by = "group")[, -"set"]
            ret <- merge(ret, ann, by = "group", sort = FALSE)
        }
    }
    
    return(ret[])
})

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
        which <- unique(ai[ai$set %in% which, "group"])
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
    checkmate::assertCharacter(which, min.len = 2, min.chars = 1, any.missing = FALSE, add = ac)
    checkmate::assertSubset(which, if (sets) mySets else replicateGroups(fGroups),
                            empty.ok = FALSE, add = ac)
    checkmate::assertFlag(exclusive, add = ac)
    checkmate::reportAssertions(ac)
    
    if (sets)
    {
        anaInfo <- analysisInfo(fGroups)
        
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
    
    # annotate and merge ions for all set objects
    usFGroups <- sapply(sets(fGroups), unset, obj = fGroups, simplify = FALSE)
    usComponents <- checkAndUnSetOther(sets(fGroups), components, "components")
    usFGroups <- Map(usFGroups, usComponents, prefAdduct, f = selectIons, MoreArgs = list(...))
    
    # and re-group with new adduct information
    return(do.call(makeSet, c(unname(usFGroups), list(groupAlgo = fGroups@groupAlgo, groupArgs = fGroups@groupArgs,
                   verbose = fGroups@groupVerbose, labels = names(usFGroups), adducts = NULL))))
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
                            groups = groupTable(fGroups), groupInfo = groupInfo(fGroups),
                            analysisInfo = analysisInfo(fGroups), features = obj, ftindex = groupFeatIndex(fGroups),
                            algorithm = makeSetAlgorithm(list(fGroups)))
    
    anaInfo <- analysisInfo(ret)
    ftind <- groupFeatIndex(ret)
    fTable <- featureTable(ret)

    if (length(obj) > 0)
    {
        ret@annotations <- rbindlist(sapply(sets(obj), function(s)
        {
            anaInds <- which(anaInfo$set == s)
            anas <- anaInfo[anaInds, "analysis"]
            grps <- names(ret)[sapply(ftind[anaInds], function(x) any(x != 0))]
            
            firstFeats <- rbindlist(lapply(ftind[anaInds, grps, with = FALSE], function(x)
            {
                firstAna <- which(x != 0)[1]
                return(featureTable(ret)[[anas[firstAna]]][x[firstAna]])
            }))
            
            return(data.table(group = grps, adduct = firstFeats$adduct))
        }, simplify = FALSE), idcol = "set", fill = TRUE) # set fill for empty objects
        ret@annotations[, neutralMass := groupInfo(ret)[ret@annotations$group, "mzs"]]
    }
    else
        ret@annotations <- data.table()
    
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
    
    # prepare features: add adducts needed for neutralization and clearout group assignments
    fGroupsList <- Map(fGroupsList, adductsChr, f = function(fGroups, add)
    {
        ann <- annotations(fGroups)
        
        fGroups@features@features <- lapply(featureTable(fGroups), function(ft)
        {
            ft <- copy(ft)
            if (nrow(ft) == 0)
                ft[, adduct := character()]
            else if (!is.null(add))
                ft[, adduct := add]
            else
                ft[, adduct := ann[match(ft$group, group)]$adduct]
            ft[, group := NULL]
            return(ft)
        })
        return(fGroups)
    })
    
    allFeats <- sapply(fGroupsList, getFeatures, simplify = FALSE)
    featSet <- doMakeFeaturesSet(allFeats, adducts)
    
    return(do.call(groupFeatures, c(list(featSet, algorithm = groupAlgo, verbose = verbose), groupArgs)))
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
    
    return(featureGroupsUnset(groups = copy(groupTable(obj)), groupInfo = gInfo,
                              analysisInfo = unSetAnaInfo(analysisInfo(obj)),
                              features = unset(getFeatures(obj), set), ftindex = copy(groupFeatIndex(obj)),
                              groupQualities = copy(groupQualities(obj)), groupScores = copy(groupScores(obj)),
                              annotations = ann, iSTDs = copy(internalStandards(fGroups)),
                              algorithm = paste0(algorithm(obj), "_unset")))
})
