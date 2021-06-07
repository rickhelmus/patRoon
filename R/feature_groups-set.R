#' @include main.R
#' @include feature_groups.R
#' @include features-set.R
NULL

minSetsFGroupsFilter <- function(fGroups, absThreshold = 0, relThreshold = 0, negate = FALSE, verbose = TRUE)
{
    threshold <- getHighestAbsValue(absThreshold, relThreshold, length(sets(fGroups)))
    if (threshold == 0)
        return(fGroups)
    
    setsAna <- analysisInfo(fGroups)$set
    return(doFGroupsFilter(fGroups, "minimum sets", c(threshold, negate), function(fGroups)
    {
        pred <- function(x) length(unique(setsAna[x > 0])) >= threshold
        if (negate)
            pred <- Negate(pred)
        
        return(fGroups[, sapply(groupTable(fGroups), pred, USE.NAMES = FALSE)])
    }, "minSets", verbose))
}

#' @param set \setsWF The name of the set.
#' @param sets \setsWF For \code{[} and \code{filter}: a \code{character} with name(s) of the sets to keep (or remove if
#'   \code{negate=TRUE}).
#'
#'   For \code{plotInt}: if \code{TRUE} then feature intensities are plot per set (order follows the
#'   \link[=analysis-information]{analysis information}).
#'
#'   For \code{plotVenn}, \code{overlap} and \code{unique}: If \code{TRUE} then the \code{which} argument changes its
#'   meaning and is used to specify the names of the sets to be compared.
#' @param absMinSets,relMinSets \setsWF Feature groups are only kept when they contain data for at least this (absolute
#'   or relative) amount of sets. Set to \code{NULL} to ignore.
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
#'   \item \code{filter} and the subset operator (\code{[}) have specific arguments to choose/filter by (feature
#'   presence in) sets. See the \code{sets} argument description.
#'
#'   \item \code{as.data.table}: normalization of intensities is performed per set.
#'
#'   \item \code{export} Only allows to export data from one set. The \code{unset} method is used prior to exporting the
#'   data.
#'
#'   \item \code{overlap}, \code{unique}, \code{plotVenn}, \code{plotInt} allow to handle data per set. See the
#'   \code{sets} argument description.
#'
#'   \item \code{selectIons} Will perform a re-grouping of features. The implications of this are discussed below.
#'
#'   \item \code{makeSet} Currently not yet supported for set objects.
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
#' @export
setMethod("adducts<-", "featureGroupsSet", function(obj, value, set, reGroup = TRUE)
{
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

#' @describeIn featureGroupsSet Obtain feature information (see \code{\link{features}}).
#' @export
setMethod("featureTable", "featureGroupsSet", function(obj) featureTable(obj@features))

#' @describeIn featureGroupsSet Subset on analyses/feature groups.
#' @param rGroups An optional \code{character} vector: if specified only keep
#'   results for the given replicate groups (equivalent to the \code{rGroups}
#'   argument to \code{\link[=filter,featureGroups-method]{filter}}).
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
setMethod("filter", "featureGroupsSet", function(obj, ..., negate = FALSE, sets = NULL, absMinSets = NULL,
                                                 relMinSets = NULL)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFlag(negate, add = ac)
    assertSets(obj, sets, TRUE, add = ac)
    aapply(checkmate::assertNumber, . ~ absMinSets + relMinSets, lower = 0, finite = TRUE, null.ok = TRUE,
           fixed = list(add = ac))
    checkmate::reportAssertions(ac)

    if (!is.null(sets) && length(sets) > 0)
    {
        if (negate)
            sets <- setdiff(get("sets", pos = 2)(obj), sets)
        obj <- obj[, sets = sets]
    }

    if (...length() > 0)
        obj <- callNextMethod(obj, ..., negate = negate)
    
    if (!is.null(absMinSets) || !is.null(relMinSets))
        obj <- minSetsFGroupsFilter(obj, absMinSets, relMinSets, negate = negate)

    return(obj)
})

#' @rdname featureGroups-class
#' @export
setMethod("plotInt", "featureGroupsSet", function(obj, average = FALSE, xnames = !sets, showLegend = sets, pch = 20,
                                                  type = "b", lty = 3, col = NULL, ..., sets = FALSE)
{
    aapply(checkmate::assertFlag, . ~ average + xnames + showLegend + sets)
    
    if (!sets)
        return(callNextMethod(obj, average, xnames, showLegend, pch, type, lty, col, ...))
    else if (xnames)
        warning("xnames option is ignored if sets=TRUE")

    if (length(obj) == 0)
    {
        noDataPlot()
        return(invisible(NULL))
    }

    anaInfo <- analysisInfo(obj)    
    if (average)
    {
        gTable <- copy(averageGroups(obj))
        gTable[, set := anaInfo[match(replicateGroups(obj), anaInfo$group), "set"]]
    }
    else
    {
        gTable <- copy(groupTable(obj))
        gTable[, set := anaInfo$set]
    }
    
    if (is.null(col))
        col <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(sets(obj)))
    
    oldp <- par(no.readonly = TRUE)
    if (showLegend)
    {
        makeLegend <- function(x, y, ...)
        {
            return(legend(x, y, sets(obj), col = col, pch = pch, text.col = col, xpd = NA, ncol = 1,
                          cex = 0.75, bty = "n", ...))
        }
        
        plot.new()
        leg <- makeLegend(0, 0, plot = FALSE)
        lw <- (grconvertX(leg$rect$w, to = "ndc") - grconvertX(0, to = "ndc"))
        par(omd = c(0, 1 - lw, 0, 1), new = TRUE)
    }
    
    maxs <- max(table(gTable$set))
    plot(x = c(0, maxs), y = c(0, max(gTable[, -"set"])), type = "n", xlab = "", ylab = "Intensity", xaxt = "n")
    axis(1, seq_len(maxs), seq_len(maxs))
    
    for (s in seq_along(sets(obj)))
    {
        gt <- gTable[set == sets(obj)[s]]
        for (g in names(obj))
            lines(x = seq_len(nrow(gt)), y = gt[[g]], type = type, pch = pch, lty = lty, col = col[s], ...)
    }
    
    if (showLegend)
        makeLegend(par("usr")[2], par("usr")[4])
    
    par(oldp)
})

#' @rdname featureGroups-class
#' @export
setMethod("plotVenn", "featureGroupsSet", function(obj, which = NULL, ..., sets = FALSE)
{
    checkmate::assertFlag(sets)
    if (sets)
    {
        mySets <- get("sets", pos = 2)(obj)
        ai <- analysisInfo(obj)
        which = sapply(mySets, function(s) ai[ai$set == s, "group"], simplify = FALSE)
    }
    callNextMethod(obj, which = which, ...)
})

#' @rdname featureGroups-class
#' @export
setMethod("unique", "featureGroupsSet", function(x, which, ..., sets = FALSE)
{
    checkmate::assertFlag(sets)
    if (sets)
    {
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

#' @rdname featureGroups-class
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

#' @rdname featureGroups-class
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
    ann <- annotations(obj)
    if (nrow(gInfo) > 0)
    {
        adducts <- sapply(unique(ann$adduct), as.adduct)
        addMZs <- sapply(adducts, adductMZDelta)
        gInfo$mzs <- gInfo$mzs + addMZs[ann$adduct]
        ann <- ann[, -"set"]
    }
    
    return(featureGroupsUnset(groups = groupTable(obj), groupInfo = gInfo,
                              analysisInfo = unSetAnaInfo(analysisInfo(obj)),
                              features = unset(getFeatures(obj), set), ftindex = groupFeatIndex(obj),
                              annotations = ann, algorithm = paste0(algorithm(obj), "_unset")))
})
