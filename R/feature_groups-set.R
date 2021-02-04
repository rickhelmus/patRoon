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

featureGroupsSet <- setClass("featureGroupsSet",
                             slots = c(groupAlgo = "character", groupArgs = "list", groupVerbose = "logical"),
                             contains = "featureGroups")

setMethod("sets", "featureGroupsSet", function(obj) sets(getFeatures(obj)))

#' @export
setMethod("adducts", "featureGroupsSet", function(obj, set)
{
    assertSets(obj, set, FALSE)
    s <- set # workaround for below to distinguish between DT and parent set
    ann <- annotations(obj)[set == s]
    return(setNames(ann$adduct, ann$group))
})

#' @export
setReplaceMethod("adducts", "featureGroupsSet", function(obj, value, set, reGroup = TRUE)
{
    assertSets(obj, set, FALSE)
    
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

setMethod("removeGroups", "featureGroupsSet", function(fGroups, indices, updateFeatures)
{
    # HACK: subset annotations here as format with sets is different
    ann <- fGroups@annotations
    if (nrow(ann) > 0)
        fGroups@annotations <- data.table() # disable subsetting in fGroups method
    
    fGroups <- callNextMethod()
    
    if (nrow(ann) > 0)
        fGroups@annotations <- ann[set %in% sets(fGroups) & group %in% names(fGroups)]
    
    return(fGroups)
})

#' @describeIn featureGroupsSet Shows summary information for this object.
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
    
    if (!is.null(sets) && length(sets) > 0)
    {
        i <- mergeAnaSubsetArgWithSets(i, sets, analysisInfo(x))
        x@annotations <- x@annotations[set %in% sets]
    }
    
    return(callNextMethod(x, i, j, ..., rGroups = rGroups))
})

# UNDONE: mention that object will be unset
#' @describeIn featureGroupsSet Exports feature groups to a \file{.csv} file that
#'   is readable to Bruker ProfileAnalysis (a 'bucket table'), Bruker TASQ (an
#'   analyte database) or that is suitable as input for the \verb{Targeted peak
#'   detection} functionality of \href{http://mzmine.github.io/}{MZmine}.
#' @param out The destination file for the exported data.
#' @export
setMethod("export", "featureGroupsSet", function(obj, type, out, set) export(unset(obj, set), type, out))

# UNDONE: mention that object will be unset
#' @export
setMethod("getXCMSSet", "featureGroupsSet", function(obj, ..., set) getXCMSSet(unset(obj, set), ...))

# UNDONE: mention that object will be unset
setMethod("getXCMSnExp", "featureGroupsSet", function(obj, ..., set) getXCMSnExp(unset(obj, set), ...))

#' @describeIn featureGroupsSet Obtain a summary table (a \code{\link{data.table}})
#'   with retention, \emph{m/z}, intensity and optionally other feature data.
#' @param features If \code{TRUE} then feature specific data will be added. If
#'   \code{average=TRUE} this data will be averaged for each feature group.
#' @param regression Set to \code{TRUE} to add regression data for each feature
#'   group. For this a linear model is created (intensity/area [depending on
#'   \code{areas} argument] \emph{vs} concentration). The model concentrations
#'   (e.g. of a set of standards) is derived from the \code{conc} column of the
#'   \link[=analysis-information]{analysis information}. From this model the
#'   intercept, slope and R2 is added to the output. In addition, when
#'   \code{features=TRUE}, concentrations for each feature are added. Note that
#'   no regression information is added when no \code{conc} column is present in
#'   the analysis information or when less than two concentrations are specified
#'   (\emph{i.e.} the minimum amount).
#' @export
setMethod("as.data.table", "featureGroupsSet", function(x, average = FALSE, areas = FALSE, features = FALSE,
                                                        qualities = FALSE, regression = FALSE, normFunc = NULL)
{
    # NOTE keep args in sync with featureGroupsScreeningSet
    
    anaInfo <- analysisInfo(x)
    
    # HACK: add annotations later as format with sets is different
    ann <- fGroups@annotations
    if (nrow(ann) > 0)
        x@annotations <- data.table()
    
    # NOTE: we normalize hereafter per set afterwards
    ret <- callNextMethod(x, average = average, areas = areas, features = features, qualities = qualities,
                          regression = regression, normFunc = NULL)
    
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
        if (features)
            ret <- merge(ret, ann, by = c("group", "set"))
        else
        {
            # collapse annotation info for each group
            ann <- copy(ann)
            ann[, adduct := paste0(adduct, collapse = ","), by = "group"]
            ann <- unique(ann, by = "group")[, -"set"]
            ret <- merge(ret, ann, by = "group")
        }
    }
    
    return(ret[])
})

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
        obj <- minSetsFGroupsFilter(absMinSets, relMinSets, negate = negate)

    return(obj)
})

#' @export
setMethod("plotChroms", "featureGroupsSet", function(obj, ...)
{
    # UNDONE: find a neat way to keep both argument order and also override title if not specified
    
    # ac <- checkmate::makeAssertCollection()
    # checkmate::assertFlag(retMin, add = ac)
    # checkmate::assertString(title, null.ok = TRUE, add = ac)
    # checkmate::reportAssertions(ac)
    # 
    # if (is.null(title) && length(obj) == 1)
    # {
    #     # override default title
    #     gInfo <- groupInfo(obj)
    #     title <- sprintf("Group '%s'\nrt: %.1f - neutralized mass: %.4f", names(obj)[1],
    #                      if (retMin) gInfo[1, "rts"] / 60 else gInfo[1, "rts"],
    #                      gInfo[1, "mzs"])
    # }
    
    callNextMethod(obj, ...)
})

#' @export
setMethod("plotVenn", "featureGroupsSet", function(obj, which = NULL, ..., sets = FALSE)
{
    checkmate::assertFlag(sets)
    if (sets)
    {
        mySets <- get("sets", pos = 2)(fGroups)
        ai <- analysisInfo(obj)
        which = sapply(mySets, function(s) ai[ai$set == s, "group"], simplify = FALSE)
    }
    callNextMethod(obj, which = which, ...)
})

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
    usComponents <- sapply(sets(components), unset, obj = components, simplify = FALSE)
    usComponents <- usComponents[sets(fGroups)]
    usFGroups <- Map(usFGroups, usComponents, prefAdduct, f = selectIons, MoreArgs = list(...))
    
    # and re-group with new adduct information
    return(do.call(makeSet, c(unname(usFGroups), list(groupAlgo = fGroups@groupAlgo, groupArgs = fGroups@groupArgs,
                   verbose = fGroups@groupVerbose, labels = names(usFGroups), adducts = NULL))))
})

setMethod("groupFeatures", "featuresSet", function(feat, algorithm, ..., verbose = TRUE)
{
    # UNDONE: xcms3 not yet supported
    checkmate::assertChoice(algorithm, c("openms", "xcms"))
    
    otherArgs <- list(...)
    if (algorithm == "xcms")
        otherArgs <- modifyList(otherArgs, list(exportedData = FALSE))
    
    fGroups <- do.call(callNextMethod, c(list(feat = feat, algorithm = algorithm, verbose = verbose), otherArgs))
    
    ret <- featureGroupsSet(groupAlgo = algorithm, groupArgs = otherArgs, groupVerbose = verbose,
                            groups = groupTable(fGroups), groupInfo = groupInfo(fGroups),
                            analysisInfo = analysisInfo(fGroups), features = feat, ftindex = groupFeatIndex(fGroups),
                            algorithm = makeSetAlgorithm(list(fGroups)))
    
    anaInfo <- analysisInfo(ret)
    ftind <- groupFeatIndex(ret)
    fTable <- featureTable(ret)
    
    ret@annotations <- rbindlist(sapply(sets(feat), function(s)
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
    }, simplify = FALSE), idcol = "set")
    ret@annotations[, neutralMass := groupInfo(ret)[ret@annotations$group, "mzs"]]
    
    return(ret)
})

#' @export
setMethod("makeSet", "featureGroups", function(obj, ..., groupAlgo, groupArgs = NULL, verbose = TRUE,
                                               adducts, labels)
{
    if (is.null(labels) && is.null(adducts))
        stop("The labels and adducts arguments are not set (NULL). ",
             "Please set the labels argument, as automatic labelling requires adducts.")
    
    fGroupsList <- list(obj, ...)
    ac <- checkmate::makeAssertCollection()
    assertMakeSetArgs(fGroupsList, "featureGroups", adducts, TRUE, labels, ac)
    checkmate::assertList(groupArgs, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
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
            if (nrow(annotations(fGroupsList[[i]])) == 0)
                stop("Missing feature ion annotations. Either set the adducts argument or run selectIons()")
        }
        adducts <- adductsChr <- setNames(rep(list(NULL), length(fGroupsList)), names(fGroupsList))
    }
    
    # prepare features: add adducts needed for neutralization and clearout group assignments
    fGroupsList <- Map(fGroupsList, adductsChr, f = function(fGroups, add)
    {
        ftindAna <- transpose(groupFeatIndex(fGroups))
        ann <- annotations(fGroups)
        
        fGroups@features@features <- lapply(featureTable(fGroups), function(ft)
        {
            ft <- copy(ft)
            if (!is.null(add))
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


featureGroupsUnset <- setClass("featureGroupsUnset", contains = "featureGroups")
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
    }
    
    return(featureGroupsUnset(groups = groupTable(obj), groupInfo = gInfo, analysisInfo = analysisInfo(obj),
                              features = unset(getFeatures(obj), set), ftindex = groupFeatIndex(obj),
                              annotations = ann[, -"set"], algorithm = paste0(algorithm(obj), "_unset")))
})
