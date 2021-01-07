#' @include main.R
#' @include feature_groups.R
#' @include features-set.R
NULL

featureGroupsSet <- setClass("featureGroupsSet",
                             contains = "featureGroups")

setMethod("sets", "featureGroupsSet", function(obj) sets(getFeatures(obj)))

setMethod("removeGroups", "featureGroupsSet", function(fGroups, indices)
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
setMethod("featureTable", "featureGroupsSet", function(obj, neutralized = TRUE) featureTable(obj@features, neutralized = neutralized))

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
                                                        regression = FALSE, normFunc = NULL)
{
    # NOTE keep args in sync with featureGroupsScreeningSet
    
    anaInfo <- analysisInfo(x)
    
    # NOTE: we normalize hereafter per set afterwards
    ret <- callNextMethod(x, average = average, areas = areas, features = features,
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
    
    return(ret[])
})

setMethod("filter", "featureGroupsSet", function(obj, ..., negate = FALSE, sets = NULL)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFlag(negate, add = ac)
    assertSets(obj, sets, TRUE, add = ac)
    checkmate::reportAssertions(ac)

    if (!is.null(sets) && length(sets) > 0)
    {
        if (negate)
            sets <- setdiff(get("sets", pos = 2)(obj), sets)
        obj <- obj[, sets = sets]
    }

    if (...length() > 0)
        return(callNextMethod(obj, ..., negate = negate))
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

setMethod("groupFeatures", "featuresSet", function(feat, algorithm, ..., verbose = TRUE)
{
    # UNDONE: xcms3 not yet supported
    checkmate::assertChoice(algorithm, c("openms", "xcms"))
    
    otherArgs <- list(...)
    if (algorithm == "xcms")
        otherArgs <- modifyList(otherArgs, list(exportedData = FALSE))
    
    fGroups <- do.call(callNextMethod, c(list(feat = feat, algorithm = algorithm, verbose = verbose), otherArgs))
    
    ret <- featureGroupsSet(groups = groupTable(fGroups), groupInfo = groupInfo(fGroups),
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
        
        return(data.table(group = grps, adduct = firstFeats$adduct, isonr = firstFeats$isonr))
    }, simplify = FALSE), idcol = "set")
    
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
                stop("Missing feature ion annotations. Either set the adducts argument or run mergeIons()")
        }
        adducts <- adductsChr <- setNames(rep(list(NULL), length(fGroupsList)), names(fGroupsList))
    }
    
    # prepare features: add adducts needed for neutralization and remove left-over features
    fGroupsList <- Map(fGroupsList, adductsChr, f = function(fGroups, add)
    {
        ftindAna <- transpose(groupFeatIndex(fGroups))
        ann <- annotations(fGroups)
        
        fGroups@features@features <- Map(featureTable(fGroups), ftindAna, f = function(ft, fti)
        {
            gInds <- which(fti != 0)
            fti <- fti[fti != 0]
            ft <- copy(ft)
            if (!is.null(add))
                ft[fti, adduct := add]
            else
                ft[fti, adduct := ann$adduct[gInds]]
            
            ft <- ft[fti] # remove features not in any group
            
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
