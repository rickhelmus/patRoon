#' @include main.R
#' @include feature_groups.R
#' @include features-set.R
NULL

featureGroupsSet <- setClass("featureGroupsSet",
                             slots = c(groupAlgorithm = "character"),
                             contains = "featureGroups")

setMethod("initialize", "featureGroupsSet",
          function(.Object, ...) callNextMethod(.Object, algorithm = "set", ...))

setMethod("sets", "featureGroupsSet", function(obj) sets(getFeatures(obj)))
setMethod("adducts", "featureGroupsSet", function(obj) adducts(getFeatures(obj)))
setMethod("groupAlgorithm", "featureGroupsSet", function(obj) obj@groupAlgorithm)

#' @describeIn featureGroupsSet Shows summary information for this object.
#' @export
setMethod("show", "featureGroupsSet", function(object)
{
    callNextMethod(object)
    printf("sets: %s\n", paste0(sets(object), collapse = ", "))
    printf("adducts: %s\n", paste0(sapply(adducts(getFeatures(object)), as.character), collapse = ", "))
    printf("grouping algorithm: %s\n", groupAlgorithm(object))
})

#' @describeIn featureGroupsSet Obtain feature information (see \code{\link{features}}).
#' @export
setMethod("featureTable", "featureGroupsSet", function(obj, neutralized = TRUE, set = NULL) featureTable(obj@features, neutralized, set))

#' @describeIn featureGroupsSet Subset on analyses/feature groups.
#' @param rGroups An optional \code{character} vector: if specified only keep
#'   results for the given replicate groups (equivalent to the \code{rGroups}
#'   argument to \code{\link[=filter,featureGroups-method]{filter}}).
#' @export
setMethod("[", c("featureGroupsSet", "ANY", "ANY", "missing"), function(x, i, j, ..., sets = NULL, drop = TRUE)
{
    assertSets(x, sets)
    
    if (!is.null(sets) && length(sets) > 0)
        i <- mergeAnaSubsetArgWithSets(i, sets, analysisInfo(x))

    return(callNextMethod(x, i, j, ...))
})

#' @describeIn featureGroupsSet Extract intensity values.
#' @export
setMethod("[[", c("featureGroupsSet", "ANY", "ANY"), function(x, i, j, sets = NULL)
{
    ac <- checkmate::makeAssertCollection()
    assertExtractArg(i, add = ac)
    assertSets(x, sets, add = ac)
    checkmate::reportAssertions(ac)

    if (!is.null(sets) && length(sets) > 0)
        i <- mergeAnaSubsetArgWithSets(i, sets, analysisInfo(x))
    
    return(callNextMethod(x, i, j))
})

# UNDONE: mention that object will be ionized
#' @describeIn featureGroupsSet Exports feature groups to a \file{.csv} file that
#'   is readable to Bruker ProfileAnalysis (a 'bucket table'), Bruker TASQ (an
#'   analyte database) or that is suitable as input for the \verb{Targeted peak
#'   detection} functionality of \href{http://mzmine.github.io/}{MZmine}.
#' @param out The destination file for the exported data.
#' @export
setMethod("export", "featureGroupsSet", function(obj, type, out, sets = NULL) callNextMethod(ionize(obj, sets), type, out))

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
setMethod("as.data.table", "featureGroupsSet", function(x, neutralized = TRUE, sets = NULL, features = FALSE, ...)
{
    # UNDONE: also support reporting ionized features with different adducts?
    
    assertSets(x, sets)
    
    if (!is.null(sets) && length(sets) > 0)
        x <- x[, sets = sets]

    anaInfo <- analysisInfo(x) # get before ionizing    
    if (!neutralized)
        x <- ionize(x)
    
    ret <- callNextMethod(x, features = features, ...)
    
    if (features) # add set column
    {
        ret[, set := anaInfo[match(analysis, anaInfo$analysis), "set"]]
        setcolorder(ret, c("group", "group_ret", "group_mz", "set", "analysis"))
    }
    
    return(ret[])
})

setMethod("filter", "featureGroupsSet", function(obj, ..., sets = NULL, negate = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    assertSets(obj, sets, add = ac)
    checkmate::assertFlag(negate, add = ac)
    checkmate::reportAssertions(ac)

    if (!is.null(sets) && length(sets) > 0)
    {
        if (negate)
            sets <- setdiff(get("sets", pos = 2)(obj), sets)
        obj <- obj[, sets = sets]
    }

    if (length(list(...)) > 0)
        return(callNextMethod(obj, ..., negate = negate))
    return(obj)
})

#' @export
setMethod("plotEIC", "featureGroupsSet", function(obj, ..., retMin = FALSE, title = NULL)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFlag(retMin, add = ac)
    checkmate::assertString(title, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (is.null(title) && length(obj) == 1)
    {
        # override default title
        gInfo <- groupInfo(obj)
        title <- sprintf("Group '%s'\nrt: %.1f - neutralized mass: %.4f", names(obj)[1],
                         if (retMin) gInfo[1, "rts"] / 60 else gInfo[1, "rts"],
                         gInfo[1, "mzs"])
    }
    callNextMethod(obj, ..., retMin = retMin, title = title)
})

setMethod("groupFeatures", "featuresSet", function(feat, algorithm, ..., verbose = TRUE)
{
    # UNDONE: xcms3 not yet supported
    checkmate::assertChoice(algorithm, c("openms", "xcms"))
    
    otherArgs <- list(...)
    if (algorithm == "xcms")
        otherArgs <- modifyList(otherArgs, list(exportedData = FALSE))
    
    fGroups <- do.call(callNextMethod, c(list(feat = feat, algorithm = algorithm, verbose = verbose), otherArgs))
    
    return(featureGroupsSet(groupAlgorithm = algorithm, groups = groups(fGroups), groupInfo = groupInfo(fGroups),
                            analysisInfo = analysisInfo(fGroups), features = feat, ftindex = groupFeatIndex(fGroups)))
})

featureGroupsSetIonized <- setClass("featureGroupsSetIonized", contains = "featureGroups")
setMethod("initialize", "featureGroupsSetIonized",
          function(.Object, ...) callNextMethod(.Object, algorithm = "set_ionized", ...))
setMethod("ionize", "featureGroupsSet", function(obj, sets)
{
    # UNDONE: mention that group names remain the same and thus represent neutral masses
    
    if (!is.null(sets) && length(sets) > 0)
        obj <- obj[, sets = sets]
    
    if (!allSame(adducts(obj)))
        stop("Selected sets for conversion must have have equal adducts")
    
    addMZ <- adductMZDelta(adducts(obj)[[1]])
    gInfo <- groupInfo(obj)
    gInfo$mzs <- gInfo$mzs + addMZ
    
    return(featureGroupsSetIonized(groups = groups(obj), groupInfo = gInfo, analysisInfo = analysisInfo(obj),
                                   features = ionize(getFeatures(obj)), ftindex = groupFeatIndex(obj)))
})
