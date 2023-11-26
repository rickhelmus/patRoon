#' @include main.R
#' @include feature_groups.R
NULL

#' @rdname featureGroups-class
#' @export
featureGroupsXCMS <- setClass("featureGroupsXCMS", slots = c(xs = "ANY"), contains = "featureGroups")

setMethod("initialize", "featureGroupsXCMS",
          function(.Object, ...) callNextMethod(.Object, algorithm = "xcms", ...))

#' Group features using XCMS (old interface)
#'
#' Group and align features with the legacy \code{\link[xcms]{xcmsSet}} function from the \pkg{xcms} package.
#'
#' @templateVar algo XCMS
#' @templateVar do group features
#' @templateVar generic groupFeatures
#' @templateVar algoParam xcms
#' @template algo_generator
#' 
#' @details Grouping of features and
#'   alignment of their retention times are performed with the \code{\link[xcms:group-methods]{xcms::group}} and
#'   \code{\link[xcms:retcor-methods]{xcms::retcor}} functions, respectively. Both functions have an extensive list of
#'   parameters to modify their behavior and may therefore be used to potentially optimize results.
#'
#' @template feat-arg
#' @template rtalign-arg
#' @template loadrawdata-arg
#' 
#' @param retcorArgs named \code{character vector} that may contain extra parameters to be used by
#'   \code{\link[xcms:retcor-methods]{xcms::retcor}}.
#' @param groupArgs named \code{character vector} that may contain extra grouping parameters to be used by
#'   \code{\link[xcms:group-methods]{xcms::group}}
#' 
#' @inheritParams groupFeatures
#'
#' @inherit groupFeatures return
#'
#' @template sets-loadrawdata-RTalign-note
#'
#' @references \addCitations{xcms}{1} \cr\cr \addCitations{xcms}{2} \cr\cr \addCitations{xcms}{3}
#'
#' @templateVar what groupFeaturesXCMS
#' @templateVar cl features
#' @template main-rd-method
#' @export
setMethod("groupFeaturesXCMS", "features", function(feat, rtalign = TRUE, loadRawData = TRUE,
                                                    groupArgs = list(mzwid = 0.015), 
                                                    retcorArgs = list(method = "obiwarp"), verbose = TRUE)
{
    # UNDONE: keep loadRawData things? Or just require that it's exported? If keep document also for OpenMS and implications.

    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertFlag, . ~ rtalign + loadRawData + verbose, fixed = list(add = ac))
    aapply(checkmate::assertList, . ~ groupArgs + retcorArgs, any.missing = FALSE, names = "unique", fixed = list(add = ac))
    checkmate::reportAssertions(ac)

    xs <- getXCMSSet(feat, verbose = verbose, loadRawData = loadRawData)
    return(doGroupFeaturesXCMS(xs, feat, rtalign, loadRawData, groupArgs, retcorArgs, verbose))
})

#' @rdname groupFeaturesXCMS
#' @export
setMethod("groupFeaturesXCMS", "featuresSet", function(feat, groupArgs = list(mzwid = 0.015), verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFlag(verbose, add = ac)
    checkmate::assertList(groupArgs, any.missing = FALSE, names = "unique", add = ac)
    checkmate::reportAssertions(ac)
    
    # HACK: force non-set features method to allow grouping of neutralized features
    # UNDONE: or simply export this functionality with a flag?
    xs <- selectMethod("getXCMSSet", "features")(feat, verbose = verbose, loadRawData = FALSE)
    
    return(doGroupFeaturesXCMS(xs, feat, rtalign = FALSE, loadRawData = FALSE, groupArgs, list(), verbose))
})

doGroupFeaturesXCMS <- function(xs, feat, rtalign, loadRawData, groupArgs, retcorArgs, verbose)
{
    if (length(feat) == 0)
        return(featureGroupsXCMS(features = feat))
    
    hash <- makeHash(feat, rtalign, loadRawData, groupArgs, retcorArgs)
    cachefg <- loadCacheData("featureGroupsXCMS", hash)
    if (!is.null(cachefg))
        return(cachefg)
    
    if (verbose)
        cat("Grouping features with XCMS...\n===========\n")
    
    if (verbose)
        xs <- do.call(xcms::group, c(list(xs), groupArgs))
    else
        suppressMessages(invisible(utils::capture.output(xs <- do.call(xcms::group, c(list(xs), groupArgs)))))
    
    if (!loadRawData && rtalign)
    {
        if (verbose)
            cat("Skipping RT alignment: no raw data\n")
    }
    else if (rtalign)
    {
        if (verbose)
        {
            xs <- do.call(xcms::retcor, c(list(xs), retcorArgs))
            xs <- do.call(xcms::group, c(list(xs), groupArgs))
        }
        else
        {
            suppressMessages(invisible(utils::capture.output(xs <- do.call(xcms::retcor, c(list(xs), retcorArgs)))))
            suppressMessages(invisible(utils::capture.output(xs <- do.call(xcms::group, c(list(xs), groupArgs)))))
        }
    }
    
    ret <- importFeatureGroupsXCMSFromFeat(xs, analysisInfo(feat), feat)
    saveCacheData("featureGroupsXCMS", ret, hash)
    
    if (verbose)
        cat("\n===========\nDone!\n")
    return(ret)
}

getFeatIndicesFromXS <- function(xs)
{
    plist <- as.data.table(xcms::peaks(xs))
    plist[, subind := seq_len(.N), by = "sample"]

    xsftidx <- xcms::groupval(xs)
    ret <- as.data.table(t(xsftidx))

    for (grp in seq_along(ret))
    {
        idx <- plist[ret[, grp, with = FALSE][[1]], subind]
        set(ret, j = grp, value = idx)
    }

    ret[is.na(ret)] <- 0

    return(ret)
}

importFeatureGroupsXCMSFromFeat <- function(xs, analysisInfo, feat)
{
    xcginfo <- xcms::groups(xs)
    gNames <- makeFGroupName(seq_len(nrow(xcginfo)), xcginfo[, "rtmed"], xcginfo[, "mzmed"])
    gInfo <- data.frame(rts = xcginfo[, "rtmed"], mzs = xcginfo[, "mzmed"], row.names = gNames, stringsAsFactors = FALSE)

    groups <- data.table(t(xcms::groupval(xs, value = "maxo")))
    setnames(groups, gNames)
    groups[is.na(groups)] <- 0

    return(featureGroupsXCMS(xs = xs, groups = groups, groupInfo = gInfo, features = feat,
                             ftindex = setnames(getFeatIndicesFromXS(xs), gNames)))
}

#' Imports feature groups from XCMS (old interface)
#'
#' Imports grouped features from a legacy \code{\link{xcmsSet}} object from the \pkg{xcms} package.
#'
#' @template analysisInfo-arg
#' @param xs An \code{\link{xcmsSet}} object.
#' 
#' @inherit groupFeaturesXCMS references
#' @inherit importFeatureGroups return
#'
#' @seealso \code{\link{importFeaturesXCMS3}} and \code{\link{groupFeatures}}
#'
#' @export
importFeatureGroupsXCMS <- function(xs, analysisInfo)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(xs, "xcmsSet", add = ac)
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, c("mzXML", "mzML"), add = ac)
    checkmate::reportAssertions(ac)

    if (length(xcms::groups(xs)) == 0)
        stop("Provided xcmsSet does not contain any grouped features!")

    feat <- importFeaturesXCMS(xs, analysisInfo)
    return(importFeatureGroupsXCMSFromFeat(xs, analysisInfo, feat))
}

#' @rdname featureGroups-class
#' @export
setReplaceMethod("analysisInfo", "featureGroupsXCMS", function(obj, value)
{
    obj <- callNextMethod()
    xcms::sampclass(obj@xs) <- analysisInfo(obj)$group # sync
    return(obj)
})

#' @rdname featureGroups-class
#' @export
setMethod("delete", "featureGroupsXCMS", function(obj, ...)
{
    # UNDONE: update individual features somehow too?
    
    old <- obj
    obj <- callNextMethod()
    
    if (length(old) > length(obj))
    {
        keep <- names(old) %chin% names(obj)
        xcms::groups(obj@xs) <- xcms::groups(obj@xs)[keep, , drop = FALSE]
        xcms::groupidx(obj@xs) <- xcms::groupidx(obj@xs)[keep]
    }

    return(obj)
})
