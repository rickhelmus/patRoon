#' @include main.R
#' @include feature_groups.R
NULL

#' @rdname featureGroups-class
#' @export
featureGroupsXCMS3 <- setClass("featureGroupsXCMS3", slots = c(xdata = "XCMSnExp"), contains = "featureGroups")

setMethod("initialize", "featureGroupsXCMS3",
          function(.Object, ...) callNextMethod(.Object, algorithm = "xcms3", ...))


#' @details \code{groupFeaturesXCMS} uses the \pkg{xcms} package for grouping of
#'   features. Grouping of features and alignment of their retention times are
#'   performed with the \code{\link[xcms:group-methods]{xcms::group}} and
#'   \code{\link[xcms:retcor-methods]{xcms::retcor}} functions, respectively.
#'   Both functions have an extensive list of parameters to modify their
#'   behaviour and may therefore be used to potentially optimize results.
#'
#' @param exportedData Set to \code{TRUE} if analyses were exported as
#'   \code{mzXML} or \code{mzML} files.
#' @param groupArgs,retcorArgs named \code{character vector} that can contain
#'   extra parameters to be used by \code{\link[xcms:group-methods]{xcms::group}} and
#'   \code{\link[xcms:retcor-methods]{xcms::retcor}}, respectively.
#'
#' @references \addCitations{xcms}{1} \cr\cr
#'   \addCitations{xcms}{2} \cr\cr
#'   \addCitations{xcms}{3}
#'
#' @rdname feature-grouping
#' @export
groupFeaturesXCMS3 <- function(feat, rtalign = TRUE,
                               groupParam = xcms::PeakDensityParam(sampleGroups = analysisInfo(feat)$group),
                               retAlignParam = xcms::ObiwarpParam(), verbose = TRUE)
{
    # UNDONE: aligning gives XCMS errors when multithreading

    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(feat, "features", add = ac)
    aapply(checkmate::assertFlag, . ~ rtalign + verbose, fixed = list(add = ac))
    assertS4(groupParam, add = ac)
    assertS4(retAlignParam, add = ac)
    checkmate::reportAssertions(ac)

    anaInfo <- analysisInfo(feat)

    if (length(feat) == 0)
        return(featureGroupsXCMS(analysisInfo = anaInfo, features = feat))

    hash <- makeHash(feat, rtalign, groupParam, retAlignParam)
    cachefg <- loadCacheData("featureGroupsXCMS3", hash)
    if (!is.null(cachefg))
        return(cachefg)

    if (verbose)
        cat("Grouping features with XCMS...\n===========\n")

    xdata <- getXCMSnExp(feat)

    if (rtalign)
    {
        # group prior to alignment when using the peaks group algorithm
        if (is(groupParam, getClass("PeakGroupsParam", where = "xcms")))
        {
            if (verbose)
                cat("Performing grouping prior to retention time alignment...\n")
            xdata <- verboseCall(xcms::groupChromPeaks, list(xdata, groupParam), verbose)
        }

        if (verbose)
            cat("Performing retention time alignment...\n")
        xdata <- verboseCall(xcms::adjustRtime, list(xdata, retAlignParam), verbose)
    }

    if (verbose)
        cat("Performing grouping...\n")
    xdata <- verboseCall(xcms::groupChromPeaks, list(xdata, groupParam), verbose)

    ret <- importFeatureGroupsXCMS3FromFeat(xdata, anaInfo, feat)
    saveCacheData("featureGroupsXCMS", ret, hash)

    if (verbose)
        cat("\n===========\nDone!\n")
    return(ret)
}

getFeatIndicesFromXCMSnExp <- function(xdata)
{
    plist <- as.data.table(xcms::chromPeaks(xdata))
    plist[, subind := seq_len(.N), by = "sample"]

    xdftidx <- xcms::featureValues(xdata)
    ret <- as.data.table(t(xdftidx))

    for (grp in seq_along(ret))
    {
        idx <- plist[ret[, grp, with = FALSE][[1]], subind]
        set(ret, j = grp, value = idx)
    }

    ret[is.na(ret)] <- 0

    return(ret)
}

importFeatureGroupsXCMS3FromFeat <- function(xdata, analysisInfo, feat)
{
    xcginfo <- xcms::featureDefinitions(xdata)
    gNames <- makeFGroupName(seq_len(nrow(xcginfo)), xcginfo[, "rtmed"], xcginfo[, "mzmed"])
    gInfo <- data.frame(rts = xcginfo[, "rtmed"], mzs = xcginfo[, "mzmed"], row.names = gNames, stringsAsFactors = FALSE)

    groups <- data.table(t(xcms::featureValues(xdata, value = "maxo")))
    setnames(groups, gNames)
    groups[is.na(groups)] <- 0

    return(featureGroupsXCMS3(xdata = xdata, groups = groups, groupInfo = gInfo, analysisInfo = analysisInfo, features = feat,
                              ftindex = setnames(getFeatIndicesFromXCMSnExp(xdata), gNames)))
}

#' @details \code{importFeatureGroupsXCMS} converts grouped features from an
#'   \code{\link{xcmsSet}} object (from the \pkg{xcms} package).
#'
#' @param xs An \code{\link{xcmsSet}} object.
#' @param analysisInfo A \code{data.frame} with \link[=analysis-information]{analysis info}.
#'
#' @rdname feature-grouping
#' @export
importFeatureGroupsXCMS3 <- function(xdata, analysisInfo)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(xdata, "XCMSnExp", add = ac)
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, c("mzXML", "mzML"), add = ac)
    checkmate::reportAssertions(ac)

    if (length(xcms::hasFeatures(xdata)) == 0)
        stop("Provided XCMS data does not contain any grouped features!")

    feat <- importFeaturesXCMS(xdata, analysisInfo)
    return(importFeatureGroupsXCMSFromFeat(xdata, anaInfo, feat))
}

setMethod("removeGroups", "featureGroupsXCMS3", function(fGroups, indices)
{
    fGroups <- callNextMethod(fGroups, indices)

    # update XCMSnExp
    if (length(indices) > 0)
    {
        xcms::groups(fGroups@xs) <- xcms::groups(fGroups@xs)[-indices, , drop = FALSE]
        groupidx(fGroups@xs) <- groupidx(fGroups@xs)[-indices]
    }

    return(fGroups)
})
