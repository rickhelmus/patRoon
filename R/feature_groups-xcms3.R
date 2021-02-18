#' @include main.R
#' @include feature_groups.R
NULL

#' @rdname featureGroups-class
#' @export
featureGroupsXCMS3 <- setClass("featureGroupsXCMS3", slots = c(xdata = "ANY"), contains = "featureGroups")

setMethod("initialize", "featureGroupsXCMS3",
          function(.Object, ...) callNextMethod(.Object, algorithm = "xcms3", ...))


#' @details \code{groupFeaturesXCMS3} uses the new interface from the \pkg{xcms}
#'   package for grouping of features. Grouping of features and alignment of
#'   their retention times are performed with the
#'   \code{\link[xcms:groupChromPeaks]{xcms::groupChromPeaks}} and
#'   \code{\link[xcms:adjustRtime]{xcms::adjustRtime}} functions, respectively.
#'   Both of these functions support an extensive amount of parameters that
#'   modify their behaviour and may therefore require optimization.
#'
#' @param groupParam,retAlignParam parameter object that is directly passed to
#'   \code{\link[xcms:groupChromPeaks]{xcms::groupChromPeaks}} and
#'   \code{\link[xcms:adjustRtime]{xcms::adjustRtime}}, respectively.
#'
#' @references \addCitations{xcms}{1} \cr\cr \addCitations{xcms}{2} \cr\cr
#'   \addCitations{xcms}{3}
#'
#' @rdname feature-grouping
#' @export
groupFeaturesXCMS3 <- function(feat, rtalign = TRUE, exportedData = TRUE,
                               groupParam = xcms::PeakDensityParam(sampleGroups = analysisInfo(feat)$group),
                               retAlignParam = xcms::ObiwarpParam(), verbose = TRUE)
{
    # UNDONE: aligning gives XCMS errors when multithreading

    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(feat, "features", add = ac)
    aapply(checkmate::assertFlag, . ~ rtalign + exportedData + verbose, fixed = list(add = ac))
    assertS4(groupParam, add = ac)
    assertS4(retAlignParam, add = ac)
    checkmate::reportAssertions(ac)

    anaInfo <- analysisInfo(feat)

    if (length(feat) == 0)
        return(featureGroupsXCMS(analysisInfo = anaInfo, features = feat))

    hash <- makeHash(feat, rtalign, exportedData, groupParam, retAlignParam)
    cachefg <- loadCacheData("featureGroupsXCMS3", hash)
    if (!is.null(cachefg))
        return(cachefg)

    if (verbose)
        cat("Grouping features with XCMS...\n===========\n")

    xdata <- getXCMSnExp(feat, verbose = verbose, exportedData = exportedData)

    if (!exportedData && rtalign)
    {
        if (verbose)
            cat("Skipping RT alignment: no raw data\n")
    }
    else if (rtalign)
    {
        # group prior to alignment when using the peaks group algorithm
        if (isXCMSClass(retAlignParam, "PeakGroupsParam"))
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
    saveCacheData("featureGroupsXCMS3", ret, hash)

    if (verbose)
        cat("\n===========\nDone!\n")
    return(ret)
}

getFeatIndicesFromXCMSnExp <- function(xdata)
{
    plist <- as.data.table(xcms::chromPeaks(xdata))
    plist[, subind := seq_len(.N), by = "sample"]

    xdftidx <- xcms::featureValues(xdata, value = "index")
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

#' @details \code{importFeatureGroupsXCMS3} converts grouped features from an
#'   \code{\link{XCMSnExp}} object (from the \pkg{xcms} package).
#'
#' @param xdata An \code{\link{XCMSnExp}} object.
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

    feat <- importFeaturesXCMS3(xdata, analysisInfo)
    return(importFeatureGroupsXCMS3FromFeat(xdata, analysisInfo, feat))
}

setMethod("delete", "featureGroupsXCMS3", function(obj, ...)
{
    # UNDONE: update individual features somehow too?
    
    old <- obj
    obj <- callNextMethod()

    if (length(old) > length(obj))
        obj@xdata <- xcms::filterFeatureDefinitions(obj@xdata, names(old) %chin% names(obj))
    
    # UNDONE: enable when we can update features (or mark dirty so that getXCMSnExp knows when to update)
    # if (length(analyses(old)) > length(analyses(obj)))
    #     obj@xdata <- xcms::filterFile(obj@xdata, analyses(obj), keepFeatures = TRUE)
    
    return(obj)
})
