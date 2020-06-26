#' @include main.R
#' @include feature_groups.R
NULL

#' @rdname suspect-screening
featureGroupsScreening <- setClass("featureGroupsScreening", contains = "featureGroups")

setMethod("initialize", "featureGroupsScreening",
          function(.Object, ...) callNextMethod(.Object, algorithm = "screening", ...))

#' @details \code{groupFeaturesScreening} uses results from \code{screenSuspects}
#'   to transform an existing \code{\link{featureGroups}} object by (1) renaming
#'   any matched feature groups by the respective name of the suspect and (2)
#'   filtering out any feature groups that were not matched. A common workflow
#'   is to first obtain and group features (with \emph{e.g.}
#'   \code{\link{findFeatures}} and \code{\link{groupFeatures}}), screen them
#'   with \code{screenSuspects}, convert the \code{featureGroups} object that was
#'   used for screening with this method function and continue any further
#'   workflow steps such as compound identification as with 'regular'
#'   \code{featureGroups}.
#'
#' @param fGroups The \code{\link{featureGroups}} object that should be
#'   transformed (and was used to obtain the screening results).
#' @param scr The screening results table returned by \code{screenSuspects}.
#'
#' @return \code{groupFeaturesScreening} returns a modified \code{featureGroups}
#'   object in which those feature groups that were not matched by any suspects
#'   are removed and others renamed by the respective suspect name. In case of
#'   duplicate suspect results, feature group names are made unique with
#'   \code{\link{make.unique}}.
#'
#' @note Please note that \code{groupFeaturesScreening} method can only
#'   transform the \code{featureGroups} object that was used to obtain the given
#'   screening results.
#'
#' @rdname suspect-screening
#' @aliases groupFeaturesScreening
#' @export
setMethod("groupFeaturesScreening", "featureGroups", function(fGroups, scr)
{
    assertScreeningResults(scr, fromFGroups = TRUE)

    cat("Converting screening results to feature groups... ")

    if (any(is.na(scr$group)))
    {
        cat("\nRemoving empty screening results\n")
        scr <- scr[!is.na(scr$group), ]
    }

    feat <- getFeatures(fGroups)
    anaInfo <- analysisInfo(fGroups)
    gNames <- make.unique(scr$name)
    fgNames <- colnames(groupTable(fGroups))

    gInfo <- data.frame(rts = scr$exp_rt, mzs = scr$exp_mz, row.names = gNames)

    rmcols <- c("name", "rt", "mz", "group", "d_rt", "d_mz", "exp_rt", "exp_mz")
    groups <- transpose(as.data.table(scr)[, -rmcols, with = FALSE])
    setnames(groups, gNames)

    ftind <- copy(groupFeatIndex(fGroups))
    ftind <- ftind[, scr$group, with = FALSE] # re-order and isolate screened groups
    setnames(ftind,  gNames)

    ret <- featureGroupsScreening(groups = groups, analysisInfo = anaInfo, groupInfo = gInfo, features = feat, ftindex = ftind)

    cat("Done!\n")
    return(ret)
})
