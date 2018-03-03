#' @include main.R
NULL

#' Base features class
#'
#' Holds information for all features present within a set of analysis.
#'
#' This class provides a way to store intensity, retention times, \emph{m/z} and
#' other data for all features in a set of analyses. The class is \code{virtual}
#' and derived objects are created by 'feature finders' such as
#' \code{findFeaturesOpenMS}, \code{findFeaturesXCMS} and
#' \code{findFeaturesBruker}.
#'
#' @param obj,x,object \code{features} object to be accessed
#'
#' @seealso \code{\link{feature-finding}}
#'
#' @slot features List of features per analysis file. Use the
#'   \code{featureTable} method for access.
#' @slot analysisInfo Analysis group information. Use the \code{analysisInfo} method
#'   for access.
#' @export
features <- setClass("features",
                     slots = c(features = "list", analysisInfo = "data.frame"),
                     prototype = list(features = list(), analysisInfo = data.frame()),
                     contains = "VIRTUAL")

#' @describeIn features Obtain total number of features.
#' @export
setMethod("length", "features", function(x) sum(sapply(x@features, nrow)))

#' @describeIn features Shows summary information for this object.
#' @export
setMethod("show", "features", function(object)
{
    ftcounts <- sapply(object@features, nrow)
    printf("A features object ('%s')\n", class(object))
    printf("Total feature count: %d\n", sum(ftcounts))
    printf("Average feature count/analysis: %.0f\n", sum(ftcounts) / nrow(analysisInfo(object)))
    printf("Least features: %s\n", names(object)[which.min(ftcounts)])
    printf("Most features: %s\n", names(object)[which.max(ftcounts)])
    showAnaInfo(analysisInfo(object))
    showObjectSize(object)
})

#' @describeIn features Get table with feature information
#'
#' @return \code{featureTable}: A \code{list} containing a
#'   \code{\link{data.table}} for each analysis with feature data
#'
#' @export
setMethod("featureTable", "features", function(obj) obj@features)

#' @describeIn features Get analysis information
#' @return \code{analysisInfo}: A \code{data.frame} containing a column with
#'   analysis name (\code{analysis}), its path (\code{path}), and other columns
#'   such as replicate group name (\code{group}) and blank reference
#'   (\code{ref}).
#' @export
setMethod("analysisInfo", "features", function(obj) obj@analysisInfo)

#' @rdname target-screening
#' @export
setMethod("screenTargets", "features", function(obj, targets, rtWindow, mzWindow)
{
    fTable <- featureTable(obj)
    anaInfo <- analysisInfo(obj)

    retlist <- lapply(seq_len(nrow(targets)), function (ti)
    {
        rbindlist(lapply(names(fTable), function(ana)
        {
            fts <- fTable[[ana]][abs(ret - targets$rt[ti]) <= rtWindow & abs(mz - targets$mz[ti]) <= mzWindow, ]

            if (nrow(fts) == 0) # no results? --> add NA result
                return(data.table(name = targets$name[ti], rt = targets$rt[ti], mz = targets$mz[ti], analysis = ana,
                                  feature = NA, d_rt = NA, d_mz = NA, intensity = NA, area = NA))

            return(rbindlist(lapply(seq_len(nrow(fts)), function(i)
            {
                data.table(name = targets$name[ti], rt = targets$rt[ti], mz = targets$mz[ti], analysis = ana,
                           feature = fts[["ID"]][i], d_rt = fts[["ret"]][i] - targets$rt[ti],
                           d_mz = fts[["mz"]][i] - targets$mz[ti], intensity = fts[["intensity"]][i],
                           area = if (is.null(fts[["area"]][i])) 0 else fts[["area"]][i])
            })))

        }))
    })

    return(rbindlist(retlist, fill = TRUE))
})

#' @templateVar func findFeatures
#' @templateVar what find features
#' @templateVar ex1 findFeaturesOpenMS
#' @templateVar ex2 findFeaturesBruker
#' @templateVar algos bruker,openms,xcms,envipick
#' @template generic-algo
#'
#' @rdname feature-finding
#' @aliases findFeatures
#' @export
findFeatures <- function(analysisInfo, algorithm, ...)
{
    f <- switch(algorithm,
                bruker = findFeaturesBruker,
                openms = findFeaturesOpenMS,
                xcms = findFeaturesXCMS,
                envipick = findFeaturesEnviPick,
                stop("Invalid algorithm! Should be: bruker, openms, xcms or envipick"))

    f(analysisInfo, ...)
}
