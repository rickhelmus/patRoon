#' @include feature_groups.R

#' @rdname suspect-screening
featureGroupsBrukerTASQ <- setClass("featureGroupsBrukerTASQ", contains = "featureGroups")

setMethod("initialize", "featureGroupsBrukerTASQ",
          function(.Object, ...) callNextMethod(.Object, algorithm = "bruker_tasq", ...))

#' Imports feature groups from Bruker TASQ
#'
#' Imports screening results from Bruker TASQ as feature groups.
#'
#' @templateVar algo Bruker TASQ
#' @templateVar generic importFeatureGroups
#' @templateVar algoParam brukertasq
#' @template algo_importer
#'
#' @details The feature groups across analyses are formed based on the name of suspects and their closeness in retention
#'   time. The latter is necessary because TASQ does not necessarily perform checks on retention times and may therefore
#'   assign a suspect to peaks with different retention times across analyses (or within a single analysis). Hence,
#'   suspects with equal names are hierarchically clustered on their retention times (using \pkg{\link{fastcluster}}) to
#'   form the feature groups. The cut-off value for this is specified by the \code{clusterRTWindow} argument. The input
#'   for this function is obtained by generating an Excel export of the 'global' results and subsequently converting the
#'   file to \file{.csv} format.
#'
#' @template analysisInfo-arg
#'
#' @param path The file path to an Excel export of the Global results table from TASQ, converted to \file{.csv} format.
#' @param clusterRTWindow This retention time window (in seconds) is used to group hits across analyses together. See
#'   also the details section.
#'
#' @return A new \code{featureGroups} object containing converted screening results from Bruker TASQ.
#'
#' @note This function uses estimated min/max values for retention times and dummy min/max \emph{m/z} values for
#'   conversion to features, since this information is not (readily) available. Hence, when plotting, for instance,
#'   extracted ion chromatograms (with \code{\link{plotChroms}}) the integrated chromatographic peak range shown is
#'   incorrect.
#'
#'   This function may use suspect names to base file names used for reporting, logging etc. Therefore, it is important
#'   that these are file-compatible names.
#'
#' @references \addCitations{fastcluster}{1}
#'
#' @export
importFeatureGroupsBrukerTASQ <- function(path, analysisInfo, clusterRTWindow = 12)
{
    ac <- checkmate::makeAssertCollection()
    assertCSVFile(path, TASQImportCols(), add = ac)
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, add = ac)
    checkmate::assertNumber(clusterRTWindow, finite = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    fts <- importFeaturesBrukerTASQ(analysisInfo, path)
    fTable <- featureTable(fts)
    analysisInfo <- fts@analysisInfo # may have been updated
    
    tExport <- loadTASQFile(path, analysisInfo)

    # If no retention times were specified in TASQ for a suspect then we cannot
    # assume that suspect hits across analyses are from the same peak. Hence, to
    # group the suspects across analyses, we need to
    # - make sure unique names are assigned to each suspect hit at different retention times
    # - the duplicates should be grouped by similar retention time
    # --> perform HCA for each suspect to group them by close retention time and then re-assign names based on cluster number
    
    tExport[, cl := {
        distm <- dist(ret)
        hc <- fastcluster::hclust(distm)
        cutree(hc, h = clusterRTWindow)
    }, by = "group"]
    
    tExport[cl > 1, group := paste0(group, "-", cl)]

    # calculate group averages
    gInfoDT <- tExport[, .(rts = mean(ret), mzs = mean(mz)), by = "group"]
    
    gInfo <- as.data.frame(gInfoDT[, c("rts", "mzs")])
    rownames(gInfo) <- gInfoDT$group

    groupNames <- unique(tExport$group)
    groups <- data.table(matrix(0, nrow = nrow(analysisInfo), ncol = length(groupNames)))
    setnames(groups, groupNames)

    ftindex <- data.table(matrix(0L, nrow = nrow(analysisInfo), ncol = length(groupNames)))
    setnames(ftindex, groupNames)

    for (grp in groupNames)
    {
        g <- tExport[group == grp, c("ID", "analysis", "intensity")]
        anai <- match(g$analysis, analysisInfo$analysis) # align analysis order

        set(groups, anai, grp, as.numeric(g$intensity))

        finds <- sapply(seq_len(nrow(g)), function(i) fTable[[g$analysis[i]]][ID == g$ID[i], which = TRUE])
        set(ftindex, anai, grp, finds)
    }

    ret <- featureGroupsBrukerTASQ(groups = groups, groupInfo = gInfo, ftindex = ftindex, features = fts)

    return(ret)
}
