#' @include feature_groups.R

#' @rdname suspect-screening
featureGroupsBrukerTASQ <- setClass("featureGroupsBrukerTASQ", contains = "featureGroups")

setMethod("initialize", "featureGroupsBrukerTASQ",
          function(.Object, ...) callNextMethod(.Object, algorithm = "bruker_tasq", ...))

#' @details \code{importFeatureGroupsBrukerTASQ} will convert screening results
#'   from Bruker TASQ to a \code{\link{featureGroups}} object. The feature
#'   groups across analyses are formed based on the name of suspects and their
#'   closeness in retention time. The latter is necessary because TASQ does not
#'   necessarily perform checks on retention times and may therefore assign a
#'   suspect to peaks with different retention times across analyses (or within
#'   a single analysis). Hence, suspects with equal names are hierarchically
#'   clustered on their retention times (using \pkg{\link{fastcluster}}) to form
#'   the feature groups. The cut-off value for this is specified by the
#'   \code{clusterRTWindow} argument. The input for this function is obtained by
#'   generating an Excel export of the 'global' results and subsequently
#'   converting the file to \file{.csv} format.
#'
#' @param path The file path to an Excel export of the Global results table from
#'   TASQ, converted to \file{.csv} format.
#' @param analysisInfo A table with \link[=analysis-information]{analysis
#'   information}.
#' @param clusterRTWindow This retention time window (in seconds) is used to
#'   group hits across analyses together. See also the details section.
#'
#' @return \code{importFeatureGroupsBrukerTASQ} returns a new
#'   \code{featureGroups} object containing converted screening results from
#'   Bruker TASQ.
#'
#' @note \code{importFeatureGroupsBrukerTASQ} will use estimated min/max values
#'   for retention times and dummy min/max \emph{m/z} values for conversion to
#'   features, since this information is not (readily) available. Hence, when
#'   plotting, for instance, extracted ion chromatograms (with
#'   \code{\link{plotChroms}}) the integrated chromatographic peak range shown
#'   is incorrect.
#'
#'
#' @references \addCitations{fastcluster}{1}
#'
#' @rdname suspect-screening
#' @export
importFeatureGroupsBrukerTASQ <- function(path, analysisInfo, clusterRTWindow = 12)
{
    selCols <- c("Row", "Data Set", "Analyte", "RT [min]", "m/z meas.", "Height")

    ac <- checkmate::makeAssertCollection()
    assertCSVFile(path, selCols, add = ac)
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, add = ac)
    checkmate::assertNumber(clusterRTWindow, finite = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    fts <- importFeaturesBrukerTASQ(analysisInfo, path)
    fTable <- featureTable(fts)
    analysisInfo <- fts@analysisInfo # may have been updated
    
    tExport <- fread(path, select = selCols)
    tExport <- tExport[!is.na(`RT [min]`) & `Data Set` %in% analysisInfo$analysis] # skip empty/other results
    setnames(tExport, selCols, c("ID", "analysis", "group", "rts", "mzs", "intensity"))

    tExport[, rts := rts * 60] # min --> s
    
    # If no retention times were specified in TASQ for a suspect then we cannot
    # assume that suspect hits across analyses are from the same peak. Hence, to
    # group the suspects across analyses, we need to
    # - make sure unique names are assigned to each suspect hit at different retention times
    # - the duplicates should be grouped by similar retention time
    # --> perform HCA for each suspect to group them by close retention time and then re-assign names based on cluster number
    
    tExport[, cl := {
        distm <- dist(rts)
        hc <- fastcluster::hclust(distm)
        cutree(hc, h = clusterRTWindow)
    }, by = "group"]
    
    tExport[cl > 1, group := paste0(group, "-", cl)]

    # calculate group averages
    gInfoDT <- tExport[, .(rts = mean(rts), mzs = mean(mzs)), by = "group"]
    
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

    ret <- featureGroupsBrukerTASQ(groups = groups, analysisInfo = analysisInfo, groupInfo = gInfo, ftindex = ftindex,
                                   features = fts)

    return(ret)
}
