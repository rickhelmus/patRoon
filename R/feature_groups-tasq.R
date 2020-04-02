#' @include feature_groups.R

#' @rdname suspect-screening
featureGroupsBrukerTASQ <- setClass("featureGroupsBrukerTASQ", contains = "featureGroups")

setMethod("initialize", "featureGroupsBrukerTASQ",
          function(.Object, ...) callNextMethod(.Object, algorithm = "bruker_tasq", ...))

# UNDONE: does this work with analytes that are matched more than once in an analysis?

#' @details \code{importFeatureGroupsBrukerTASQ} will convert screening results
#'   from Bruker TASQ to a \code{\link{featureGroups}} object. The groups across
#'   analyses are formed by the name of suspects. However, for suspects that
#'   were found >1 in the same analysis ambiguity exists as the same name occurs
#'   multiple times for these analyses. For this situation, grouping is
#'   performed by clustering on closeness of retention times (using
#'   \pkg{\link{fastcluster}}). The cut-off value for this is specified by the
#'   \code{clusterRTWindow} argument. The input for this function is obtained by
#'   generating an Excel export of the 'global' results and subsequently
#'   converting the file to \file{.csv} format. Similar to
#'   \code{groupFeaturesScreening}, this method will return an object that is
#'   suitable for any further workflow processes.
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
#'   \code{\link{plotEIC}}) the integrated chromatographic peak range shown is
#'   incorrect.
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
    
    tExport <- fread(path, select = selCols)
    tExport <- tExport[!is.na(`RT [min]`) & `Data Set` %in% analysisInfo$analysis] # skip empty/other results
    setnames(tExport, selCols, c("ID", "analysis", "group", "rts", "mzs", "intensity"))

    tExport[, rts := rts * 60] # min --> s
    
    # if TASQ detects >1 of the same suspect in a single analysis it will return
    # each result as a row. Since they all have the same name we cannot easily
    # group therm across analyses. Thus,
    # - the duplicate hits must have unique names
    # - the duplicates should be grouped by similar retention time
    # --> perform HCA for each suspect to group them by close retention time and then re-assign names based on cluster number
    
    tExport[, cl := {
        if (anyDuplicated(.SD))
        {
            distm <- dist(rts)
            hc <- fastcluster::hclust(distm)
            cutree(hc, h = 12) # UNDONE: make h configurable
        }
        else
            1
    }, by = "group", .SDcols = c("analysis", "group")]
    
    tExport[cl > 1, group := paste0(group, "-", cl)]

    gInfoDT <- tExport[!duplicated(group), c("group", "rts", "mzs")]
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
