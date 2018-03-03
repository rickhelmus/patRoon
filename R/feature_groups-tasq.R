#' @include feature_groups.R

featureGroupsBrukerTASQ <- setClass("featureGroupsBrukerTASQ", contains = "featureGroups")

# UNDONE: does this work with analytes that are matched more than once in an analysis?

#' @details \code{importFeatureGroupsBrukerTASQ} will convert screening results
#'   from Bruker TASQ to a \code{\link{featureGroups}} object. Groups are made
#'   based on target analyte names and individual results per analysis are used
#'   to generate the features. The input for this function is obtained by
#'   generating an Excel export of the 'global' results and subsequently
#'   converting the file to \file{.csv} format. Similar to
#'   \code{groupFeaturesScreening}, this method will return an object that is
#'   suitable for any further workflow processes.
#'
#' @param path The file path to an Excel export of the Global results table from
#'   TASQ, converted to \file{.csv} format.
#' @param analysisInfo A table with \link[=analysis-information]{analysis
#'   information}.
#'
#' @return \code{importFeatureGroupsBrukerTASQ} returns a new
#'   \code{featureGroups} object containing converted screening results from
#'   Bruker TASQ.
#'
#' @note \code{importFeatureGroupsBrukerTASQ} will use dummy values for
#'   retention time and \emph{m/z} values for conversion to features, since this
#'   information is not available. Hence, when plotting, for instance, extracted
#'   ion chromatograms (with \code{\link{plotEIC}}) the integrated
#'   chromatographic peak range shown is incorrect.
#'
#' @rdname target-screening
#' @export
importFeatureGroupsBrukerTASQ <- function(path, analysisInfo)
{
    selCols <- c("Row", "Sample", "Analyte", "RT [min]", "m/z meas.", "Height")
    tExport <- fread(path, select = selCols)
    tExport <- tExport[!is.na(`RT [min]`) & Sample %in% analysisInfo$analysis] # skip empty/other results
    setnames(tExport, selCols, c("ID", "analysis", "group", "rts", "mzs", "intensity"))

    fts <- importFeaturesBrukerTASQ(analysisInfo, path)
    fTable <- featureTable(fts)

    gInfoDT <- tExport[!duplicated(group), c("group", "rts", "mzs")]
    gInfoDT[, rts := rts * 60] # min --> s
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
