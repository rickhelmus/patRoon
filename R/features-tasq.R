#' @include features.R

#' @rdname target-screening
featuresBrukerTASQ <- setClass("featuresBrukerTASQ", contains = "features")

setMethod("initialize", "featuresBrukerTASQ",
          function(.Object, ...) callNextMethod(.Object, algorithm = "bruker_tasq", ...))


# internally used by TASQ feature groups
importFeaturesBrukerTASQ <- function(analysisInfo, TASQExportFile)
{
    cat("Importing features from TASQ...")

    selCols <- c("Row", "Sample", "RT [min] exp.", "m/z exp.", "Height", "Area")
    tExport <- fread(TASQExportFile, select = selCols)
    setnames(tExport, selCols, c("ID", "analysis", "ret", "mz", "intensity", "area"))

    tExport <- tExport[!is.na(ret)] # skip empty results
    tExport[, ret := ret * 60] # min --> s

    tAnalyses <- unique(tExport$analysis)
    tAnalyses <- tAnalyses[tAnalyses %in% analysisInfo$analysis]

    ret <- featuresBrukerTASQ(analysisInfo = analysisInfo)

    ret@features <- sapply(tAnalyses, function(ana)
    {
        ft <- tExport[analysis == ana]
        ft[, analysis := NULL]

        # dummy ranges
        ft[, retmin := ret - 3]
        ft[, retmax := ret + 3]
        ft[, mzmin := mz - 0.0025]
        ft[, mzmax := mz + 0.0025]

        return(ft)
    }, simplify = FALSE)

    cat("Done!\n")

    return(ret)
}
