#' @include features.R

#' @rdname suspect-screening
featuresBrukerTASQ <- setClass("featuresBrukerTASQ", contains = "features")

setMethod("initialize", "featuresBrukerTASQ",
          function(.Object, ...) callNextMethod(.Object, algorithm = "bruker_tasq", ...))


# internally used by TASQ feature groups
importFeaturesBrukerTASQ <- function(analysisInfo, TASQExportFile)
{
    cat("Importing features from TASQ...")

    selCols <- c("Row", "Data Set", "RT [min]", "m/z meas.", "Height", "Area", "FWHM [s]")
    tExport <- fread(TASQExportFile, select = selCols)
    setnames(tExport, selCols, c("ID", "analysis", "ret", "mz", "intensity", "area", "FWHM"))

    tExport <- tExport[!is.na(ret)] # skip empty results
    tExport[, ret := ret * 60] # min --> s

    tAnalyses <- unique(tExport$analysis)
    tAnalyses <- tAnalyses[tAnalyses %in% analysisInfo$analysis]
    analysisInfo <- analysisInfo[analysisInfo$analysis %in% tAnalyses, ]

    ret <- featuresBrukerTASQ(analysisInfo = analysisInfo)

    ret@features <- sapply(analysisInfo$analysis, function(ana)
    {
        ft <- tExport[analysis == ana]
        ft[, analysis := NULL]

        # estimate min/max ret from fwhm (https://en.wikipedia.org/wiki/Full_width_at_half_maximum)
        ft[, s := FWHM / 2.355]
        ft[, retmin := ret - 2*s]
        ft[, retmax := ret + 2*s]
        
        # dummy ranges
        ft[, mzmin := mz - 0.0025]
        ft[, mzmax := mz + 0.0025]

        return(ft[, -c("FWHM", "s")])
    }, simplify = FALSE)
    
    cat("Done!\n")

    return(ret)
}
