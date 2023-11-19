#' @include features.R
NULL

# NOTE: first column is ID and unnamed (V1)
TASQImportCols <- function() c("V1", "Data Set", "Analyte Name", "RT [min]", "m/z meas.", "Height of PI", "Area of PI",
                               "FWHM [s]")

loadTASQFile <- function(path, analysisInfo)
{
    # NOTE: first line is empty
    
    tExport <- fread(path, select = TASQImportCols())
    setnames(tExport, c("ID", "analysis", "group", "ret", "mz", "intensity", "area", "FWHM"))
    
    tExport <- tExport[!is.na(ret) & analysis %in% analysisInfo$analysis] # skip empty/other results
    tExport[, ret := ret * 60] # min --> s
    
    return(tExport)
}

#' @rdname suspect-screening
featuresBrukerTASQ <- setClass("featuresBrukerTASQ", contains = "features")

setMethod("initialize", "featuresBrukerTASQ",
          function(.Object, ...) callNextMethod(.Object, algorithm = "bruker_tasq", ...))


# internally used by TASQ feature groups
importFeaturesBrukerTASQ <- function(analysisInfo, TASQExportFile)
{
    cat("Importing features from TASQ...")

    tExport <- loadTASQFile(TASQExportFile, analysisInfo)
    tAnalyses <- unique(tExport$analysis)
    tAnalyses <- tAnalyses[tAnalyses %in% analysisInfo$analysis]
    analysisInfo <- analysisInfo[analysis %in% tAnalyses]

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
