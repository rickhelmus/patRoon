#' @include features.R
#' @include main.R
NULL

#' @rdname features-class
#' @export
featuresQAlgorithms <- setClass("featuresQAlgorithms", contains = "features")

setMethod("initialize", "featuresQAlgorithms",
          function(.Object, ...) callNextMethod(.Object, algorithm = "qAlgorithms", ...))


#' Find features using qAlgorithms
#'
#' uses the \href{https://github.com/odea-project/qAlgorithms}{qAlgorithms} software to find features in LC-MS data.
#'
#' @templateVar algo qalgorithms
#' @templateVar do automatically find features
#' @templateVar generic findFeatures
#' @templateVar algoParam qalgorithms
#' @template algo_generator
#' @param ppm The mass accuracy in parts per million (ppm) for the qPeaks when using centroid data. Default is 5 ppm.
#'
#' @details This functionality has been tested with qAlgorithms version >= 0.1.1. Please make sure it is installed and
#'   configured, e.g. by configuring the path of the binaries with the \code{patRoon.path.qAlgorithms} option.
#'
#'   The file format of analyses must be \file{mzML}. Both profile and centroid data are supported. Yet, the profile 
#'   data is recommended for better results.
#'
#' @inheritParams findFeatures
#'
#' @templateVar what \code{findFeaturesQAlgorithms}
#'
#' @inherit findFeatures return
#'
#' @export
findFeaturesQAlgorithms <- function(analysisInfo, ppm = 5, verbose = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, "mzML", verifyCentroided = FALSE, add = ac)
    
    if (verbose)
        printf("Finding features with qAlgorithms for %d analyses ...\n", nrow(analysisInfo))
    
    outDir <- paste0(getwd(), "/log/qAlgorithms")
    
    if (!dir.exists(outDir))
        dir.create(outDir, recursive = TRUE)
    
    qAlgoExe <- getExtDepPath("qalgorithms", "qAlgorithms")

    cmdQueue <- lapply(seq_len(nrow(analysisInfo)), function(anai)
    {
        dfile <- getMzMLAnalysisPath(analysisInfo$analysis[anai], analysisInfo$path[anai])
        hash <- makeHash(makeFileHash(dfile), list("qAlgorithms"))
        return(list(hash = hash, dataFile = dfile))
    })
    names(cmdQueue) <- analysisInfo$analysis

    fList <- list()
    if (length(cmdQueue) > 0)
    {
        fList <- lapply(cmdQueue, function(cmd, ...)
        {
            pksOut <- loadCacheData("featuresQAlgorithms", cmd$hash)
            if (!is.null(pksOut))
                return(pksOut)
            else
            {
                fl <- cmd$dataFile
                
                if (verbose)
                    printf("Finding features for '%s' ...\n", fl)
                
                # deletes previous results for analysis file, needed?
                res <- list.files(outDir, full.names = TRUE)
                res <- res[grep(tools::file_path_sans_ext(basename(fl)), res)]
                res <- lapply(res, file.remove)
                
                # ppm <- sprintf("-ppm %f", ppm)
                
                # how to silence qBinning log output and warning and how to use -ppm argument
                executeCommand(qAlgoExe, args = c("-i", fl, "-o", outDir, "-e", "-ppm", ppm), stdout = verbose)
                
                # polarity switching generates two result files
                res <- list.files(outDir, full.names = TRUE)
                res <- res[grep(tools::file_path_sans_ext(basename(fl)), res)]
                
                if (!any(vapply(res, file.exists, FALSE)))
                {
                    warning("No results file found for ", fl)
                    return()
                }
                else
                {
                    # handles polarity switching by adding adduct and mass columns and amending ID string
                    # check with Rick future plans for handling polarity switching
                    # also, shall I create an importFeaturesqAlgorithms function?
                    pksList <- list()
                    for (i in seq_along(res))
                    {
                        pks <- data.table::fread(res[i])
                        
                        if (grepl("positive", res[i]))
                        {
                            pks$adduct <- "[M+H]+"
                            pks$mass <- pks$mz - 1.007276
                            pks$ID <- paste0("pos_", pks$ID)
                        }
                        
                        if (grepl("negative", res[i]))
                        {
                            pks$adduct <- "[M-H]-"
                            pks$mass <- pks$mz + 1.007276
                            pks$ID <- paste0("neg_", pks$ID)
                        }
                        
                        pks$mzmin <- pks$mz - pks$mzUncertainty
                        pks$mzmax <- pks$mz + pks$mzUncertainty
                        pks$ret <- pks$retentionTime
                        pks$retentionTime <- NULL
                        pks$retmin <- pks$lowestRetentionTime
                        pks$lowestRetentionTime <- NULL
                        pks$retmax <- pks$highestRetentionTime
                        pks$highestRetentionTime <- NULL
                        pks$retUncertainty <- pks$retentionTimeUncertainty
                        pks$retentionTimeUncertainty <- NULL
                        pks$binID <- NULL
                        pks$intensity <- pks$height
                        pks$height <- NULL
                        pks$intensityUncertainty <- pks$heightUncertainty
                        pks$heightUncertainty <- NULL
                        
                        data.table::setcolorder(pks, c("ID", "mz", "ret", "mzmin", "mzmax", "retmin", "retmax", 
                                                       "intensity", "area", "mzUncertainty", "retUncertainty", 
                                                       "areaUncertainty"))
                        
                        pksList[[i]] <- pks
                    }
                    
                    pksOut <- data.table::rbindlist(pksList)
                    saveCacheData("featuresQAlgorithms", pksOut, cmd$hash)
                    pksOut
                }
            }
        })
    }
    
    if (verbose)
    {
        printf("Done!\n")
        printFeatStats(fList)
    }

    return(featuresQAlgorithms(analysisInfo = analysisInfo, features = fList))
}
