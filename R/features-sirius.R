#' @include features.R
NULL

SIRFeatMPFinishHandler <- function(cmd)
{
    pattern <- paste0("^[[:digit:]]+_", cmd$analysis, "_[[:digit:]]+$")
    resDirs <- list.files(cmd$outPath, pattern = pattern, full.names = TRUE)
    ret <- rbindlist(lapply(resDirs, function(dir)
    {
        # assumptions: only one analysis, first isotope is our target monoisotopic mass
        json <- jsonlite::fromJSON(file.path(dir, "lcms.json.gz"), FALSE)
        iso <- json[["traceSets"]][[1]][["ionTrace"]][["isotopes"]][[1]]
        # NOTE: SIRIUS retention times are in milliseconds
        eic <- data.table(time = unlist(json[["traceSets"]][[1]][["retentionTimes"]]) / 1000,
                          intensity = unlist(iso[["intensities"]]))
        mzs <- unlist(iso[["masses"]])
        area <- json[["abundance"]][[1]]
        ret <- eic$time[which.max(eic$intensity)] # UNDONE: verify
        return(data.table(ret = ret, mz = mean(mzs), mzmin = min(mzs), mzmax = max(mzs),
                          retmin = min(eic$time), retmax = max(eic$time), area = area,
                          intensity = max(eic$intensity)))
    }))
    
    ret[, ID := seq_len(nrow(..ret))]
    setcolorder(ret, "ID")
    return(ret[])
}

SIRFeatMPPrepareHandler <- function(cmd)
{
    command <- patRoon:::getCommandWithOptPath(patRoon:::getSiriusBin(), "SIRIUS")
    outPath <- tempfile("sirius_out")
    args <- c("-i", cmd$dataFile, "-o", outPath, "lcms-align")
    return(utils::modifyList(cmd, list(command = command, args = args, outPath = outPath)))
}

#' @rdname features-class
#' @export
featuresSIRIUS <- setClass("featuresSIRIUS", contains = "features")

setMethod("initialize", "featuresSIRIUS",
          function(.Object, ...) callNextMethod(.Object, algorithm = "sirius", ...))

#' @details \code{findFeaturesSIRIUS} uses the
#'   \code{\link[enviPick]{enviPickwrap}}. function from the \pkg{enviPick} R
#'   package to extract features.
#'
#' @note \code{findFeaturesSIRIUS} Requires analysis files to be in the
#'   \code{mzXML} format.
#'
#' @rdname feature-finding
#' @export
findFeaturesSIRIUS <- function(analysisInfo, verbose = TRUE)
{
    # UNDONE: docs
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, "mzML", add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    anaCount <- nrow(analysisInfo)
    
    if (verbose)
        printf("Finding features with SIRIUS for %d analyses ...\n", anaCount)

    cmdQueue <- Map(analysisInfo$analysis, analysisInfo$path, f = function(ana, path)
    {
        dfile <- getMzMLAnalysisPath(ana, path)
        hash <- makeHash(makeFileHash(dfile))
        logf <- paste0(ana, ".txt")
        return(list(hash = hash, dataFile = dfile, analysis = ana, logFile = logf))
    })
    
    fList <- list()
    if (length(cmdQueue) > 0)
    {
        fList <- executeMultiProcess(cmdQueue, patRoon:::SIRFeatMPFinishHandler,
                                     prepareHandler = patRoon:::SIRFeatMPPrepareHandler,
                                     showProgress = verbose, logSubDir = "sirius_features",
                                     cacheName = "featuresSIRIUS")
    }
    
    if (verbose)
    {
        printf("Done!\n")
        printFeatStats(fList)
    }
    
    return(featuresSIRIUS(analysisInfo = analysisInfo, features = fList))
}
