#' @include features.R
NULL

SAFDMPFinishHandler <- function(cmd)
{
    # UNDONE
    fread(file.path(cmd$outPath, cmd$fileName))
}

#' @rdname features-class
#' @export
featuresSAFD <- setClass("featuresSAFD", contains = "features")

setMethod("initialize", "featuresSAFD",
          function(.Object, ...) callNextMethod(.Object, algorithm = "safd", ...))

makeSAFDCommand <- function(inPath, outPath, fileName, mzRange, maxNumbIter, maxTPeakW, resolution,
                            minMSW, RThreshold, minInt, sigIncThreshold, S2N, minPeakWS)
{
    # UNDONE: check if julia exists? allow to configure path?
    return(list(command = "julia", args = c(inPath, outPath, fileName, mzRange[1], mzRange[2],
                                            maxNumbIter, maxTPeakW, resolution,
                                            minMSW, RThreshold, minInt, sigIncThreshold, S2N,
                                            minPeakWS),
                outPath = outPath, fileName = fileName))
}

importSAFDPeakList <- function(peaklist)
{
    # peaklist is a single number (zero) when no results
    if (length(peaklist) == 1 || nrow(peaklist) == 0)
        return(data.table(ID = character(), ret = numeric(), mz = numeric(), intensity = numeric(),
                          area = numeric(), retmin = numeric(), retmax = numeric(), mzmin = numeric(),
                          mzmax = numeric()))

    ft <- as.data.table(peaklist)

    setnames(ft, c("m/z", "max_int", "sum_int", "RT", "minRT", "maxRT", "peak_ID"),
             c("mz", "intensity", "area", "ret", "retmin", "retmax", "ID"))

    # Estimate mzrange from variance
    s <- sqrt(ft$`var_m/z`)
    ft[, mzmin := mz - 2*s]
    ft[, mzmax := mz + 2*s]

    return(ft[, c("ID", "ret", "mz", "intensity", "area", "retmin", "retmax", "mzmin", "mzmax")])
}


#' @details \code{findFeaturesSAFD} uses the
#'   \code{\link[enviPick]{enviPickwrap}}. function from the \pkg{enviPick} R
#'   package to extract features.
#'
#' @note \code{findFeaturesSAFD} Requires analysis files to be in the
#'   \code{mzXML} format.
#'
#' @rdname feature-finding
#' @export
findFeaturesSAFD <- function(analysisInfo, profPath, mzRange = c(0, 400), 
                             maxNumbIter = 1000, maxTPeakW = 300, resolution = 30000,
                             minMSW = 0.02, RThreshold = 0.75, minInt = 2000,
                             sigIncThreshold = 5, S2N = 2, minPeakWS = 3, verbose = TRUE)
{
    # UNDONE: more assertions
    # UNDONE: docs
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, add = ac)
    checkmate::assertCharacter(profPath, min.chars = 1, min.len = 1, null.ok = TRUE, add = ac)
    assertCanCreateDirs(profPath, add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    anaCount <- nrow(analysisInfo)
    profPath <- rep(profPath, length.out = anaCount)
    
    params <- list(mzRange, maxNumbIter, maxTPeakW, resolution, minMSW, RThreshold, minInt,
                   sigIncThreshold, S2N, minPeakWS)
    baseHash <- makeHash(params)
    
    if (verbose)
        printf("Finding features with SAFD for %d analyses ...\n", anaCount)

    outp <- tempfile("safd")
    mkdirp(outp)
    
    cmdQueue <- Map(analysisInfo$analysis, profPath, f = function(ana, path)
    {
        fpMZXML <- getMzXMLAnalysisPath(ana, path)
        fpNCDF <- getAnalysisPath(ana, path, "cdf")
        fp <- NULL
        
        if (file.exists(fpMZXML))
            fp <- fpMZXML
        else if (file.exists(fpNCDF))
            fp <- fpNCDF
        else
            stop(sprintf("Cannot find %s or %s\n", fpMZXML, fpNCDF))
        
        # UNDONE: really want to hash big files?
        hash <- makeHash(makeFileHash(fp), params)
        
        cmd <- do.call(makeSAFDCommand, c(list(dirname(fp), outp, basename(fp)), params))
        
        return(c(cmd, list(inputPath = fp, hash = hash, logFile = paste0(ana, ".txt"))))
    })
    
    fts <- executeMultiProcess(cmdQueue, patRoon:::SAFDMPFinishHandler,
                               cacheName = "featuresSAFD", logSubDir = "safd",
                               showProgress = verbose)
        
    if (verbose)
    {
        printf("Done!\n")
        printFeatStats(fts)
    }
    
    return(featuresSAFD(analysisInfo = analysisInfo, features = fts))
}
