#' @include features.R
NULL

SAFDMPFinishHandler <- function(cmd)
{
    results <- fread(file.path(cmd$outPath, paste0(cmd$fileName, "_report.csv")))
    setnames(results,
             c("Nr", "Rt", "MeasMass", "MinMass", "MaxMass", "Area", "Int",
               "FeatPurity", "MediRes"),
             c("ID", "ret", "mz", "mzmin", "mzmax", "area", "intensity",
               "purity", "mediRes"))
    results[, ret := ret * 60] # min --> sec
    # UNDONE
    results[, retmin := ret - (SecInPeak / 2)]
    results[, retmax := ret + (SecInPeak / 2)]
    results[, c("SecInPeak", "ScanNum", "ScanInPeak") := NULL]
    
    return(results[])
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
    return(list(command = "julia", args = c(system.file("misc", "runSAFD.jl", package = "patRoon"),
                                            inPath, outPath, fileName, mzRange[1], mzRange[2],
                                            maxNumbIter, maxTPeakW, resolution,
                                            minMSW, RThreshold, minInt, sigIncThreshold, S2N,
                                            minPeakWS),
                outPath = outPath, fileName = fileName))
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
    # UNDONE: docs
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, add = ac)
    checkmate::assertCharacter(profPath, min.chars = 1, min.len = 1, null.ok = TRUE, add = ac)
    assertCanCreateDirs(profPath, add = ac)
    checkmate::assertNumeric(mzRange, lower = 0, finite = TRUE, any.missing = FALSE, len = 2, add = ac)
    aapply(checkmate::assertCount, . ~ maxNumbIter + maxTPeakW + resolution + sigIncThreshold +
               S2N, positive = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertCount, . ~ minInt + minPeakWS, fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ minMSW + RThreshold, fixed = list(add = ac))
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    if (mzRange[1] > mzRange[2])
        stop("First element of mzRange should be smaller than second.")
    
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
