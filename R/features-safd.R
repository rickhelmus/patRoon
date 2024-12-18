# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include features.R
NULL

SAFDMPFinishHandler <- function(cmd)
{
    fExt <- if (cmd$cent) "_Cent_report.csv" else "_report.csv"
    results <- fread(file.path(cmd$outPath, paste0(cmd$fileName, fExt)))
    setnames(results,
             c("Nr", "Rt", "MeasMass", "RtStart", "RtEnd", "MinMass", "MaxMass", "Area", "Int",
               "FeatPurity", "MediRes"),
             c("ID", "ret", "mz", "retmin", "retmax", "mzmin", "mzmax", "area", "intensity",
               "purity", "mediRes"))
    results[, c("ret", "retmin", "retmax") := .(ret * 60, retmin * 60, retmax * 60)] # min --> sec
    results[, c("MinInPeak", "ScanNum", "ScanInPeak") := NULL]
    
    return(results[])
}

SAFDMPPrepareHandler <- function(cmd)
{
    outp <- tempfile("safd")
    mkdirp(outp)
    cmd$args <- c(cmd$args, outp)
    return(c(cmd, list(outPath = outp)))
}

#' @rdname features-class
#' @export
featuresSAFD <- setClass("featuresSAFD", contains = "features")

setMethod("initialize", "featuresSAFD",
          function(.Object, ...) callNextMethod(.Object, algorithm = "safd", ...))

makeSAFDCommand <- function(inPath, fileName, cent, mzRange, maxNumbIter, maxTPeakW, resolution,
                            minMSW, RThreshold, minInt, sigIncThreshold, S2N, minPeakWS)
{
    # UNDONE: check if julia exists? allow to configure path?
    return(list(command = "julia", args = c(system.file("misc", "runSAFD.jl", package = "patRoon"),
                                            inPath, fileName, cent, mzRange[1], mzRange[2],
                                            maxNumbIter, maxTPeakW, resolution,
                                            minMSW, RThreshold, minInt, sigIncThreshold, S2N,
                                            minPeakWS),
                cent = cent, fileName = fileName))
}

#' Find features using SAFD
#'
#' Uses \href{https://bitbucket.org/SSamanipour/safd.jl/src/master/}{SAFD} to obtain features. This functionality is
#' still experimental. Please see the details below.
#'
#' @templateVar algo SAFD
#' @templateVar do automatically find features
#' @templateVar generic findFeatures
#' @templateVar algoParam safd
#' @template algo_generator
#'
#' @details The support for SAFD is still experimental, and its interface might change in the future.
#'
#'   In order to use SAFD, please make sure that its \code{julia} packages are installed and you have verified that
#'   everything works, \emph{e.g.} by running the test data.
#'
#'   This algorithm supports profile and centroided MS data. If the use of profile data is desired, centroided data
#'   must still be available for other functionality of \code{patRoon}. The centroided data is specified through the
#'   'regular' \link[=analysis-information]{analysis info} mechanism. The location to any profile data is specified
#'   through the \code{profPath} argument (\code{NULL} for no profile data). The base file names (\emph{i.e.} the file
#'   name without path and extension) of both centroid and profile data must be the same. Furthermore, the format of the
#'   profile data must be \file{mzXML}.
#'
#' @inheritParams findFeatures
#'
#' @param profPath A \code{character} vector with paths to the profile MS data for each analysis (will be re-cycled if
#'   necessary). See the \verb{Using SAFD} section for more details.
#' @param mzRange The \emph{m/z} window to be imported (passed to the \code{import_files_MS1} function).
#' @param maxNumbIter,maxTPeakW,resolution,minMSW,RThreshold,minInt,sigIncThreshold,S2N,minPeakWS Parameters directly
#'   passed to the \code{safd_s3D} function.
#'
#' @templateVar what \code{findFeaturesSAFD}
#' @template uses-multiProc
#' 
#' @template parallelization-cache_input
#'
#' @references \insertRef{Samanipour2019}{patRoon}
#' 
#' @inherit findFeatures return
#' 
#' @export
findFeaturesSAFD <- function(analysisInfo, profPath = NULL, mzRange = c(0, 400), 
                             maxNumbIter = 1000, maxTPeakW = 300, resolution = 30000,
                             minMSW = 0.02, RThreshold = 0.75, minInt = 2000,
                             sigIncThreshold = 5, S2N = 2, minPeakWS = 3, verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, "mzXML", verifyCentroided = TRUE, add = ac)
    if (!is.null(profPath))
    {
        checkmate::assertCharacter(profPath, min.chars = 1, min.len = 1, add = ac)
        assertCanCreateDirs(profPath, add = ac)
    }
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
    if (!is.null(profPath))
        profPath <- rep(profPath, length.out = anaCount)
    
    params <- list(is.null(profPath), mzRange, maxNumbIter, maxTPeakW, resolution, minMSW, RThreshold, minInt,
                   sigIncThreshold, S2N, minPeakWS)
    baseHash <- makeHash(params)
    
    if (verbose)
        printf("Finding features with SAFD for %d analyses ...\n", anaCount)

    anaPaths <- if (is.null(profPath)) analysisInfo$path else profPath
    cmdQueue <- Map(analysisInfo$analysis, anaPaths, f = function(ana, path)
    {
        # fpNCDF <- getAnalysisPath(ana, path, "cdf") UNDONE: also support netcdf?
        fp <- getMzXMLAnalysisPath(ana, path)
        
        if (length(fp) == 0)
            stop(sprintf("Cannot find %s or %s\n", fpMZXML, fpNCDF))
        
        # UNDONE: really want to hash big files?
        hash <- makeHash(makeFileHash(fp), params)
        
        cmd <- do.call(makeSAFDCommand, c(list(dirname(fp), basename(fp)), params))
        
        return(c(cmd, list(inputPath = fp, hash = hash, logFile = paste0(ana, ".txt"))))
    })
    
    fts <- executeMultiProcess(cmdQueue, patRoon:::SAFDMPFinishHandler,
                               prepareHandler = patRoon:::SAFDMPPrepareHandler,
                               cacheName = "featuresSAFD", logSubDir = "safd",
                               showProgress = verbose)
        
    if (verbose)
    {
        printf("Done!\n")
        printFeatStats(fts)
    }
    
    return(featuresSAFD(analysisInfo = analysisInfo, features = fts))
}
