# SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include features.R
NULL

SAFDMPFinishHandler <- function(cmd)
{
    resFile <- if (cmd$cent)
        file.path(cmd$outPath, paste0(cmd$fileName, "_Cent_report.csv")) # NOTE: here we need the full name...
    else
        file.path(cmd$outPath, paste0(tools::file_path_sans_ext(cmd$fileName), "_report.csv"))
    results <- fread(resFile)
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
    CSVPaths <- applyMSData(cmd$anaInfoRow, needTypes = if (cmd$cent) "centroid", showProgress = FALSE, func = function(ana, path, backend)
    {
        openMSReadBackend(backend, path)
        paths <- list(mzCSV = tempfile("mz", fileext = ".csv"), intCSV = tempfile("int", fileext = ".csv"),
                      rtCSV = tempfile("rt", fileext = ".csv"))
        patRoon:::makeSAFDInput(backend, paths$mzCSV, paths$intCSV, paths$rtCSV, cmd$mzRange[1], cmd$mzRange[2])
        return(paths)
    })[[1]]
    
    outp <- tempfile("safd")
    mkdirp(outp)
    cmd$args <- c(cmd$args, CSVPaths$mzCSV, CSVPaths$intCSV, CSVPaths$rtCSV, outp)
    
    return(c(cmd, list(outPath = outp)))
}

#' @rdname features-class
#' @export
featuresSAFD <- setClass("featuresSAFD", contains = "features")

setMethod("initialize", "featuresSAFD",
          function(.Object, ...) callNextMethod(.Object, algorithm = "safd", ...))

makeSAFDCommand <- function(inPath, fileName, cent, mzRange, maxNumbIter, maxTPeakW, resolution,
                            minMSW, RThreshold, minInt, sigIncThreshold, S2N, minPeakWS, centroidMethod, centroidDM)
{
    # UNDONE: check if julia exists? allow to configure path?
    return(list(command = "julia", args = c(system.file("misc", "runSAFD.jl", package = "patRoon"),
                                            inPath, fileName, cent, mzRange[1], mzRange[2],
                                            maxNumbIter, maxTPeakW, resolution,
                                            minMSW, RThreshold, minInt, sigIncThreshold, S2N,
                                            minPeakWS, centroidMethod, centroidDM),
                cent = cent, mzRange = mzRange, fileName = fileName))
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
#'   In order to use SAFD, please make sure that its \command{Julia} packages are installed and you have verified that
#'   everything works, \emph{e.g.} by running the test data with \command{SAFD}.
#'
#'   As of \pkg{patRoon} \samp{3.0}, \code{findFeaturesSAFD} uses the \link{msdata} interface instead of the
#'   \pkg{MS_Import.jl} \command{Julia} package to read HRMS data. This means that \pkg{MS_Import.jl} does not need to
#'   be installed, and all file formats supported by \link{msdata} are also supported for \command{SAFD} feature
#'   detection. This includes IMS-HRMS data, however, in that case IMS resolved spectra are summed and the IMS dimension
#'   is removed to make the data compatible for \command{SAFD}.
#'
#'   The \command{SAFD} algorithm was primarily developed to detect features in profile \emph{m/z} data, but centroided
#'   data is also supported. To use profile data, ensure that the paths are correctly set up in the
#'   \link[=analysis-information]{analysisInfo}. Furthermore, when using profile data you probably also need to specify
#'   centroided data in the \code{analysisInfo}, as \emph{e.g.} \code{\link{generateMSPeakLists}} currently does not
#'   support profile data. If IMS-HRMS data is used it is treated as profile data, as this data is typically not or
#'   partially centroided (\code{generateMSPeakLists} supports IMS-HRMS data directly).
#'
#' @inheritParams findFeatures
#'
#' @param prefCentroid Set to \code{TRUE} to prefer centroided data over other MS data specified in \code{analysisInfo}.
#'
#'   \strong{NOTE}: if \code{prefCentroid=FALSE} but the package option \option{patRoon.MS.preferIMS=TRUE} (see
#'   \code{\link{msdata}}), then centroided data will still be preferred over IMS data.
#' @param mzRange The \emph{m/z} window to be imported.
#' @param maxNumbIter,maxTPeakW,resolution,minMSW,RThreshold,minInt,sigIncThreshold,S2N,minPeakWS Parameters directly
#'   passed to the \code{safd_s3D} function.
#' @param centroidMethod,centroidDM Passed to the \code{safd_s3d_cent} function (\code{method} and \code{mdm} arguments,
#'   respectively).
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
findFeaturesSAFD <- function(analysisInfo, prefCentroid = FALSE, mzRange = c(0, 400), 
                             maxNumbIter = 1000, maxTPeakW = 300, resolution = 30000,
                             minMSW = 0.02, RThreshold = 0.75, minInt = 2000,
                             sigIncThreshold = 5, S2N = 2, minPeakWS = 3, centroidMethod = "RFM", centroidDM = 0.005,
                             verbose = TRUE)
{
    # UNDONE: for now don't support centroidMethod = "MDM" as it's a bit of a hassle to pass on the data
    # UNDONE: centroidDM default is just a dummy value, not sure what would be good...
    
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, add = ac)
    checkmate::assertFlag(prefCentroid, add = ac)
    checkmate::assertNumeric(mzRange, lower = 0, finite = TRUE, any.missing = FALSE, len = 2, add = ac)
    aapply(checkmate::assertCount, . ~ maxNumbIter + maxTPeakW + resolution + sigIncThreshold +
               S2N, positive = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertCount, . ~ minInt + minPeakWS, fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ minMSW + RThreshold, fixed = list(add = ac))
    checkmate::assertChoice(centroidMethod, c("RFM", "BG", "CT"), add = ac)
    checkmate::assertNumber(centroidDM, lower = 0, finite = TRUE, add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    if (mzRange[1] > mzRange[2])
        stop("First element of mzRange should be smaller than second.")

    filePaths <- NULL
    takeCent <- FALSE
    
    filePathsCentroid <- getCentroidedMSFilesFromAnaInfo(analysisInfo, mustExist = FALSE)
    if (prefCentroid && !is.null(filePathsCentroid))
    {
        filePaths <- filePathsCentroid
        takeCent <- TRUE
    }
    else
    {
        filePaths <- getMSFilesFromAnaInfo(analysisInfo, getMSFileTypes(), lapply(getMSFileTypes(), getMSFileFormats))
        takeCent <- isTRUE(all.equal(filePaths, filePathsCentroid))
    }
    
    anaCount <- nrow(analysisInfo)
    
    params <- list(takeCent, mzRange, maxNumbIter, maxTPeakW, resolution, minMSW, RThreshold, minInt,
                   sigIncThreshold, S2N, minPeakWS, centroidMethod, centroidDM)
    
    if (verbose)
        printf("Finding features with SAFD for %d analyses ...\n", anaCount)

    anaHashes <- getMSFileHashesFromAvailBackend(analysisInfo, needTypes = if (takeCent) "centroid")
    
    cmdQueue <- Map(split(analysisInfo, seq_len(nrow(analysisInfo))), filePaths, anaHashes, f = function(ai, fp, hash)
    {
        hash <- makeHash(hash, params)
        cmd <- do.call(makeSAFDCommand, c(list(dirname(fp), basename(fp)), params))
        return(c(cmd, list(anaInfoRow = ai, inputPath = fp, hash = hash, logFile = paste0(ai$analysis, ".txt"))))
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
