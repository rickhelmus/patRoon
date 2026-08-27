# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include features.R
NULL

#' @rdname features-class
#' @export
featuresSIRIUS <- setClass("featuresSIRIUS", contains = "features")

setMethod("initialize", "featuresSIRIUS",
          function(.Object, ...) callNextMethod(.Object, algorithm = "sirius", ...))

#' Find features using SIRIUS
#'
#' Uses \href{https://bio.informatik.uni-jena.de/software/sirius/}{SIRIUS} to find features.
#'
#' @templateVar algo SIRIUS
#' @templateVar do automatically find features
#' @templateVar generic findFeatures
#' @templateVar algoParam sirius
#' @template algo_generator
#'
#' @param noiseIntensity,alignMaxRTDev,minSNR Parameters for the SIRIUS feature finding algorithm. See the SIRIUS
#'   documentation for details. Set to \code{NULL} to use the default values.
#'
#' @inheritParams findFeatures
#'
#' @details The MS files should be in the \file{mzML} or \file{mzXML} format.
#'
#' @template sirius-args
#' @template centroid_note_mandatory
#'
#' @inherit findFeatures return
#'
#' @export
findFeaturesSIRIUS <- function(analysisInfo, noiseIntensity = NULL, alignMaxRTDev = NULL,
                               minSNR = NULL, login = "check", alwaysLogin = FALSE, projectPath = NULL,
                               runMode = "execute", SIRIUSAPI = NULL, SIRIUSPath = NULL, verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    # UNDONE: API docs say that mzXML is also supported?
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, fileTypes = "centroid", allowedFormats = "mzML", add = ac)
    aapply(checkmate::assertNumber, . ~ noiseIntensity + alignMaxRTDev + minSNR,
           lower = 0, finite = TRUE, null.ok = TRUE, fixed = list(add = ac))
    assertCommonSIRIUSArgs(login, alwaysLogin, projectPath, runMode, SIRIUSAPI, add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    filePaths <- getCentroidedMSFilesFromAnaInfo(analysisInfo, "mzML")
    
    hash <- makeHash(analysisInfo[, c("analysis", "path_centroid"), with = FALSE], lapply(filePaths, makeFileHash),
                     noiseIntensity, alignMaxRTDev, minSNR)
    
    cachefg <- loadCacheData("featuresSIRIUS", hash)
    if (!is.null(cachefg))
        return(cachefg)
    
    if (is.null(SIRIUSAPI))
        SIRIUSAPI <- startSIRIUS(SIRIUSPath)
    
    doSIRIUSLogin(login, alwaysLogin, SIRIUSAPI)
    projectID <- openSIRIUSProject(projectPath, SIRIUSAPI, runMode)
    
    if (runMode == "execute")
    {
        if (verbose)
            printf("Running SIRIUS job to get features...\n")
        
        params <- RSirius::LcmsSubmissionParameters$new()
        if (!is.null(noiseIntensity))
            params$noiseIntensity <- noiseIntensity
        if (!is.null(alignMaxRTDev))
            params$alignMaxRetentionTimeDeviation <- alignMaxRTDev
        if (!is.null(minSNR))
            params$minSNR <- minSNR
        
        job <- SIRIUSAPI$projects_api$ImportMsRunDataAsJob(projectID, as.list(unname(filePaths)), parameters = params)
        
        # NOTE: maxProgress can change during the job execution, so we normalize the current progress to it at each update
        prog <- if (verbose) openProgBar(0, 1)
        repeat
        {
            Sys.sleep(1)
            jp <- SIRIUSAPI$jobs_api$GetJob(projectID, job$id)$progress
            if (jp$state %in% c("CANCELED", "FAILED", "DONE"))
                break
            if (verbose)
                setTxtProgressBar(prog, jp$currentProgress / jp$maxProgress)
        }
        if (verbose)
        {
            setTxtProgressBar(prog, 1)
            close(prog)
        }
    }
    
    if (verbose)
        printf("Importing SIRIUS features...\n")
    
    SIRQT <- SIRIUSAPI$features_api$GetFeatureQuantTable(projectID, opt_fields = "columnSources")
    SIRAlignedFeatIDs <- unlist(SIRQT$rowIds)
    SIRAnas <- baseName(tools::file_path_sans_ext(unlist(SIRQT$columnSources)))
    names(SIRAnas) <- unlist(SIRQT$columnIds)
    allFeats <- rbindlist(doMap(FALSE, SIRAlignedFeatIDs, f = function(SIRFeatID)
    {
        # UNDONE: no qualities yet?
        res <- getSIRIUSPagedResults(SIRIUSAPI$features_api$GetFeaturesPage, projectID, SIRFeatID,
                                     opt_fields = "qualities", showProgress = FALSE)
        rbindlist(lapply(res, \(x) x$toList()))
    }, stripEnv = FALSE, prog = verbose))
    
    allFeats[, analysis := SIRAnas[runId]]
    allFeats[, runId := NULL]
    setnames(allFeats,
             c("featureId", "alignedFeatureId", "averageMz", "apexMz", "rtStartSeconds", "rtEndSeconds",
               "rtApexSeconds", "rtFwhmSeconds", "apexIntensity", "areaUnderCurve"),
             c("ID", "SIRAlignedFeatureID", "mz", "mzAPEX", "retmin", "retmax", "ret", "FWHM", "intensity",
               "area"))
    allFeats[, c("mzmin", "mzmax") := list(mz - 0.005, mz + 0.005)] # UNDONE: change this whenever this data becomes available
    setcolorder(allFeats, c("ID", "SIRAlignedFeatureID", "ret", "retmin", "retmax", "mz", "mzmin", "mzmax"))
    
    allFeatsList <- split(allFeats, by = "analysis", keep.by = FALSE)
    allFeatsList <- allFeatsList[intersect(analysisInfo$analysis, names(allFeatsList))] # re-order
    
    ret <- featuresSIRIUS(analysisInfo = analysisInfo, features = allFeatsList)
    saveCacheData("featuresSIRIUS", ret, hash)
    
    if (verbose)
        cat("\n===========\nDone!\n")
    
    return(ret)
}
