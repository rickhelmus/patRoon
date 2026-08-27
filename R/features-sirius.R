# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include features.R
NULL

loadSIRFeat <- function(json, fileIndex)
{
    # assumptions: first isotope is our target monoisotopic mass
    
    iso <- json[["traceSets"]][[fileIndex]][["ionTrace"]][["isotopes"]][[1]]
    
    featOffset <- iso$detectedFeatureOffset
    featLen <- iso$detectedFeatureLength
    
    # NOTE: +1 to make it 1 based for R vectors
    featRange <- seq(featOffset + 1, featOffset + featLen)
    rawRange <- seq(iso$indexOffset + featOffset + 1, iso$indexOffset + featOffset + featLen)
    
    # NOTE: SIRIUS retention times are in milliseconds
    eic <- data.table(time = unlist(json[["traceSets"]][[fileIndex]][["retentionTimes"]])[rawRange] / 1000,
                      intensity = unlist(iso[["intensities"]])[featRange])
    mzs <- unlist(iso[["masses"]])[featRange]
    area <- json[["abundance"]][[fileIndex]]
    ret <- eic$time[which.max(eic$intensity)] # UNDONE: verify
    
    return(data.table(ret = ret, mz = mean(mzs), mzmin = min(mzs), mzmax = max(mzs),
                      retmin = min(eic$time), retmax = max(eic$time), area = area,
                      intensity = max(eic$intensity)))
}

SIRFeatMPFinishHandler <- function(cmd)
{
    pattern <- paste0("^[[:digit:]]+_\\Q", cmd$analysis, "\\E_[[:digit:]]+$")
    resDirs <- list.files(cmd$outPath, pattern = pattern, full.names = TRUE)
    ret <- rbindlist(lapply(resDirs, function(dir)
    {
        # assumptions: only one analysis
        json <- jsonlite::fromJSON(file.path(dir, "lcms.json.gz"), FALSE)
        return(patRoon:::loadSIRFeat(json, 1))
    }))
    
    setorderv(ret, "mz") # order is inconsistent between runs --> fix order by sorting
    
    ret[, ID := as.character(seq_len(nrow(..ret)))]
    setcolorder(ret, "ID")
    return(ret[])
}

SIRFeatMPPrepareHandler <- function(cmd)
{
    command <- patRoon:::getExtDepPath("sirius")
    outPath <- tempfile("sirius_out")
    args <- c("-i", cmd$dataFile, "-o", outPath, "lcms-align")
    return(utils::modifyList(cmd, list(command = command, args = args, outPath = outPath)))
}

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
#' @inheritParams findFeatures
#'
#' @details The features are collected by running the \command{lcms-align} \command{SIRIUS} command for every analysis.
#'
#'   The MS files should be in the \file{mzML} or \file{mzXML} format. Furthermore, this algorithms requires the
#'   presence of (data-dependent) MS/MS data.
#'
#' @template centroid_note_mandatory
#'
#' @templateVar what \code{findFeaturesSIRIUS}
#' @template uses-multiProc
#' @template parallelization-cache_input
#'
#' @references \insertRef{Dhrkop2019}{patRoon}
#'
#' @inherit findFeatures return
#'
#' @export
findFeaturesSIRIUS <- function(analysisInfo, verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, fileTypes = "centroid", allowedFormats = "mzML", add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    anaCount <- nrow(analysisInfo)
    
    if (verbose)
        printf("Finding features with SIRIUS for %d analyses ...\n", anaCount)

    filePaths <- getCentroidedMSFilesFromAnaInfo(analysisInfo, "mzML")
    cmdQueue <- Map(analysisInfo$analysis, filePaths, f = function(ana, fp)
    {
        hash <- makeHash(makeFileHash(fp))
        logf <- paste0(ana, ".txt")
        return(list(hash = hash, dataFile = fp, analysis = ana, logFile = logf))
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

findFeaturesSIRIUSNew <- function(analysisInfo, noiseIntensity = NULL, traceMaxMassDeviation = NULL,
                                  minSNR = NULL, login = "check", alwaysLogin = FALSE, projectPath = NULL,
                                  runMode = "execute", SIRIUSAPI = NULL, SIRIUSPath = NULL, verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    # UNDONE: API docs say that mzXML is also supported?
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, fileTypes = "centroid", allowedFormats = "mzML", add = ac)
    aapply(checkmate::assertNumber, . ~ noiseIntensity + traceMaxMassDeviation + minSNR,
           lower = 0, finite = TRUE, null.ok = TRUE, fixed = list(add = ac))
    assertCommonSIRIUSArgs(login, alwaysLogin, projectPath, runMode, SIRIUSAPI, add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    filePaths <- getCentroidedMSFilesFromAnaInfo(analysisInfo, "mzML")
    
    hash <- makeHash(analysisInfo[, c("analysis", "path_centroid"), with = FALSE], lapply(filePaths, makeFileHash),
                     noiseIntensity, traceMaxMassDeviation, minSNR, projectPath)
    
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
        if (!is.null(traceMaxMassDeviation))
            params$traceMaxMassDeviation <- traceMaxMassDeviation
        if (!is.null(minSNR))
            params$minSNR <- minSNR
        
        job <- SIRIUSAPI$projects_api$ImportMsRunDataAsJob(projectID, as.list(unname(filePaths)), parameters = params)
        
        # NOTE: maxProgress can change during the job execution, so we normalize the current progress to it at each update
        prog <- openProgBar(0, 1)
        repeat
        {
            Sys.sleep(1)
            jp <- SIRIUSAPI$jobs_api$GetJob(projectID, job$id)$progress
            if (jp$state %in% c("CANCELED", "FAILED", "DONE"))
                break
            setTxtProgressBar(prog, jp$currentProgress / jp$maxProgress)
        }
        setTxtProgressBar(prog, 1)
        close(prog)
    }
    
    alignedFeatTab <- getSIRIUSAlignedFeatTab(SIRIUSAPI, projectID)
    SIRQuantTabInt <- getSIRIUSQuantTab(SIRIUSAPI, projectID, "intensity")
    SIRQuantTabArea <- getSIRIUSQuantTab(SIRIUSAPI, projectID, "area")
    
    allFeats <- merge(SIRQuantTabInt, alignedFeatTab, by = "SIRID", sort = FALSE)
    allFeats[SIRQuantTabArea, area := i.area, on = c("SIRID", "analysis")]
    setcolorder(allFeats, "intensity", before = "area")
    allFeats <- allFeats[intensity > 0]
    allFeats[, ID := as.character(seq_len(.N)), by = "analysis"]
    setcolorder(allFeats, "ID")
    allFeatsList <- split(allFeats, by = "analysis", keep.by = FALSE)
    allFeatsList <- allFeatsList[intersect(analysisInfo$analysis, names(allFeatsList))] # re-order
    
    ret <- featuresSIRIUS(analysisInfo = analysisInfo, features = allFeatsList)
    saveCacheData("featuresSIRIUS", ret, hash)
    
    if (verbose)
        cat("\n===========\nDone!\n")
    
    return(ret)
}
