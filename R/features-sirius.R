# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
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
    pattern <- paste0("^[[:digit:]]+_", cmd$analysis, "_[[:digit:]]+$")
    resDirs <- list.files(cmd$outPath, pattern = pattern, full.names = TRUE)
    ret <- rbindlist(lapply(resDirs, function(dir)
    {
        # assumptions: only one analysis
        json <- jsonlite::fromJSON(file.path(dir, "lcms.json.gz"), FALSE)
        return(patRoon:::loadSIRFeat(json, 1))
    }))
    
    setorderv(ret, "mz") # order is inconsistent between runs --> fix order by sorting
    
    ret[, ID := seq_len(nrow(..ret))]
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
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, verifyCentroided = TRUE, "mzML", add = ac)
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
