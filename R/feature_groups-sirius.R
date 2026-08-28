# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include features.R
#' @include feature_groups.R
NULL

#' @rdname featureGroups-class
#' @export
featureGroupsSIRIUS <- setClass("featureGroupsSIRIUS", slots = c(groupInfoSIR = "data.table"),
                                contains = "featureGroups")

setMethod("initialize", "featureGroupsSIRIUS",
          function(.Object, ...) callNextMethod(.Object, algorithm = "sirius", ...))


#' Group features using SIRIUS
#'
#' Uses \href{https://bio.informatik.uni-jena.de/software/sirius/}{SIRIUS} to find \emph{and} group features.
#'
#' @templateVar algo SIRIUS
#' @templateVar do group features
#' @templateVar generic groupFeatures
#' @templateVar algoParam sirius
#' @template algo_generator
#'
#' @details This algorithm always first finds features with \command{SIRIUS} and can therefore not group features from
#'   other algorithms. The MS files should be in the \file{mzML} or \file{mzXML} format.
#'
#' @param \dots Further arguments passed to \code{\link{findFeaturesSIRIUS}} (\code{groupFeaturesSIRIUS}) or
#'   \code{\link{groupFeaturesSIRIUS}} (\code{importFeatureGroupsSIRIUS}).
#'
#' @template sirius-args
#' @template centroid_note_mandatory
#'
#' @template analysisInfo-arg
#' @inheritParams groupFeatures
#'
#' @inherit groupFeatures return
#'
#' @return Returns a \code{featureGroupsSIRIUS} object dervied from \code{\link{featureGroups}}. This object contains an
#'   additional slot \code{groupInfoSIR} with information about the SIRIUS aligned features.
#'
#' @export
groupFeaturesSIRIUS <- function(analysisInfo, ..., login = "check", alwaysLogin = FALSE, projectPath = NULL,
                                runMode = "execute", SIRIUSAPI = NULL, SIRIUSPath = NULL, verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    # UNDONE: API docs say that mzXML is also supported?
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, fileTypes = "centroid", allowedFormats = "mzML", add = ac)
    assertCommonSIRIUSArgs(login, alwaysLogin, projectPath, runMode, SIRIUSAPI, add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    filePaths <- getCentroidedMSFilesFromAnaInfo(analysisInfo, "mzML")
    hash <- makeHash(analysisInfo[, c("analysis", "path_centroid"), with = FALSE], lapply(filePaths, makeFileHash), ...,
                     if (!is.null(projectPath) && file.exists(projectPath)) makeFileHash(projectPath))
    
    cachefg <- loadCacheData("featureGroupsSIRIUS", hash)
    if (!is.null(cachefg))
        return(cachefg)

    # setup here so it can be shared with findFeaturesSIRIUS()
    if (is.null(SIRIUSAPI))
        SIRIUSAPI <- startSIRIUS(SIRIUSPath)
    
    # NOTE: setup projectPath/ID here so it can be shared with findFeaturesSIRIUS()
    if (is.null(projectPath))
        projectPath <- tempfile("patRoonSIRIUS", fileext = ".sirius")
    if (length(names(projectPath)) == 0)
        names(projectPath) <- "patRoonProjectID"
    
    # HACK: we always need to find features with SIRIUS to get group info from its project file, so temporarily disable
    # loading from cache
    features <- withOpt(
        cache.mode = "save",
        findFeaturesSIRIUS(analysisInfo, ..., login = login, alwaysLogin = alwaysLogin,
                           projectPath = projectPath, runMode = runMode, SIRIUSAPI = SIRIUSAPI,
                           SIRIUSPath = SIRIUSPath, verbose = verbose)
    )
    
    # no need to login
    # doSIRIUSLogin(login, alwaysLogin, SIRIUSAPI)
    
    projectID <- openSIRIUSProject(projectPath, SIRIUSAPI, "read")
    
    if (verbose)
        printf("Importing SIRIUS aligned features...\n")
    
    SIRAlignedFeats <- getSIRIUSPagedResults(SIRIUSAPI$features_api$GetAlignedFeaturesPage, projectID,
                                             opt_fields = "qualities", showProgress = verbose, pageSize = 100)
    SIRAlignedFeats <- rbindlist(lapply(SIRAlignedFeats, function(f)
    {
        feat <- f$toList()[c("alignedFeatureId", "rtApexSeconds", "ionMass", "quality")]
        feat$qualityIsotope = f$qualities$ISOTOPE_QUALITY
        feat$qualityPeak = f$qualities$PEAK_QUALITY
        return(feat)
    }))
    
    allFeats <- as.data.table(features)
    
    # NOTE: somehow SIRIUS can return empty aligned features
    SIRAlignedFeats <- SIRAlignedFeats[alignedFeatureId %chin% unique(allFeats$SIRAlignedFeatureID)]
    
    gInfo <- SIRAlignedFeats[, .(ret = rtApexSeconds, mz = ionMass, SIRAlignedFeatureID = alignedFeatureId)]
    setorderv(gInfo, c("ret", "mz"))
    gInfo[, group := mapply(seq_len(.N), ret, mz, FUN = makeFGroupName)]
    setcolorder(gInfo, "group")
    
    gTab <- dcast(allFeats, analysis ~ SIRAlignedFeatureID, value.var = "intensity", fill = 0)
    gTab <- gTab[match(analysisInfo$analysis, analysis)]
    gTab[, analysis := NULL]
    setnames(gTab, gInfo$group[match(names(gTab), gInfo$SIRAlignedFeatureID)])
    gTab <- gTab[, gInfo$group, with = FALSE] # sync order
    
    allFeats[, row := seq_len(.N), by = "analysis"]
    ftind <- dcast(allFeats, analysis ~ SIRAlignedFeatureID, value.var = "row", fill = 0)
    ftind <- ftind[match(analysisInfo$analysis, analysis)]
    ftind[, analysis := NULL]
    ftind[, (names(ftind)) := lapply(.SD, as.integer), .SDcols = names(ftind)]
    setnames(ftind, gInfo$group[match(names(ftind), gInfo$SIRAlignedFeatureID)])
    ftind <- ftind[, gInfo$group, with = FALSE] # sync order
    
    gInfoSIR <- gInfo[, c("group", "SIRAlignedFeatureID"), with = FALSE]
    gInfoSIR[SIRAlignedFeats, c("quality", "qualityIsotope", "qualityPeak") :=
                 .(i.quality, i.qualityIsotope, i.qualityPeak), on = c("SIRAlignedFeatureID" = "alignedFeatureId")]
    gInfo[, SIRAlignedFeatureID := NULL]
    
    ret <- featureGroupsSIRIUS(groups = gTab, groupInfo = gInfo, features = features, ftindex = ftind,
                               groupInfoSIR = gInfoSIR)
    
    saveCacheData("featureGroupsSIRIUS", ret, hash)
    
    if (verbose)
        cat("\n===========\nDone!\n")
    
    return(ret)
}

#' @details \code{importFeatureGroupsSIRIUS} is a simple wrapper around \code{groupFeaturesSIRIUS} to import feature
#'   groups from an existing SIRIUS project. It will set \code{runMode="read"} and \code{projectPath} to the provided
#'   \code{input} path.
#'
#' @param input Sets \code{projectPath}.
#'
#' @rdname groupFeaturesSIRIUS
#' @export
importFeatureGroupsSIRIUS <- function(input, analysisInfo, ...)
{
    checkmate::assertFileExists(input)
    groupFeaturesSIRIUS(analysisInfo, projectPath = input, runMode = "read", ...)
}

#' @rdname featureGroups-class
#' @export
setMethod("delete", "featureGroupsSIRIUS", function(obj, ...)
{
    obj <- callNextMethod()
    obj@groupInfoSIR <- obj@groupInfoSIR[group %chin % names(obj)]
    return(obj)
})
