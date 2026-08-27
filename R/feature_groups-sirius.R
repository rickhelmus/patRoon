# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include features.R
#' @include feature_groups.R
NULL

doSIRIUSFGroups <- function(inputFiles, verbose)
{
    command <- getExtDepPath("sirius")
    outPath <- tempfile("sirius_out")
    args <- c("-i", paste0(inputFiles, collapse = ","),
              "-o", outPath,
              "lcms-align")
    
    executeCommand(command, args, stdout = if (verbose) "" else FALSE, stderr = if (verbose) "" else FALSE)
    
    return(outPath)
}

processSIRIUSFGroups <- function(outPath, anaInfo)
{
    resDirs <- list.files(outPath, pattern = "^[[:digit:]]+_.+_[[:digit:]]+$", full.names = TRUE)
    
    resTbl <- rbindlist(Map(resDirs, seq_along(resDirs), f = function(dir, grpi)
    {
        json <- jsonlite::fromJSON(file.path(dir, "lcms.json.gz"), FALSE)
        anas <- tools::file_path_sans_ext(unlist(json[["sampleNames"]]))
        feats <- setNames(lapply(seq_along(anas), loadSIRFeat, json = json), anas)
        feats <- rbindlist(feats, idcol = "analysis")
        feats[, group := grpi]
        return(feats)
    }))

    if (nrow(resTbl) > 0)
    {
        resTbl[, ID := seq_len(.N), by = "analysis"]
        fList <- split(resTbl, by = "analysis", keep.by = FALSE)
        fList <- fList[intersect(anaInfo$analysis, names(fList))] # re-order
        # no need anymore, and clashes with group assignments in fGroups constructor
        fList <- lapply(fList, function(fl)
        {
            set(fl, j = c("ID", "group"), value = list(as.character(fl$ID), NULL))
        })
        features <- featuresSIRIUS(analysisInfo = anaInfo, features = fList)
        
        ngrp <- max(resTbl$group)
        gTab <- data.table(matrix(0, nrow = nrow(anaInfo), ncol = ngrp))
        ftind <- copy(gTab)
        gInfo <- data.table(ret = numeric(ngrp), mz = numeric(ngrp))
        
        for (grpi in seq_len(ngrp))
        {
            grpRes <- resTbl[group == grpi]
            ainds <- match(grpRes$analysis, anaInfo$analysis)
            set(gTab, ainds, j = grpi, value = grpRes$intensity)
            set(ftind, ainds, j = grpi, value = grpRes$ID)
            
            # UNDONE: does SIRIUS report group rets/mzs?
            gInfo[grpi, c("ret", "mz") := .(mean(grpRes$ret), mean(grpRes$mz))]
        }

        # group order is not consistent between runs --> sort
        ord <- order(gInfo$mz)
        gInfo <- gInfo[ord]
        gTab <- gTab[, ord, with = FALSE]; ftind <- ftind[, ord, with = FALSE]

        gNames <- mapply(seq_len(ngrp), gInfo$ret, gInfo$mz, FUN = makeFGroupName)
        gInfo[, group := gNames]
        setcolorder(gInfo, "group")
        setnames(gTab, gNames)
        setnames(ftind, gNames)

        return(featureGroupsSIRIUS(groups = gTab, groupInfo = gInfo, features = features, ftindex = ftind))
    }

    return(featureGroupsSIRIUS(groups = data.table(), groupInfo = data.table(),
                               features = featuresSIRIUS(analysisInfo = anaInfo, features = list()),
                               ftindex = data.table()))
}

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
#' @details Finding and grouping features is done by running the \command{lcms-align} command on every analyses at once.
#'   For this reason, grouping feature data from other algorithms than \command{SIRIUS} is not supported.
#'
#'   The MS files should be in the \file{mzML} or \file{mzXML} format. Furthermore, this algorithms requires the
#'   presence of (data-dependent) MS/MS data.
#'
#' @template centroid_note_mandatory
#'
#' @template analysisInfo-arg
#' @inheritParams groupFeatures
#'
#' @inherit groupFeatures return
#'
#' @references \insertRef{Dhrkop2019}{patRoon}
#'
#' @export
groupFeaturesSIRIUS <- function(analysisInfo, verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, fileTypes = "centroid", allowedFormats = "mzML", add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    inputFiles <- getCentroidedMSFilesFromAnaInfo(analysisInfo, "mzML")
    
    hash <- makeHash(analysisInfo[, c("analysis", "path_centroid"), with = FALSE], lapply(inputFiles, makeFileHash))
    
    cachefg <- loadCacheData("featureGroupsSIRIUS", hash)
    if (!is.null(cachefg))
        return(cachefg)

    if (verbose)
        cat("Grouping features with SIRIUS...\n===========\n")

    outPath <- doSIRIUSFGroups(inputFiles, verbose)
    ret <- processSIRIUSFGroups(outPath, analysisInfo)
    
    saveCacheData("featureGroupsSIRIUS", ret, hash)

    if (verbose)
        cat("\n===========\nDone!\n")

    return(ret)
}

groupFeaturesSIRIUSNew <- function(analysisInfo, ..., login = "check", alwaysLogin = FALSE, projectPath = NULL,
                                   runMode = "execute", SIRIUSAPI = NULL, SIRIUSPath = NULL, verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    # UNDONE: API docs say that mzXML is also supported?
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, fileTypes = "centroid", allowedFormats = "mzML", add = ac)
    assertCommonSIRIUSArgs(login, alwaysLogin, projectPath, runMode, SIRIUSAPI, add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    filePaths <- getCentroidedMSFilesFromAnaInfo(analysisInfo, "mzML")
    hash <- makeHash(analysisInfo[, c("analysis", "path_centroid"), with = FALSE], lapply(filePaths, makeFileHash), ...)
    
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
        findFeaturesSIRIUSNew(analysisInfo, ..., login = login, alwaysLogin = alwaysLogin,
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
    setorderv(gInfo, "mz")
    gInfo[, group := mapply(seq_len(.N), ret, mz, FUN = makeFGroupName)]
    setcolorder(gInfo, "group")
    
    gTab <- dcast(allFeats, analysis ~ SIRAlignedFeatureID, value.var = "intensity", fill = 0)
    gTab <- gTab[match(analysisInfo$analysis, analysis)]
    gTab[, analysis := NULL]
    setnames(gTab, gInfo$group)
    
    allFeats[, row := seq_len(.N), by = "analysis"]
    ftind <- dcast(allFeats, analysis ~ SIRAlignedFeatureID, value.var = "row", fill = 0)
    ftind <- ftind[match(analysisInfo$analysis, analysis)]
    ftind[, analysis := NULL]
    ftind[, (names(ftind)) := lapply(.SD, as.integer), .SDcols = names(ftind)]
    setnames(ftind, gInfo$group)
    
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

#' @rdname featureGroups-class
#' @export
setMethod("delete", "featureGroupsSIRIUS", function(obj, ...)
{
    obj <- callNextMethod()
    obj@groupInfoSIR <- obj@groupInfoSIR[group %chin % names(obj)]
    return(obj)
})
