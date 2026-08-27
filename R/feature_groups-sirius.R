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
featureGroupsSIRIUS <- setClass("featureGroupsSIRIUS", slots = c(SIRIDMapping = "data.table"),
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

groupFeaturesSIRIUSNew <- function(analysisInfo, ..., verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    # UNDONE: API docs say that mzXML is also supported?
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, fileTypes = "centroid", allowedFormats = "mzML", add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    features <- findFeaturesSIRIUSNew(analysisInfo, ..., verbose = verbose)
    hash <- makeHash(features)
        
    cachefg <- loadCacheData("featureGroupsSIRIUS", hash)
    if (!is.null(cachefg))
        return(cachefg)

    allFeats <- as.data.table(features)
    
    # NOTE: feat properties are actually group properties at the moment, so just pick them from the features data
    gInfo <- unique(allFeats, by = "SIRID")[, c("ret", "mz", "SIRID"), with = FALSE]
    setorderv(gInfo, "mz")
    gInfo[, group := mapply(seq_len(.N), ret, mz, FUN = makeFGroupName)]
    setcolorder(gInfo, "group")
    
    gTab <- dcast(allFeats, analysis ~ SIRID, value.var = "intensity", fill = 0)
    gTab <- gTab[match(analysisInfo$analysis, analysis)]
    gTab[, analysis := NULL]
    setnames(gTab, gInfo$group)
    
    ftind <- dcast(allFeats, analysis ~ SIRID, value.var = "ID", fill = 0)
    ftind <- ftind[match(analysisInfo$analysis, analysis)]
    ftind[, analysis := NULL]
    ftind[, (names(ftind)) := lapply(.SD, as.integer), .SDcols = names(ftind)]
    setnames(ftind, gInfo$group)
    
    SIRIDMapping <- gInfo[, c("group", "SIRID"), with = FALSE]
    gInfo[, SIRID := NULL]
    
    ret <- featureGroupsSIRIUS(groups = gTab, groupInfo = gInfo, features = features, ftindex = ftind,
                               SIRIDMapping = SIRIDMapping)
    
    
    
    if (F)
    {
        features <- findFeaturesSIRIUSNew(analysisInfo, login = login, alwaysLogin = alwaysLogin,
                                          projectPath = projectPath, runMode = runMode, SIRIUSAPI = SIRIUSAPI,
                                          SIRIUSPath = SIRIUSPath, ..., verbose = verbose)
        
        if (verbose)
            printf("Running SIRIUS job to find and group features...\n")
        
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
        
        SIRFeats <- SIRIUSAPI$features_api$GetAlignedFeatures(projectID)
        alignedFeatTab <- rbindlist(lapply(SIRFeats, function(f)
        {
            return(data.table(SIRID = f$alignedFeatureId, ret = f$rtApexSeconds, mz = f$ionMass,
                              mzmin = f$ionMass - 0.005, mzmax = f$ionMass + 0.005, # UNDONE
                              retmin = f$rtStartSeconds, retmax = f$rtEndSeconds))
        }))
        
        getSirQuant <- function(type)
        {
            sirtype <- if (type == "intensity") "APEX_INTENSITY" else "AREA_UNDER_CURVE"
            # based on https://github.com/sirius-ms/sirius-client-openAPI/issues/188#issuecomment-5423823645
            sirq <- SIRIUSAPI$features_api$GetFeatureQuantTable(projectID, type = sirtype, opt_fields = "columnSources")
            tab <- rbindlist(lapply(sirq$values, \(v) as.list(unlist(v))))
            anas <- baseName(tools::file_path_sans_ext(unlist(sirq$columnSources)))
            setnames(tab, anas)
            setnafill(tab, fill = 0)
            tab[, SIRID := unlist(sirq$rowIds)]
            tab <- melt(tab, id.vars = "SIRID", variable.name = "analysis", variable.factor = FALSE, value.name = type)
            return(tab)
        }
        
        SIRQuantTabInt <- getSirQuant("intensity")
        SIRQuantTabArea <- getSirQuant("area")
        
        allFeats <- merge(SIRQuantTabInt, alignedFeatTab, by = "SIRID", sort = FALSE)
        allFeats[SIRQuantTabArea, area := i.area, on = c("SIRID", "analysis")]
        setcolorder(allFeats, "intensity", before = "area")
        allFeats <- allFeats[intensity > 0]
        allFeats[, ID := as.character(seq_len(.N)), by = "analysis"]
        setcolorder(allFeats, "ID")
        allFeatsList <- split(allFeats, by = "analysis", keep.by = FALSE)
        allFeatsList <- allFeatsList[intersect(analysisInfo$analysis, names(allFeatsList))] # re-order
        
        gInfo <- alignedFeatTab[, c("ret", "mz", "SIRID"), with = FALSE]
        gInfo[, group := mapply(seq_len(.N), ret, mz, FUN = makeFGroupName)]
        setcolorder(gInfo, "group")
        setorderv(gInfo, "mz")
        
        gTab <- dcast(allFeats, analysis ~ SIRID, value.var = "intensity", fill = 0)
        gTab <- gTab[match(analysisInfo$analysis, analysis)]
        gTab[, analysis := NULL]
        setnames(gTab, gInfo$group)
        
        ftind <- dcast(allFeats, analysis ~ SIRID, value.var = "ID", fill = 0)
        ftind <- ftind[match(analysisInfo$analysis, analysis)]
        ftind[, analysis := NULL]
        ftind[, (names(ftind)) := lapply(.SD, as.integer), .SDcols = names(ftind)]
        setnames(ftind, gInfo$group)
        
        SIRIDMapping <- gInfo[, c("group", "SIRID"), with = FALSE]
        gInfo[, SIRID := NULL]
        
        ret <- featureGroupsSIRIUS(groups = gTab, groupInfo = gInfo, features = allFeatsList, ftindex = ftind,
                                   SIRIDMapping = SIRIDMapping)
    }
        
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
    obj@SIRIDMapping <- obj@SIRIDMapping[group %chin % names(obj)]
    return(obj)
})
