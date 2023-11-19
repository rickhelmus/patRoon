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
        fList <- lapply(fList, set, j = "group", value = NULL)
        features <- featuresSIRIUS(analysisInfo = anaInfo, features = fList)
        
        ngrp <- max(resTbl$group)
        gTab <- data.table(matrix(0, nrow = nrow(anaInfo), ncol = ngrp))
        ftind <- copy(gTab)
        gInfo <- data.frame(rts = numeric(ngrp), mzs = numeric(ngrp))
        
        for (grpi in seq_len(ngrp))
        {
            grpRes <- resTbl[group == grpi]
            ainds <- match(grpRes$analysis, anaInfo$analysis)
            set(gTab, ainds, j = grpi, value = grpRes$intensity)
            set(ftind, ainds, j = grpi, value = grpRes$ID)
            
            # UNDONE: does SIRIUS report group rets/mzs?
            gInfo[grpi, c("rts", "mzs")] <- list(mean(grpRes$ret), mean(grpRes$mz))
        }

        # group order is not consistent between runs --> sort
        ord <- order(gInfo$mzs)
        gInfo <- gInfo[ord, ]
        gTab <- gTab[, ord, with = FALSE]; ftind <- ftind[, ord, with = FALSE]

        gNames <- mapply(seq_len(ngrp), gInfo$rts, gInfo$mzs, FUN = makeFGroupName)
        rownames(gInfo) <- gNames
        setnames(gTab, gNames)
        setnames(ftind, gNames)

        return(featureGroupsSIRIUS(groups = gTab, groupInfo = gInfo, features = features, ftindex = ftind))
    }

    return(featureGroupsSIRIUS(groups = data.table(), groupInfo = data.frame(),
                               features = featuresSIRIUS(analysisInfo = anaInfo, features = list()),
                               ftindex = data.table()))
}

#' @rdname featureGroups-class
#' @export
featureGroupsSIRIUS <- setClass("featureGroupsSIRIUS", contains = "featureGroups")

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
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, "mzML", verifyCentroided = TRUE, add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    inputFiles <- mapply(analysisInfo$analysis, analysisInfo$path, FUN = getMzMLAnalysisPath)
    
    hash <- makeHash(analysisInfo, lapply(inputFiles, makeFileHash))
    
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
