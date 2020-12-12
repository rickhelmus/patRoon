#' @include features.R
#' @include feature_groups.R
NULL

doSIRIUSFGroups <- function(inputFiles, verbose)
{
    command <- getCommandWithOptPath(getSiriusBin(), "SIRIUS")
    outPath <- tempfile("sirius_out")
    args <- c("-i", paste0(inputFiles, collapse = ","),
              "-o", outPath,
              "lcms-align")
    
    executeCommand(command, args, stdout = if (verbose) "" else FALSE,
                   stderr = if (verbose) "" else FALSE)
    
    return(outPath)
}

processSIRIUSFGroups <- function(outPath, anaInfo)
{
    resDirs <- list.files(outPath, pattern = "^[[:digit:]]+_[[:alnum:]]+_[[:digit:]]+$",
                          full.names = TRUE)
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
        
        gNames <- mapply(seq_len(ngrp), gInfo$rts, gInfo$mzs, FUN =  makeFGroupName)
        rownames(gInfo) <- gNames
        setnames(gTab, gNames)
        setnames(ftind, gNames)
        
        return(featureGroupsSIRIUS(groups = gTab, groupInfo = gInfo, analysisInfo = anaInfo,
                                   features = features, ftindex = ftind))
    }

    return(featureGroupsSIRIUS(groups = data.table(), groupInfo = data.frame(), analysisInfo = anaInfo,
                               features = featuresSIRIUS(analysisInfo = anaInfo, features = list()),
                               ftindex = data.table()))
}

#' @rdname featureGroups-class
#' @export
featureGroupsSIRIUS <- setClass("featureGroupsSIRIUS", contains = "featureGroups")

setMethod("initialize", "featureGroupsSIRIUS",
          function(.Object, ...) callNextMethod(.Object, algorithm = "sirius", ...))


#' @rdname feature-grouping
#' @export
groupFeaturesSIRIUS <- function(analysisInfo, verbose = TRUE)
{
    # UNDONE: docs
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, "mzML", add = ac)
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
