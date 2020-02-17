#' @include main.R
#' @include compounds.R
NULL

generateMATFile <- function(gName, precMz, adduct, plistsFG, out)
{
    # format: see https://mtbinfo-team.github.io/mtbinfo.github.io/MS-FINDER/tutorial#section-3
    
    fileName <- tempfile("spec", fileext = ".mat", tmpdir = out)
    fileCon <- withr::local_connection(file(fileName, "w"))
    writeParam <- function(param, val) cat(sprintf("%s: %s\n", param, val), file = fileCon)
    
    writeParam("NAME", gName)
    writeParam("PRECURSORMZ", precMz)
    writeParam("PRECURSORTYPE", as.character(adduct)) # UNDONE: check format compat
    writeParam("IONTYPE", if (adduct@charge >= 1) "Positive" else "Negative")
    
    writeParam("MSTYPE", "MS1")
    writeParam("Num Peaks", nrow(plistsFG[["MS"]]))
    cat(paste(plistsFG[["MS"]]$mz, plistsFG[["MS"]]$intensity, sep = "\t"), sep = "\n", file = fileCon)
    
    # UNDONE: MS2 is required?
    writeParam("MSTYPE", "MS2")
    writeParam("Num Peaks", nrow(plistsFG[["MSMS"]]))
    cat(paste(plistsFG[["MSMS"]]$mz, plistsFG[["MSMS"]]$intensity, sep = "\t"), sep = "\n", file = fileCon)
    
    return(fileName) # UNDONE: needed?
}

#' @rdname compound-generation
#' @export
generateCompoundsMSFINDER <- function(fGroups, MSPeakLists, adduct = "[M+H]+", elements = "CHNOP",
                                      extraOpts = NULL,
                                      logPath = file.path("log", "msfinder"),
                                      maxProcAmount = getOption("patRoon.maxProcAmount"))
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertString(elements, add = ac)
    checkmate::assertList(extraOpts, null.ok = TRUE, add = ac)
    assertMultiProcArgs(logPath, maxProcAmount, add = ac)
    checkmate::reportAssertions(ac)
    
    adduct <- checkAndToAdduct(adduct)
    
    anaInfo <- analysisInfo(fGroups)
    gInfo <- groupInfo(fGroups)
    gNames <- names(fGroups)
    gCount <- length(fGroups)
    
    cacheDB <- openCacheDBScope()
    setHash <- makeHash(fGroups, MSPeakLists, adduct, elements, extraOpts)
    cachedSet <- loadCacheSet("compoundsSirius", setHash, cacheDB)
    resultHashes <- vector("character", gCount)
    names(resultHashes) <- gNames
    
    printf("Processing %d feature groups with MS-FINDER...\n", gCount)

    # Prepare INI
    
    outDir <- tempfile("msfinder")
    mkdirp(outDir)
    file.copy(system.file("misc", "MSFINDER.INI", package = "patRoon"), outDir)
    # Modify INI
    
    
    for (gName in gNames)
    {
        # Prepare input files
        
        if (is.null(MSPeakLists[[gName]][["MS"]]) || is.null(MSPeakLists[[gName]][["MSMS"]]))
            next

        generateMATFile(gName, gInfo[gName, "mzs"], adduct, MSPeakLists[[gName]], outDir)
    }
    
    # Run command (or split in multiproc batches?)
    # withr::with_dir(normalizePath("~/Rproj/MSFINDER ver 3.30"),
    #                 executeCommand("MsfinderConsoleApp.exe",
    #                c("predict", "-i", outDir, "-o", file.path(outDir, "results"), "-m", file.path(outDir, "MSFINDER.INI"))))
    executeMultiProcess(list(list(command = normalizePath("~/Rproj/MSFINDER ver 3.30/MsfinderConsoleApp.exe"),
                                  args = c("predict", "-i", outDir, "-o", file.path(outDir, "results"), "-m", file.path(outDir, "MSFINDER.INI")),
                                  logFile = "msf-log.txt")),
                        function(cmd) NULL)
    
    browser()
    # Read in results
    
    
    ngrp <- length(ret)
    printf("Loaded %d compounds from %d features (%.2f%%).\n", sum(unlist(lapply(ret, function(r) nrow(r$comptab)))),
           ngrp, if (gCount == 0) 0 else ngrp * 100 / gCount)
    
    # if (is.null(cachedSet))
    #     saveCacheSet("compoundsMSFINDER", resultHashes[resultHashes != ""], setHash, cacheDB)
    
    return(compounds(compounds = lapply(ret, "[[", "comptab"), scoreTypes = "score",
                     scoreRanges = lapply(ret, "[[", "scRanges"),
                     algorithm = "sirius"))
}
