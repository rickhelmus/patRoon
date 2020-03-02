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
    
    return(fileName)
}

unifyMSFNames <- function(resTab)
{
    unNames <- c(NAME = "compoundName",
                 ID = "identifier",
                 INCHIKEY = "InChIKey",
                 SMILES = "SMILES",
                 RESOURCES = "resources",
                 SubstructureInChIKeys = "substructureInChIKeys",
                 SubstructureOntologies = "substructureOntologies",
                 Ontology = "ontology",
                 OntologyID = "ontologyID",
                 TotalBondEnergy = "totalBondEnergy",
                 TotalScore = "score",
                 TotalHrRulesScore = "totalHrRulesScore",
                 TotalBondCleavageScore = "totalBondCleavageScore",
                 TotalMassAccuracyScore = "totalMassAccuracyScore",
                 TotalFragmentLinkageScore = "totalFragmentLinkageScore",
                 TotalBondDissociationEnergyScore = "totalBondDissociationEnergyScore",
                 DatabaseScore = "databaseScore",
                 SubstructureAssignmentScore = "substructureAssignmentScore",
                 
                 fragInfo = "fragInfo" # keep
    )
    
    unNames <- unNames[names(unNames) %in% names(resTab)] # filter out missing
    setnames(resTab, names(unNames), unNames)
    
    return(resTab[, unNames, with = FALSE]) # filter out any other columns
}

readSFD <- function(file)
{
    lines <- readLines(file)
    lines <- lines[nzchar(lines)]
    
    # split in result chunks, where it's assumed each result starts with NAME field
    linesResults <- split(lines, cumsum(grepl("^NAME:", lines)))
    
    # assume file names are named by their candidate formula
    candidateFormula <- tools::file_path_sans_ext(basename(file))
    
    results <- lapply(linesResults, function(lr)
    {
        # split by variable/value by colon, see https://stackoverflow.com/a/26247455
        parsed <- regmatches(lr, regexpr(":", lr), invert = TRUE)
        # omit non-key/value pairs (these are usually the annotated fragments)
        parsed <- parsed[lengths(parsed) == 2]
        parsed <- lapply(parsed, function(p) trimws(p))
        
        keys <- sapply(parsed, "[", 1)
        vals <- lapply(parsed, "[", 2) # lapply: keep it as a list
        
        ret <- data.table()
        ret[, (keys) := vals]
        
        annColName <- grep("^Num Fragment", keys, value = TRUE)
        annRows <- as.integer(ret[[annColName]])
        
        if (annRows > 0)
        {
            annStartRow <- which(grepl(annColName, lr, fixed = TRUE)) + 1
            # parse raw text from annotated fragment as a table
            fragInfo <- fread(text = paste0(lr[seq(annStartRow, annStartRow + (annRows-1))], collapse = "\n"))
            # figure out column names from key value (white space separated between parenthesis)
            cols <- sub(".+\\((.+)\\)", "\\1", annColName)
            cols <- unlist(strsplit(cols, " ", fixed = TRUE))
            setnames(fragInfo, cols)
            
            # unify names
            setnames(fragInfo,
                     c("M/Z", "Intensity", "Formula", "TotalScore"),
                     c("mz", "intensity", "formula", "score"))
            
            # UNDONE: keep more info? There's plenty...
            fragInfo <- fragInfo[, c("mz", "intensity", "formula", "score", "SMILES"), with = FALSE]
            fragInfo[, neutral_loss := sapply(formula, subtractFormula, formula1 = candidateFormula)]
            
            ret[, fragInfo := list(list(fragInfo))]
        }
        
        
        ret <- unifyMSFNames(ret)
        ret[, formula := candidateFormula]
        
        browser()
    })
    
    
}

readMSFResults <- function(outDir, MATFiles)
{
    # resDirs <- list.dirs(outDir, recursive = FALSE)
    # resDirs <- resDirs[grepl(tools::file_path_sans_ext(basename(MATFiles)), resDirs, fixed = TRUE)]
    resDirs <- tools::file_path_sans_ext(MATFiles)
    resDirs <- resDirs[dir.exists(resDirs)]
    
    for (rd in resDirs)
    {
        # find *.sfd files
        # process with readSFD and merge
    }
    
    ret <- rbindlist
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
    
    
    MATFiles <- sapply(gNames, function(gName)
    {
        # Prepare input files
        
        if (is.null(MSPeakLists[[gName]][["MS"]]) || is.null(MSPeakLists[[gName]][["MSMS"]]))
            return(character())

        return(generateMATFile(gName, gInfo[gName, "mzs"], adduct, MSPeakLists[[gName]], outDir))
    })
    
    # Run command (or split in multiproc batches?)
    # withr::with_dir(normalizePath("~/Rproj/MSFINDER ver 3.30"),
    #                 executeCommand("MsfinderConsoleApp.exe",
    #                c("predict", "-i", outDir, "-o", file.path(outDir, "results"), "-m", file.path(outDir, "MSFINDER.INI"))))
    
    # work-around from https://stackoverflow.com/a/50512416
    withr::with_dir(normalizePath("~/Rproj/MSFINDER-ver-3_30"),
                    shell(sprintf("start cmd /C MsfinderConsoleApp.exe predict -i %s -o %s -m %s",
                                  outDir, file.path(outDir, "results"), file.path(outDir, "MSFINDER.INI"))))
    
    # executeMultiProcess(list(list(command = normalizePath("~/Rproj/MSFINDER ver 3.30/MsfinderConsoleApp.exe"),
    #                               args = c("predict", "-i", outDir, "-o", file.path(outDir, "results"), "-m", file.path(outDir, "MSFINDER.INI")),
    #                               logFile = "msf-log.txt")),
    #                     function(cmd) NULL)
    
    # Read in results

    
        
    browser()
    
    ngrp <- length(ret)
    printf("Loaded %d compounds from %d features (%.2f%%).\n", sum(unlist(lapply(ret, function(r) nrow(r$comptab)))),
           ngrp, if (gCount == 0) 0 else ngrp * 100 / gCount)
    
    # if (is.null(cachedSet))
    #     saveCacheSet("compoundsMSFINDER", resultHashes[resultHashes != ""], setHash, cacheDB)
    
    return(compounds(compounds = lapply(ret, "[[", "comptab"), scoreTypes = "score",
                     scoreRanges = lapply(ret, "[[", "scRanges"),
                     algorithm = "sirius"))
}
