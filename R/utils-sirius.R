#' @include main.R

getSiriusBin <- function()
{
    # NOTE: this seems to fluctuate every other SIRIUS version...
    return("sirius")
}

isSIRIUS5 <- function()
{
    out <- executeCommand(patRoon:::getCommandWithOptPath(patRoon:::getSiriusBin(), "SIRIUS"), "--version",
                          stdout = TRUE)
    return(any(grepl("^(SIRIUS 5\\.)", out)))
}

getSIRIUSCmpName <- function() "unknownCompound"

getSiriusResultPath <- function(outPath, msFName)
{
    # format is resultno_specname_compoundname, older versions start with 1, newer with 0
    msFName <- basename(tools::file_path_sans_ext(msFName))
    return(list.files(outPath, pattern = sprintf("[0-9]+_%s_%s", msFName, getSIRIUSCmpName()), full.names = TRUE))
}

getAndPrepareSIRIUSFragFiles <- function(resultPath)
{
    # NOTE: SIRIUS 5 packs spectra --> unzip them
    spPath <- file.path(resultPath, "spectra")
    if (file.exists(spPath) && !file.info(spPath, extra_cols = FALSE)$isdir)
    {
        exDir <- paste0(spPath, "-unz")
        unzip(spPath, exdir = exDir)
        spPath <- exDir
    }
    
    pat <- "([A-Za-z0-9]+).*\\.tsv"
    return(list.files(spPath, full.names = TRUE, pattern = pat))
}

getFormulaFromSiriusFragFile <- function(ffile)
{
    pat <- "([A-Za-z0-9]+).*\\.tsv"
    return(gsub(pat, "\\1", basename(ffile)))
}

makeSirMSFile <- function(plistMS, plistMSMS, parentMZ, adduct, out)
{
    msFile <- file(out, "w")

    writeMeta <- function(var, data) cat(sprintf(">%s %s\n", var, data), file = msFile)

    writeMeta("compound", getSIRIUSCmpName())
    writeMeta("parentmass", parentMZ)
    writeMeta("ionization", as.character(adduct, format = "sirius"))

    cat(">ms1peaks\n", file = msFile)
    write.table(plistMS[, c("mz", "intensity")], msFile, row.names = FALSE, col.names = FALSE)

    cat("\n>ms2peaks\n", file = msFile)
    write.table(plistMSMS[, c("mz", "intensity")], msFile, row.names = FALSE, col.names = FALSE)

    close(msFile)
}

unifySirNames <- function(sir)
{
    unNames <- c(# MonoisotopicMass = "neutralMass", UNDONE
                 smiles = "SMILES",
                 inchikey2D = "InChIKey1",
                 inchi = "InChI",
                 pubchemids = "identifier",
                 PubChemNumberPatents = "numberPatents",
                 score = "score",
                 molecularFormula = "neutral_formula",
                 xlogp = "XlogP",
                 name = "compoundName",
                 links = "libraryLinks",
                 
                 # some names were changed in 4.4 and new columns were added
                 # UNDONE: there is also a compound_identifications.csv file with slightly different columns, use that?
                 formulaRank = "formulaRank",
                 InChI = "InChI",
                 InChIkey2D = "InChIKey1",
                 "CSI:FingerIDScore" = "score",
                 TreeIsotopeScore = "SIR_formulaScore" # UNDONE: better name?
                 )

    unNames <- unNames[names(unNames) %in% names(sir)] # filter out missing
    setnames(sir, names(unNames), unNames)

    return(sir[, unNames, with = FALSE]) # filter out any other columns
}

SIRMPFinishHandler <- function(cmd)
{
    if (tools::file_ext(cmd$outPath) == "sirius")
    {
        # project directory was zipped, unzip to temp directory and process that instead
        uzpath <- tempfile()
        unzip(cmd$outPath, exdir = uzpath)
        cmd$outPath <- uzpath
    }
    
    pArgs <- list(adduct = cmd$adduct)
    if (!is.null(cmd[["processArgs"]]))
        pArgs <- c(pArgs, cmd$processArgs)
    res <- mapply(cmd$msFNames, cmd$MSMSPL, SIMPLIFY = FALSE,
                  FUN = function(n, m) do.call(cmd$processFunc, c(list(outPath = cmd$outPath, msFName = n, MSMS = m), pArgs)))
    return(res)
}

SIRMPPrepareHandler <- function(cmd)
{
    command <- patRoon:::getCommandWithOptPath(patRoon:::getSiriusBin(), "SIRIUS")
    
    inPath <- tempfile("sirius_in")
    outPath <- if (is.null(cmd[["projectPath"]])) tempfile("sirius_out") else cmd$projectPath
    # unlink(outPath, TRUE) # start with fresh output directory (otherwise previous results are combined)
    stopifnot(!file.exists(inPath) || !file.exists(outPath))
    dir.create(inPath)
    
    msFNames <- mapply(cmd$precMZs, cmd$MSPL, cmd$MSMSPL, cmd$resNames, FUN = function(pmz, mspl, msmspl, n)
    {
        ret <- file.path(inPath, paste0(n, ".ms"))
        patRoon:::makeSirMSFile(mspl, msmspl, pmz, cmd$adduct, ret)
        return(ret)
    })
    
    bArgs <- character()
    if (!cmd$dryRun)
        bArgs <- c("-i", inPath)
    bArgs <- c(bArgs, "-o", outPath, cmd$args)
    
    return(utils::modifyList(cmd, list(command = command, args = bArgs, outPath = outPath, msFNames = msFNames)))
}

runSIRIUS <- function(precursorMZs, MSPLists, MSMSPLists, resNames, profile, adducts, adductsChr, ppmMax, elements,
                      database, noise, cores, withFingerID, fingerIDDatabase, topMost, projectPath,
                      extraOptsGeneral, extraOptsFormula, verbose, processFunc, processArgs, splitBatches, dryRun)
{
    mainArgs <- character()
    if (!is.null(cores))
        mainArgs <- c("--cores", cores)
    if (!is.null(extraOptsGeneral))
        mainArgs <- c(mainArgs, extraOptsGeneral)
    if (dryRun)
        mainArgs <- c(mainArgs, "--no-project-check") # internal option, see https://github.com/boecker-lab/sirius/issues/42
    
    formArgs <- c("formula",
                  "-p", profile,
                  "-e", elements,
                  "--ppm-max", ppmMax,
                  "-c", topMost)
    
    if (!is.null(database))
        formArgs <- c(formArgs, "-d", database)
    if (!is.null(noise))
        formArgs <- c(formArgs, "-n", noise)
    if (!is.null(extraOptsFormula))
        formArgs <- c(formArgs, extraOptsFormula)

    isV5 <- isSIRIUS5() # UNDONE: what if SIRUS is only available on the workers?
    cmpArgs <- if (withFingerID && isV5)
        c("fingerprint", "structure", "--database", fingerIDDatabase)
    else if (withFingerID)
        c("structure", "--database", fingerIDDatabase)
    else
        character()
    
    batchn <- 1
    if (splitBatches) 
    {
        mpm <- getOption("patRoon.MP.method", "classic")
        batchn <- if (mpm == "classic") getOption("patRoon.MP.maxProcs") else future::nbrOfWorkers()
    }
    
    # perform SIRIUS batches per adduct: we cannot specify multiple adducts per run
    retAdduct <- sapply(unique(adductsChr), function(addChr)
    {
        doWhich <- which(adductsChr == addChr)
        batches <- splitInNBatches(doWhich, batchn)
        add <- adducts[[match(addChr, adductsChr)]]
        
        printf("Annotating %d features with adduct %s...\n", length(doWhich), addChr)
        
        cmdQueue <- lapply(seq_along(batches), function(bi)
        {
            batch <- batches[[bi]]
            fArgs <- c(formArgs, "-i", addChr)
            allArgs <- c(mainArgs, fArgs, cmpArgs)
            if (isV5)
                allArgs <- c(allArgs, "write-summaries")
            return(list(args = allArgs, precMZs = precursorMZs[batch], MSPL = MSPLists[batch],
                        MSMSPL = MSMSPLists[batch], adduct = add, projectPath = projectPath, resNames = resNames[batch],
                        processFunc = processFunc, processArgs = processArgs, dryRun = dryRun,
                        logFile = paste0("sirius-batch_", bi, "-", addChr, ".txt")))
        })
        
        singular <- length(cmdQueue) == 1
        ret <- executeMultiProcess(cmdQueue, finishHandler = SIRMPFinishHandler,
                                   prepareHandler = SIRMPPrepareHandler, printOutput = verbose && singular,
                                   printError = verbose && singular, showProgress = !singular,
                                   logSubDir = paste0("sirius_", if (withFingerID) "compounds" else "formulas"))
        
        return(setNames(unlist(ret, recursive = FALSE, use.names = FALSE), resNames[doWhich]))
    }, simplify = FALSE)
    
    results <- Reduce(modifyList, retAdduct) # remove adduct dimension
    results <- results[resNames] # and re-order
    return(results)
}

doSIRIUS <- function(fGroups, MSPeakLists, doFeatures, profile, adduct, relMzDev, elements,
                     database, noise, cores, withFingerID, fingerIDDatabase, topMost, projectPath,
                     extraOptsGeneral, extraOptsFormula, verbose, cacheName, processFunc, processArgs,
                     splitBatches, dryRun)
{
    if (length(MSPeakLists) == 0)
        return(list())
    
    # only do relevant feature groups
    MSPeakLists <- MSPeakLists[, intersect(names(fGroups), groupNames(MSPeakLists))]
    if (length(MSPeakLists) == 0)
        return(list())
    
    cacheDB <- openCacheDBScope() # open manually so caching code doesn't need to on each R/W access
    baseHash <- makeHash(profile, relMzDev, elements, database, noise,
                         withFingerID, fingerIDDatabase, topMost, extraOptsGeneral,
                         extraOptsFormula, processArgs)
    setHash <- makeHash(MSPeakLists, baseHash, doFeatures)
    cachedSet <- loadCacheSet(cacheName, setHash, cacheDB)
    
    if (doFeatures)
    {
        pLists <- peakLists(MSPeakLists)
        flattenedPLists <- unlist(pLists, recursive = FALSE)
        
        # important: assign before flattenedPLists subset steps below
        flPLMeta <- data.table(name = names(flattenedPLists),
                               group = unlist(lapply(pLists, names), use.names = FALSE),
                               analysis = rep(names(pLists), times = lengths(pLists)))
        
        # ensure only present features are done
        ftind <- groupFeatIndex(fGroups)
        flPLMeta <- flPLMeta[mapply(group, analysis, FUN = function(grp, ana)
        {
            anai <- match(ana, analyses(fGroups))
            return(!is.na(anai) && ftind[[grp]][anai] != 0)
        })]
        flattenedPLists <- flattenedPLists[flPLMeta$name]
    }
    else
    {
        flattenedPLists <- averagedPeakLists(MSPeakLists)
        flPLMeta <- data.table(name = names(flattenedPLists), group = names(flattenedPLists))
    }
    
    validPL <- function(pl) !is.null(pl[["MS"]]) && !is.null(pl[["MSMS"]]) && any(pl[["MS"]]$precursor)
    flattenedPLists <- flattenedPLists[sapply(flattenedPLists, validPL)]
    flPLMeta <- flPLMeta[name %in% names(flattenedPLists)]
    
    gNamesTBD <- unique(flPLMeta$group)
    fgAdd <- getFGroupAdducts(gNamesTBD, annotations(fGroups)[match(gNamesTBD, group)], adduct, "sirius")
    
    flPLMeta[, c("adduct", "adductChr") := .(fgAdd$grpAdducts[group], fgAdd$grpAdductsChr[group])]
    
    flPLMeta[, hash := mapply(flattenedPLists, flPLMeta$adductChr, FUN = makeHash, MoreArgs = list(baseHash))]
    if (is.null(cachedSet))
        saveCacheSet(cacheName, flPLMeta$hash, setHash, cacheDB)
    
    if (length(flattenedPLists) > 0)        
    {
        cachedResults <- pruneList(sapply(flPLMeta$hash, function(h)
        {
            res <- NULL
            if (!is.null(cachedSet))
                res <- cachedSet[[h]]
            if (is.null(res))
                res <- loadCacheData(cacheName, h, cacheDB)
            return(res)
        }, simplify = FALSE))
        
        flPLMeta[, cached := hash %in% names(cachedResults)]
        doWhich <- which(!flPLMeta$cached)
        
        if (length(doWhich) > 0)
        {
            doPLists <- flattenedPLists[doWhich]
            plmzs <- lapply(doPLists, function(pl) pl[["MS"]][precursor == TRUE, mz])
            mspls <- lapply(doPLists, "[[", "MS")
            msmspls <- lapply(doPLists, "[[", "MSMS")
            allResults <- runSIRIUS(plmzs, mspls, msmspls, flPLMeta$name[doWhich], profile, flPLMeta$adduct[doWhich],
                                    flPLMeta$adductChr[doWhich], relMzDev, elements, database, noise, cores,
                                    withFingerID, fingerIDDatabase, topMost, projectPath,
                                    extraOptsGeneral, extraOptsFormula, verbose, processFunc, processArgs, splitBatches,
                                    dryRun)
        }
        else
            allResults <- list()
        
        mergeCachedGroupResults <- function(meta, res)
        {
            if (length(cachedResults) > 0)
            {
                metaCached <- meta[cached == TRUE]
                res <- c(res, setNames(cachedResults[metaCached$hash], metaCached$group))
                res <- res[intersect(meta$group, names(res))] # ensure correct order
            }
            return(res)
        }
        
        if (doFeatures)
        {
            ret <- sapply(unique(flPLMeta$analysis), function(ana)
            {
                meta <- flPLMeta[analysis == ana]
                metaNotCached <- meta[cached == FALSE]
                res <- setNames(allResults[metaNotCached$name], metaNotCached$group)
                res <- mergeCachedGroupResults(meta, res)
                return(res)
            }, simplify = FALSE)
        }
        else
            ret <- mergeCachedGroupResults(flPLMeta, allResults)
        
        metaNotCached <- flPLMeta[cached == FALSE]
        for (i in seq_len(nrow(metaNotCached)))
            saveCacheData(cacheName, allResults[[metaNotCached$name[i]]], metaNotCached$hash[i], cacheDB)
    }
    else
        ret <- list()
    
    return(ret)
}
