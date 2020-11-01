#' @include main.R

getSiriusBin <- function()
{
    # UNDONE: check bin/ subdir?
    return(if (checkmate::testOS("windows")) "sirius" else "sirius.sh")
}

isSIRIUSPre44 <- function()
{
    # SIRIUS 4.4 returns a version string when running with --version, older
    # versions don't actually report version info...
    return(!any(grepl("^SIRIUS 4\\.", executeCommand(getCommandWithOptPath(getSiriusBin(), "SIRIUS"),
                                                     "--version", stdout = TRUE, stderr = FALSE))))
}

getSiriusResultPath <- function(outPath, msFName, cmpName, isPre44)
{
    # format is resultno_specname_compoundname, older versions start with 1, newer with 0
    msFName <- basename(tools::file_path_sans_ext(msFName))
    return(list.files(outPath, pattern = sprintf("[0-9]+_%s_%s", msFName, cmpName), full.names = TRUE))
}

getSiriusFragFiles <- function(resultPath, isPre44)
{
    if (isPre44)
        pat <- "[:0-9:]+_([A-Za-z0-9]+).*\\.ms"
    else
        pat <- "([A-Za-z0-9]+).*\\.tsv"
    return(list.files(file.path(resultPath, "spectra"), full.names = TRUE, pattern = pat))
}

getFormulaFromSiriusFragFile <- function(ffile, isPre44)
{
    if (isPre44)
        pat <- "[:0-9:]+_([A-Za-z0-9]+).*\\.ms"
    else
        pat <- "([A-Za-z0-9]+).*\\.tsv"
    return(gsub(pat, "\\1", basename(ffile)))
}

makeSirMSFile <- function(plistMS, plistMSMS, parentMZ, compound, ionization, out)
{
    msFile <- file(out, "w")

    writeMeta <- function(var, data) cat(sprintf(">%s %s\n", var, data), file = msFile)

    writeMeta("compound", compound)
    writeMeta("parentmass", parentMZ)
    writeMeta("ionization", ionization)

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
                 molecularFormula = "formula",
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

# get a command queue list that can be used with executeMultiProcess()
getSiriusCommand <- function(precursorMZ, MSPList, MSMSPList, profile, adduct, ppmMax, elements,
                             database, noise, withFingerID, fingerIDDatabase, topMost, extraOpts,
                             isPre44)
{
    outPath <- tempfile("sirius")
    # unlink(outPath, TRUE) # start with fresh output directory (otherwise previous results are combined)

    stopifnot(!file.exists(outPath))

    msFName <- tempfile("spec", fileext = ".ms")    
    ionization <- as.character(adduct, format = "sirius")
    mainArgs <- c("-p", profile,
                  "-e", elements,
                  "--ppm-max", ppmMax,
                  "-c", topMost)

    if (!is.null(database))
        mainArgs <- c(mainArgs, "-d", database)
    if (!is.null(noise))
        mainArgs <- c(mainArgs, "-n", noise)
    if (!is.null(extraOpts))
        mainArgs <- c(mainArgs, extraOpts)

    if (isPre44)
    {
        if (withFingerID)
            mainArgs <- c(mainArgs, "--fingerid", "--fingerid-db", fingerIDDatabase)
        args <- c(mainArgs, "-o", outPath, msFName)
    }
    else
    {
        args <- c("-o", outPath, "-i", msFName, "formula", mainArgs)
        if (withFingerID)
            args <- c(args, "structure", "--database", fingerIDDatabase)
    }

    cmpName <- "unknownCompound"
    makeSirMSFile(MSPList, MSMSPList, precursorMZ, cmpName, ionization, msFName)

    return(list(command = getCommandWithOptPath(getSiriusBin(), "SIRIUS"), args = args,
                outPath = outPath, msFName = msFName, cmpName = cmpName, isPre44 = isPre44))
}

runSIRIUS <- function(precursorMZs, MSPLists, MSMSPLists, resNames, profile, adduct, ppmMax, elements,
                      database, noise, cores, withFingerID, fingerIDDatabase, topMost,
                      extraOptsGeneral, extraOptsFormula, verbose, isPre44,
                      processFunc, processArgs, SIRBatchSize)
{
    ionization <- as.character(adduct, format = "sirius")
    cmpName <- "unknownCompound"
    
    mainArgs <- character()
    if (!is.null(cores))
        mainArgs <- c("--cores", cores)
    if (!is.null(extraOptsGeneral))
        mainArgs <- c(mainArgs, extraOptsGeneral)
    
    formArgs <- c("-p", profile,
                  "-e", elements,
                  "--ppm-max", ppmMax,
                  "-c", topMost,
                  "-i", ionization)
    
    if (!is.null(database))
        formArgs <- c(formArgs, "-d", database)
    if (!is.null(noise))
        formArgs <- c(formArgs, "-n", noise)
    if (!is.null(extraOptsFormula))
        formArgs <- c(formArgs, extraOptsFormula)
    
    if (isPre44)
    {
        if (withFingerID)
            formArgs <- c(formArgs, "--fingerid", "--fingerid-db", fingerIDDatabase)
        args <- c(mainArgs, formArgs)
    }
    else
    {
        args <- c(mainArgs, "formula", formArgs)
        if (withFingerID)
            args <- c(args, "structure", "--database", fingerIDDatabase)
    }
    
    if (SIRBatchSize == 0)
        batches <- list(seq_along(precursorMZs))
    else
        batches <- splitInBatches(seq_along(precursorMZs), SIRBatchSize)
    
    command <- getCommandWithOptPath(getSiriusBin(), "SIRIUS")
    
    cmdQueue <- lapply(seq_along(batches), function(bi)
    {
        batch <- batches[[bi]]
        return(list(command = command, precMZs = precursorMZs[batch], MSPL = MSPLists[batch], MSMSPL = MSMSPLists[batch],
                    logFile = paste0("sirius-batch_", bi, ".txt")))
    })
    
    singular <- length(cmdQueue) == 1
    ret <- executeMultiProcess(cmdQueue, finishHandler = function(cmd)
    {
        pArgs <- list(cmpName = cmpName, adduct = adduct, isPre44 = isPre44)
        if (!is.null(processArgs))
            pArgs <- c(pArgs, processArgs)
        
        res <- mapply(cmd$msFNames, cmd$MSMSPL, SIMPLIFY = FALSE,
                      FUN = function(n, m) do.call(processFunc, c(list(outPath = cmd$outPath, msFName = n, MSMS = m), pArgs)))
        return(res)
    }, prepareHandler = function(cmd)
    {
        inPath <- tempfile("sirius_in")
        outPath <- tempfile("sirius_out")
        # unlink(outPath, TRUE) # start with fresh output directory (otherwise previous results are combined)
        stopifnot(!file.exists(inPath) || !file.exists(outPath))
        dir.create(inPath)
        
        msFNames <- mapply(cmd$precMZs, cmd$MSPL, cmd$MSMSPL, FUN = function(pmz, mspl, msmspl)
        {
            ret <- tempfile("spec", fileext = ".ms", tmpdir = inPath)
            makeSirMSFile(mspl, msmspl, pmz, cmpName, ionization, ret)
            return(ret)
        })
        
        bArgs <- if (isPre44) c(args, "-o", outPath, inPath) else c("-i", inPath, "-o", outPath, args)
        
        return(c(cmd, list(args = bArgs, outPath = outPath, msFNames = msFNames)))
    }, printOutput = verbose && singular, printError = verbose && singular,
    showProgress = !singular, logSubDir = paste0("sirius_", if (withFingerID) "compounds" else "formulas"))

    return(setNames(unlist(ret, recursive = FALSE, use.names = FALSE), resNames))
}

doSIRIUS <- function(fGroups, MSPeakLists, doFeatures, profile, adduct, relMzDev, elements,
                     database, noise, cores, withFingerID, fingerIDDatabase, topMost,
                     extraOptsGeneral, extraOptsFormula, verbose, cacheName, processFunc, processArgs,
                     SIRBatchSize)
{
    isPre44 <- isSIRIUSPre44()
    gNames <- names(fGroups)
    
    # only do relevant feature groups
    MSPeakLists <- MSPeakLists[, intersect(gNames, groupNames(MSPeakLists))]
    
    if (length(MSPeakLists) == 0)
        return(list())
    
    cacheDB <- openCacheDBScope() # open manually so caching code doesn't need to on each R/W access
    baseHash <- makeHash(profile, adduct, relMzDev, elements, database, noise,
                         withFingerID, fingerIDDatabase, topMost, extraOptsGeneral,
                         extraOptsFormula, isPre44, processArgs)
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
    
    flPLMeta[, hash := sapply(flattenedPLists, makeHash, baseHash)]
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
        doPLists <- flattenedPLists[!flPLMeta$cached]
        
        if (length(doPLists) > 0)
        {
            plmzs <- lapply(doPLists, function(pl) pl[["MS"]][precursor == TRUE, mz])
            mspls <- lapply(doPLists, "[[", "MS")
            msmspls <- lapply(doPLists, "[[", "MSMS")
            
            allResults <- runSIRIUS(plmzs, mspls, msmspls, flPLMeta[cached == FALSE]$name, profile, adduct,
                                    relMzDev, elements, database, noise, cores, withFingerID, fingerIDDatabase, topMost,
                                    extraOptsGeneral, extraOptsFormula, verbose, isPre44, processFunc, processArgs, SIRBatchSize)
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
