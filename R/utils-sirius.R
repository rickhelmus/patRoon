#' @include main.R

getSiriusBin <- function()
{
    # UNDONE: check bin/ subdir?
    return(if (checkmate::testOS("windows")) "sirius" else "sirius.sh")
}

getSIRIUSCmpName <- function() "unknownCompound"

getSiriusResultPath <- function(outPath, msFName)
{
    # format is resultno_specname_compoundname, older versions start with 1, newer with 0
    msFName <- basename(tools::file_path_sans_ext(msFName))
    return(list.files(outPath, pattern = sprintf("[0-9]+_%s_%s", msFName, getSIRIUSCmpName()), full.names = TRUE))
}

getSiriusFragFiles <- function(resultPath)
{
    pat <- "([A-Za-z0-9]+).*\\.tsv"
    return(list.files(file.path(resultPath, "spectra"), full.names = TRUE, pattern = pat))
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

SIRMPFinishHandler <- function(cmd)
{
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
    outPath <- tempfile("sirius_out")
    # unlink(outPath, TRUE) # start with fresh output directory (otherwise previous results are combined)
    stopifnot(!file.exists(inPath) || !file.exists(outPath))
    dir.create(inPath)
    
    msFNames <- mapply(cmd$precMZs, cmd$MSPL, cmd$MSMSPL, FUN = function(pmz, mspl, msmspl)
    {
        ret <- tempfile("spec", fileext = ".ms", tmpdir = inPath)
        patRoon:::makeSirMSFile(mspl, msmspl, pmz, cmd$adduct, ret)
        return(ret)
    })
    
    bArgs <- c("-i", inPath, "-o", outPath, cmd$args)
    return(utils::modifyList(cmd, list(command = command, args = bArgs, outPath = outPath, msFNames = msFNames)))
}

runSIRIUS <- function(precursorMZs, MSPLists, MSMSPLists, resNames, profile, adduct, ppmMax, elements,
                      database, noise, cores, withFingerID, fingerIDDatabase, topMost,
                      extraOptsGeneral, extraOptsFormula, verbose,
                      processFunc, processArgs, splitBatches)
{
    ionization <- as.character(adduct, format = "sirius")
    
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
    
    args <- c(mainArgs, "formula", formArgs)
    if (withFingerID)
        args <- c(args, "structure", "--database", fingerIDDatabase)

    batchn <- 1
    if (splitBatches) 
    {
        mpm <- getOption("patRoon.multiProcMethod", "classic")
        batchn <- if (mpm == "classic") getOption("patRoon.maxProcAmount") else future::nbrOfWorkers()
    }
    batches <- splitInNBatches(seq_along(precursorMZs), batchn)
    
    cmdQueue <- lapply(seq_along(batches), function(bi)
    {
        batch <- batches[[bi]]
        return(list(args = args, precMZs = precursorMZs[batch], MSPL = MSPLists[batch], MSMSPL = MSMSPLists[batch],
                    adduct = adduct, processFunc = processFunc, processArgs = processArgs,
                    logFile = paste0("sirius-batch_", bi, ".txt")))
    })
    
    singular <- length(cmdQueue) == 1
    ret <- executeMultiProcess(cmdQueue, finishHandler = SIRMPFinishHandler,
                               prepareHandler = SIRMPPrepareHandler, printOutput = verbose && singular,
                               printError = verbose && singular, showProgress = !singular,
                               logSubDir = paste0("sirius_", if (withFingerID) "compounds" else "formulas"))

    return(setNames(unlist(ret, recursive = FALSE, use.names = FALSE), resNames))
}

doSIRIUS <- function(fGroups, MSPeakLists, doFeatures, profile, adduct, relMzDev, elements,
                     database, noise, cores, withFingerID, fingerIDDatabase, topMost,
                     extraOptsGeneral, extraOptsFormula, verbose, cacheName, processFunc, processArgs,
                     splitBatches)
{
    # only do relevant feature groups
    MSPeakLists <- MSPeakLists[, intersect(gNames, groupNames(MSPeakLists))]
    
    if (length(MSPeakLists) == 0)
        return(list())
    
    cacheDB <- openCacheDBScope() # open manually so caching code doesn't need to on each R/W access
    baseHash <- makeHash(profile, adduct, relMzDev, elements, database, noise,
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
                                    extraOptsGeneral, extraOptsFormula, verbose, processFunc, processArgs, splitBatches)
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
