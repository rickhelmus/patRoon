#' @include main.R

getSiriusBin <- function()
{
    si <- Sys.info()

    if (si[["sysname"]] == "Linux")
        return("sirius")

    # windows
    if (si[["machine"]] == "x86-64")
        return("sirius-console-64")

    return("sirius-console-32")
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
        pat <- "([A-Za-z0-9]+).*\\.csv"
    return(list.files(file.path(resultPath, "spectra"), full.names = TRUE, pattern = pat))
}

getFormulaFromSiriusFragFile <- function(ffile, isPre44)
{
    if (isPre44)
        pat <- "[:0-9:]+_([A-Za-z0-9]+).*\\.ms"
    else
        pat <- "([A-Za-z0-9]+).*\\.csv"
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
                 # UNDONE: there is also a compound_dentifications.csv file with slightly different columns, use that?
                 formulaRank = "formulaRank",
                 InChI = "InChI",
                 InChIkey2D = "InChIKey1",
                 "CSI:FingerID_Score" = "score",
                 TreeIsotope_Score = "SIR_formulaScore" # UNDONE: better name?
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

runSIRIUS <- function(precursorMZs, MSPLists, MSMSPLists, profile, adduct, ppmMax, elements,
                      database, noise, cores, withFingerID, fingerIDDatabase, topMost,
                      extraOptsGeneral, extraOptsFormula, verbose, isPre44)
{
    inPath <- tempfile("sirius_in")
    outPath <- tempfile("sirius_out")
    # unlink(outPath, TRUE) # start with fresh output directory (otherwise previous results are combined)
    stopifnot(!file.exists(outPath) && !file.exists(outPath))
    dir.create(inPath)
    
    ionization <- as.character(adduct, format = "sirius")
    cmpName <- "unknownCompound"
    
    msFNames <- mapply(precursorMZs, MSPLists, MSMSPLists, FUN = function(pmz, mspl, msmspl)
    {
        ret <- tempfile("spec", fileext = ".ms", tmpdir = inPath)
        makeSirMSFile(mspl, msmspl, pmz, cmpName, ionization, ret)
        return(ret)
    })

    mainArgs <- character()
    if (!is.null(cores))
        mainArgs <- c("--cores", cores)
    if (!is.null(extraOptsGeneral))
        mainArgs <- c(mainArgs, extraOptsGeneral)
    
    formArgs <- c("-p", profile,
                  "-e", elements,
                  "--ppm-max", ppmMax,
                  "-c", topMost)
    
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
        args <- c(mainArgs, formArgs, "-o", outPath, inPath)
    }
    else
    {
        args <- c(mainArgs, "-o", outPath, "-i", inPath, "formula", formArgs)
        if (withFingerID)
            args <- c(args, "structure", "--database", fingerIDDatabase)
    }
    
    verbose <- if (verbose) "" else FALSE
    executeCommand(getCommandWithOptPath(getSiriusBin(), "SIRIUS"), args, stdout = verbose, stderr = verbose)
    return(list(outPath = outPath, msFNames = msFNames, cmpName = cmpName))
}

doSIRIUS <- function(allGNames, featMZs, groupPeakLists, profile, adduct, relMzDev, elements,
                     database, noise, cores, withFingerID, fingerIDDatabase, topMost,
                     extraOptsGeneral, extraOptsFormula, verbose, isPre44, cacheName, cacheDB, processFunc)
{
    doFGroups <- names(featMZs)
    
    # skip fgroups without proper ms peak lists
    if (length(doFGroups) > 0)
        doFGroups <- doFGroups[sapply(doFGroups,
                                      function(grp) !is.null(groupPeakLists[[grp]][["MS"]]) &&
                                          any(groupPeakLists[[grp]][["MS"]]$precursor) &&
                                          !is.null(groupPeakLists[[grp]][["MSMS"]]))]
    
    featMZs <- featMZs[doFGroups]; groupPeakLists <- groupPeakLists[doFGroups]
    
    baseHash <- makeHash(profile, adduct, relMzDev, elements, database, noise, topMost, extraOptsGeneral, extraOptsFormula)
    setHash <- makeHash(featMZs, groupPeakLists, baseHash)
    cachedSet <- loadCacheSet(cacheName, setHash, cacheDB)
    
    hashes <- sapply(doFGroups, function(grp) makeHash(featMZs[grp], groupPeakLists[[grp]], baseHash))
    
    if (is.null(cachedSet))
        saveCacheSet(cacheName, hashes, setHash, cacheDB)
    
    cachedResults <- pruneList(sapply(hashes, function(h)
    {
        res <- NULL
        if (!is.null(cachedSet))
            res <- cachedSet[[h]]
        if (is.null(res))
            res <- loadCacheData(cacheName, h, cacheDB)
        return(res)
    }, simplify = FALSE))
    doFGroups <- setdiff(doFGroups, names(cachedResults))
    
    if (length(doFGroups) > 0)        
    {        
        plmzs <- lapply(doFGroups, function(g) groupPeakLists[[g]][["MS"]][precursor == TRUE, mz])
        mspls <- lapply(groupPeakLists[doFGroups], "[[", "MS")
        msmspls <- lapply(groupPeakLists[doFGroups], "[[", "MSMS")
        runData <- runSIRIUS(plmzs, mspls, msmspls, profile, adduct, relMzDev, elements,
                             database, noise, cores, withFingerID, fingerIDDatabase, topMost,
                             extraOptsGeneral, extraOptsFormula, verbose, isPre44)
        ret <- processFunc(doFGroups, adduct, hashes, runData, isPre44, cacheDB)
        ngrp <- length(ret)
    }
    else
    {
        ret <- list()
        ngrp <- 0
    }
    
    if (length(cachedResults) > 0)
    {
        ngrp <- ngrp + length(cachedResults)
        ret <- c(ret, cachedResults)
        ret <- ret[intersect(allGNames, names(ret))] # re-order
    }
    
    # prune after combining with cached results: these may also contain zero row results
    ret <- pruneList(ret, checkZeroRows = TRUE)
    
    return(ret)
}

doSIRIUS2 <- function(allGNames, MSPeakLists, doFeatures, profile, adduct, relMzDev, elements,
                      database, noise, cores, withFingerID, fingerIDDatabase, topMost,
                      extraOptsGeneral, extraOptsFormula, verbose, cacheName, processFunc)
{
    isPre44 <- isSIRIUSPre44()
    
    # only do relevant feature groups
    MSPeakLists <- MSPeakLists[, intersect(allGNames, groupNames(MSPeakLists))]
    
    cacheDB <- openCacheDBScope() # open manually so caching code doesn't need to on each R/W access
    baseHash <- makeHash(profile, adduct, relMzDev, elements, database, noise,
                         withFingerID, fingerIDDatabase, topMost, extraOptsGeneral,
                         extraOptsFormula, isPre44)
    setHash <- makeHash(MSPeakLists, baseHash, doFeatures)
    cachedSet <- loadCacheSet(cacheName, setHash, cacheDB)

    if (doFeatures)
    {
        pLists <- peakLists(MSPeakLists)
        flattenedPLists <- unlist(pLists, recursive = FALSE)
        
        # important: assign before subset step below
        flPLMeta <- data.table(name = names(flattenedPLists),
                               group = unlist(lapply(pLists, names), use.names = FALSE),
                               analysis = rep(names(pLists), times = lengths(pLists)))
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
        
        plmzs <- lapply(doPLists, function(pl) pl[["MS"]][precursor == TRUE, mz])
        mspls <- lapply(doPLists, "[[", "MS")
        msmspls <- lapply(doPLists, "[[", "MSMS")
        
        flPLMeta[, msFName := character()]
        if (length(doPLists) > 0)
        {
            runData <- runSIRIUS(plmzs, mspls, msmspls, profile, adduct, relMzDev, elements,
                                 database, noise, cores, withFingerID, fingerIDDatabase, topMost,
                                 extraOptsGeneral, extraOptsFormula, verbose, isPre44)
            flPLMeta[cached == FALSE, msFName := runData$msFNames]
        }
        else
            runData <- list()

        processResultSet <- function(meta)
        {
            if (any(!meta$cached))
            {
                procArgs <- list(outPath = runData$outPath, cmpName = runData$cmpName, adduct = adduct,
                                 isPre44 = isPre44, cacheDB = cacheDB)
                
                metaNew <- meta[cached == FALSE]
                res <- mapply(metaNew$msFName, metaNew$hash, SIMPLIFY = FALSE,
                              FUN = function(n, h) do.call(processFunc, c(list(msFName = n, hash = h),
                                                                          procArgs)))
                names(res) <- metaNew$group
            }
            else
                res <- list()
            
            if (length(cachedResults) > 0)
            {
                metaCached <- meta[cached == TRUE]
                res <- c(res, setNames(cachedResults[metaCached$hash], metaCached$group))
            }
            
            res <- pruneList(res, checkZeroRows = TRUE)
            res <- res[intersect(allGNames, names(res))] # ensure correct order
            
            return(res)
        }

        if (doFeatures)
            ret <- sapply(unique(flPLMeta$analysis), function(ana) processResultSet(flPLMeta[analysis == ana]),
                          simplify = FALSE)
        else
            ret <- processResultSet(flPLMeta)
    }
    else
        ret <- list()
    
    return(ret)
}
