# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R

getSiriusBin <- function()
{
    # NOTE: this seems to fluctuate every other SIRIUS version...
    return("sirius")
}

isSIRIUS5 <- function()
{
    out <- executeCommand(patRoon:::getExtDepPath("sirius"), "--version", stdout = TRUE)
    return(any(grepl("^(SIRIUS 5\\.)", out)))
}

getSIRIUSCmpName <- function() "unknownCompound"

getSiriusResultPath <- function(outPath, msFName)
{
    # format is resultno_specname_compoundname, older versions start with 1, newer with 0
    msFName <- basename(tools::file_path_sans_ext(msFName))
    return(list.files(outPath, pattern = sprintf("[0-9]+_%s_%s", msFName, getSIRIUSCmpName()), full.names = TRUE))
}

getAndPrepareSIRIUSResFiles <- function(resultPath, subDir, ext)
{
    spPath <- file.path(resultPath, subDir)
    if (file.exists(spPath) && !file.info(spPath, extra_cols = FALSE)$isdir)
    {
        # NOTE: SIRIUS 5 packs most result files --> unzip them
        exDir <- paste0(spPath, "-unz")
        unzip(spPath, exdir = exDir)
        spPath <- exDir
    }
    
    pat <- paste0("([A-Za-z0-9]+).*\\.", ext)
    return(list.files(spPath, full.names = TRUE, pattern = pat))
}

getFormulaFromSIRIUSResFile <- function(ffile, ext)
{
    pat <- paste0("([A-Za-z0-9]+).*\\.", ext)
    return(gsub(pat, "\\1", basename(ffile)))
}

loadSIRIUSFingerprints <- function(resultPath, formulas, adduct)
{
    fingerprints <- data.table()
    if (file.exists(file.path(resultPath, "fingerprints")))
    {
        fpFiles <- getAndPrepareSIRIUSResFiles(resultPath, "fingerprints", "fpt")
        fpForms <- getFormulaFromSIRIUSResFile(fpFiles, "fpt")
        # obtain neutral formula of FP results from candidate list
        formCandidates <- fread(file.path(resultPath, "formula_candidates.tsv"))
        fpForms <- formCandidates[match(fpForms, precursorFormula)]$molecularFormula
        for (i in seq_along(fpFiles))
        {
            if (fpForms[i] %chin% formulas) # only consider FPs of candidates
                fingerprints[, (fpForms[i]) := fread(fpFiles[i])]
        }
        # add absoluteIndices
        fpMD <- file.path(resultPath, "..", if (adduct@charge > 0) "csi_fingerid.tsv" else "csi_fingerid_neg.tsv")
        fingerprints[, absoluteIndex := fread(file.path(fpMD), select = "absoluteIndex")[[1]]]
    }
    return(fingerprints)
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

hasSIRIUSLogin <- function()
{
    out <- system2(getExtDepPath("sirius"), c("login", "--show"), stdout = TRUE, stderr = FALSE)
    notLoggedIn <- any(out == "Not logged in.")
    isLoggedIn <- any(grepl("^Logged in as:", out))
    if (notLoggedIn == isLoggedIn)
    {
        warning("Could not determine if SIRIUS is currently logged in", call. = FALSE)
        return(FALSE)
    }
    return(isLoggedIn)
}

doSIRIUSLogin <- function(login, force)
{
    if (isFALSE(login))
        return(invisible(NULL)) # no need to do anything
    
    if (force || !hasSIRIUSLogin())
    {
        if (length(login) == 1 && login == "check")
            stop("There is no active SIRIUS login. Please consult the SIRIUS documentation and patRoon handbook for details.")
        
        if (length(login) == 1 && login == "interactive")
        {
            if (!interactive())
                stop("Cannot perform interactive login in non-interactive R sessions!", call. = FALSE)
            
            # NOTE: if "username" is part of the prompt then the RStudio backend of getPass won't hide the text input
            login <- c(username = getPass::getPass("Please enter your SIRIUS username", noblank = TRUE),
                       password = getPass::getPass("Please enter your SIRIUS password", noblank = TRUE))
        }
        
        if (!"username" %in% names(login) || !nzchar(login["username"]))
            stop("Please provide the username of your SIRIUS account", call. = FALSE)
        if (!"password" %in% names(login) || !nzchar(login["password"]))
            stop("Please provide the password of your SIRIUS account", call. = FALSE)

        # NOTE: processx::run() is used as it allows correctly setting the environment, which doesn't seem to work very well
        # with base::system2()
        runv <- processx::run(getExtDepPath("sirius"), c("login", "--user-env=SIRUSER", "--password-env=SIRPW"),
                              env = c("current", SIRUSER = login[["username"]], SIRPW = login[["password"]]))

        if ((!is.na(runv$status) && runv$status != 0) || !grepl("Login successful!", runv$stdout))
        {
            cat(runv$stderr)
            stop("Failed to perform a SIRIUS login! See error output above for details.", call. = FALSE)
        }
    }
    invisible(NULL)
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
    command <- patRoon:::getExtDepPath("sirius")
    
    # UNDONE: it seems we would only need to log in once per worker, is this adding a lot of overhead?
    doSIRIUSLogin(cmd$login, cmd$alwaysLogin)

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
                      database, noise, cores, withFingerID, fingerIDDatabase, topMost, projectPath, login, alwaysLogin,
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
    cmpArgs <- character()
    if (withFingerID != "none")
    {
        if (withFingerID == "fingerprint")
        {
            if (!isV5)
                stop("Can only obtain fingerprints with SIRIUS5", call. = FALSE)    
            cmpArgs <- "fingerprint"
        }
        else
        {
            cmpArgs <- if (isV5)
                c("fingerprint", "structure", "--database", fingerIDDatabase)
            else
                c("structure", "--database", fingerIDDatabase)
        }
    }
    
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
                        MSMSPL = MSMSPLists[batch], adduct = add, projectPath = projectPath, login = login,
                        alwaysLogin = alwaysLogin, verbose = verbose, resNames = resNames[batch],
                        processFunc = processFunc, processArgs = processArgs, dryRun = dryRun,
                        logFile = paste0("sirius-batch_", bi, "-", addChr, ".txt")))
        })
        
        singular <- length(cmdQueue) == 1
        ret <- executeMultiProcess(cmdQueue, finishHandler = SIRMPFinishHandler,
                                   prepareHandler = SIRMPPrepareHandler, printOutput = verbose && singular,
                                   printError = verbose && singular, showProgress = !singular,
                                   logSubDir = paste0("sirius_", if (withFingerID == "structure") "compounds" else "formulas"))
        
        return(setNames(unlist(ret, recursive = FALSE, use.names = FALSE), resNames[doWhich]))
    }, simplify = FALSE)
    
    results <- Reduce(modifyList, retAdduct) # remove adduct dimension
    results <- results[resNames] # and re-order
    return(results)
}

doSIRIUS <- function(fGroups, MSPeakLists, doFeatures, profile, adduct, relMzDev, elements,
                     database, noise, cores, withFingerID, fingerIDDatabase, topMost, projectPath, login, alwaysLogin,
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
                                    withFingerID, fingerIDDatabase, topMost, projectPath, login, alwaysLogin,
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

getMS2QTFPs <- function(featAnnSIR)
{
    # get all SIRIUS fingerprints in a table that is compatible with MS2Quant/MS2Tox
    
    FPs <- featAnnSIR@fingerprints[sapply(featAnnSIR@fingerprints, nrow) > 0]
    if (length(FPs) == 0)
        return(data.table())
    
    allFPs <- rbindlist(lapply(FPs, transpose, keep.names = "neutral_formula", make.names = "absoluteIndex"),
                        idcol = "group")
    fpColRange <- seq(3, ncol(allFPs))
    setnames(allFPs, fpColRange, paste0("Un", names(allFPs)[fpColRange]))
    allFPs[, id := as.character(seq_len(nrow(allFPs)))]
    allFPs[, ionization := mapply(group, neutral_formula, FUN = function(g, f)
    {
        # UNDONE: need to check for empty fragInfo?
        return(featAnnSIR[[g]][neutral_formula == f]$fragInfo[[1]]$ionization[1])
    })]
    
    allFPs[, predion := paste0(neutral_formula, "_", ionization)]
    
    # dummies for MS2Tox
    allFPs[, foldernumber := 0]
    allFPs[, predform := ""]
    
    return(allFPs[])
}

predictRespFactorsSIRFPs <- function(featAnnSIR, gInfo, calibrants, eluent, organicModifier, pHAq, concUnit)
{
    featAnnSIR <- featAnnSIR[rownames(gInfo)]
    
    allFPs <- getMS2QTFPs(featAnnSIR)
    
    if (any(!allFPs$ionization %chin% c("[M+H]+", "[M]+")))
    {
        warning("One or more features are with adducts other than [M+H]+/[M]+. ",
                "These are not (yet) supported by MS2Quant and will be ignored.", call. = FALSE)
        allFPs <- allFPs[ionization %chin% c("[M+H]+", "[M]+")]
    }
    
    if (nrow(allFPs) == 0)
        return(list(RFs = data.table(), MD = list()))
    
    # NOTE: we set the area to one to easily get the response factor below
    unknowns <- data.table(identifier = allFPs$id, retention_time = gInfo[allFPs$group, "rts"],
                           SMILES = NA_character_, conc_M = NA_real_, area = 1)
    
    baseHash <- makeHash(calibrants, eluent, organicModifier, pHAq)
    splFPs <- split(allFPs[, -c("id", "group")], seq_len(nrow(allFPs)))
    hashes <- sapply(seq_len(nrow(allFPs)), function(i) makeHash(baseHash, splFPs[[i]], unknowns$retention_time[i]))
    
    cachedData <- loadCacheData("RF_SIRFP", hashes, simplify = FALSE)
    indsTODO <- if (!is.null(cachedData)) which(!hashes %in% names(cachedData)) else seq_along(hashes)
    hashesTODO <- hashes[indsTODO]
    
    MS2QRes <- NULL
    if (length(indsTODO) > 0)
    {
        MS2QRes <- getMS2QuantRes(calibrants, unknowns[indsTODO], eluent, organicModifier, pHAq, allFPs[indsTODO])
        saveCacheData("MS2QMD", MS2QRes$MD, baseHash)
        MS2QRes$RFs <- MS2QRes$RFs[, c("identifier", "RF_pred"), with = FALSE]
        setnames(MS2QRes$RFs, "RF_pred", "RF_SIRFP")
        for (i in seq_len(nrow(MS2QRes$RFs)))
            saveCacheData("RF_SIRFP", MS2QRes$RFs$RF_SIRFP[i], hashesTODO[i])
        
        MS2QRes$RFs <- merge(allFPs[, c("group", "neutral_formula", "id"), with = FALSE], MS2QRes$RFs, by.x = "id",
                             by.y = "identifier", sort = FALSE)
    }
    
    if (!is.null(cachedData))
    {
        cachedRFs <- rbindlist(lapply(cachedData, function(cd) data.table(RF_SIRFP = cd)), idcol = "hash")
        cachedRFs[, neutral_formula := allFPs$neutral_formula[match(hash, hashes)]]
        cachedRFs[, group := allFPs$group[match(hash, hashes)]]
        cachedRFs[, hash := NULL]
        
        if (is.null(MS2QRes))
        {
            MD <- loadCacheData("MS2QMD", baseHash)
            if (is.null(MD))
            {
                warning("Could not find cached calibration data! You may have an old cache file. ",
                        "Please clear any cached data, eg by running: clearCache(\"RF_SMILES\")", call. = FALSE)
                MD <- list()
            }
            MS2QRes <- list(RFs = cachedRFs, MD = MD)
        }
        else
        {
            MS2QRes$RFs <- rbind(MS2QRes$RFs, cachedRFs)
            # restore order
            boundHashes <- c(hashesTODO, names(cachedData))
            MS2QRes$RFs <- MS2QRes$RFs[match(hashes, boundHashes)]
        }
    }
    
    if (!is.null(MS2QRes$RFs[["id"]]))
        MS2QRes$RFs[, id := NULL]

    # NOTE: do unit conversion the last thing, so we can still use cached data if the user merely changed the unit
    # NOTE: need to take the inverse before conversion
    MS2QRes$RFs[, RF_SIRFP := 1/convertConc(1/RF_SIRFP[1], "M", concUnit, formulaMW(neutral_formula[1])),
                by = "neutral_formula"][]
    
    return(MS2QRes)
}

predictLC50SIRFPs <- function(featAnnSIR, LC50Mode, concUnit)
{
    # UNDONE: check supported adducts?
    
    allFPs <- getMS2QTFPs(featAnnSIR)
    if (nrow(allFPs) == 0)
        return(data.table())

    baseHash <- makeHash(LC50Mode)
    hashes <- sapply(split(allFPs[, -c("id", "group")], seq_len(nrow(allFPs))), function(s) makeHash(baseHash, s))
    
    cachedData <- loadCacheData("LC50_SIRFP", hashes, simplify = FALSE)
    indsTODO <- if (!is.null(cachedData)) which(!hashes %in% names(cachedData)) else seq_along(hashes)
    hashesTODO <- hashes[indsTODO]
    
    LC50s <- NULL
    if (length(indsTODO) > 0)
    {
        allFPsTODO <- allFPs[indsTODO]
        allFPsTODO[, exactMass := sapply(neutral_formula, getFormulaMass)]
        suppressMessages(utils::capture.output(LC50s <- MS2Tox::FishLC50Prediction(allFPsTODO, LC50Mode)))
    
        LC50s <- merge(allFPsTODO[, c("group", "neutral_formula", "id"), with = FALSE],
                       LC50s[, c("id", "LC50_predicted")], by = "id", sort = FALSE)
        setnames(LC50s, "LC50_predicted", "LC50_SIRFP")
        LC50s[, id := NULL]
        
        for (i in seq_len(nrow(LC50s)))
            saveCacheData("LC50_SIRFP", LC50s$LC50_SIRFP[i], hashesTODO[i])
    }
    
    if (!is.null(cachedData))
    {
        cachedLC50s <- rbindlist(lapply(cachedData, function(cd) data.table(LC50_SIRFP = cd)), idcol = "hash")
        cachedLC50s[, neutral_formula := allFPs$neutral_formula[match(hash, hashes)]]
        cachedLC50s[, group := allFPs$group[match(hash, hashes)]]
        cachedLC50s[, hash := NULL]
        
        if (is.null(LC50s))
            LC50s <- cachedLC50s
        else
        {
            LC50s <- rbind(LC50s, cachedLC50s)
            # restore order
            boundHashes <- c(hashesTODO, names(cachedData))
            LC50s <- LC50s[match(hashes, boundHashes)]
        }
    }

    # NOTE: do unit conversion the last thing, so we can still use cached data if the user merely changed the unit
    # NOTE: need to take the inverse before conversion
    LC50s[, LC50_SIRFP := convertConc(LC50_SIRFP[1], "log mM", concUnit, formulaMW(neutral_formula[1])),
        by = "neutral_formula"]

    return(LC50s[])
}

syncSIRFPs <- function(obj)
{
    # sync fingerprints
    obj@fingerprints <- obj@fingerprints[names(obj@fingerprints) %chin% groupNames(obj)]
    obj@fingerprints <- pruneList(Map(obj@fingerprints, obj@groupAnnotations[names(obj@fingerprints)], f = function(fp, ann)
    {
        fpForms <- intersect(names(fp), ann$neutral_formula)
        fp <- fp[, c(fpForms, "absoluteIndex"), with = FALSE]
        return(if (ncol(fp) == 1) NULL else fp) # nullify if no candidates left
    }))
    return(obj)
}
