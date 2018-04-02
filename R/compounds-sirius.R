#' @include main.R
#' @include compounds.R
#' @include utils-sirius.R
NULL

processSiriusCompounds <- function(cmd, exitStatus, retries)
{
    # format is resultno_specname_compoundname
    resultPath <- file.path(cmd$outPath, sprintf("1_%s_%s", basename(tools::file_path_sans_ext(cmd$msFName)), cmd$cmpName))

    summary <- file.path(resultPath, "summary_csi_fingerid.csv")
    if (file.exists(summary)) # csi:fingerid got any results?
    {
        results <- fread(summary)

        if (!is.null(cmd$topMost))
        {
            if (nrow(results) > cmd$topMost)
                results <- results[seq_len(cmd$topMost)] # results should already be sorted on score
        }

        results <- unifySirNames(results)

        # NOTE: fragment info is based on SIRIUS results, ie from formula prediction and not by compounds!
        pat <- "[:0-9:]+_([A-Za-z0-9]+).*\\.ms"
        fragFiles <- list.files(file.path(resultPath, "spectra"), full.names = TRUE, pattern = pat)
        for (ff in fragFiles)
        {
            precursor <- gsub(pat, "\\1", basename(ff))

            if (precursor %in% results$formula) # may not be there if filtered out or no compound was found
            {
                fragInfo <- fread(ff)
                fragInfo[, c("rel.intensity", "exactmass") := NULL]
                setnames(fragInfo, "explanation", "formula")
                fragInfo[, PLIndex := sapply(mz, function(omz) which.min(abs(omz - cmd$MSMSSpec$mz)))]

                set(results, which(results$formula == precursor), "fragInfo", list(list(fragInfo)))
                set(results, which(results$formula == precursor), "explainedPeaks", nrow(fragInfo))
            }
        }
        
        if (is.null(results[["fragInfo"]]))
        {
            warning("no fragment info for %s")
            results[, fragInfo := list(rep(list(data.table()), nrow(results)))]
            results[, explainedPeaks := 0]
        }

        results[, analysis := cmd$analysis]
        results[, database := cmd$database]
    }
    else
        results <- data.table()

    saveCacheData("identifySirius", results, cmd$hash, cmd$cacheDB)
    return(results)
}

#' @details \code{generateCompoundsSirius} uses
#'   \href{https://bio.informatik.uni-jena.de/software/sirius/}{SIRIUS} in
#'   combination with \href{https://www.csi-fingerid.uni-jena.de/}{CSI:FingerID}
#'   for compound identification. Similar to
#'   \code{\link{generateFormulasSirius}}, candidate formulae are generated with
#'   SIRIUS. These results are then feed to CSI:FingerID to acquire candidate
#'   structures. Candidates are then scored on similarity of measured and
#'   in-silico predicted MS/MS fragments and if the PubChem database is used, on
#'   the number of references in patents. This method requires the availability
#'   of MS/MS data, and feature groups without it will be ignored.
#'
#' @templateVar ident TRUE
#' @template sirius-args
#'
#' @templateVar genForm FALSE
#' @template form-args
#'
#' @param fingerIDDatabase Database specifically used for
#'   \command{CSI:FingerID}. If \code{NULL}, the value of the
#'   \code{formulaDatabase} parameter will be used or \code{"pubchem"} when that
#'   is also \code{NULL}. Sets the \option{--fingerid-db} option.
#'
#' @references \insertRef{Duhrkop2015}{patRoon} \cr\cr
#'   \insertRef{Duhrkop2015-2}{patRoon} \cr\cr
#'   \insertRef{Bcker2008}{patRoon}
#'
#' @rdname compound-generation
#' @export
generateCompoundsSirius <- function(fGroups, MSPeakLists, maxMzDev = 5, adduct = "[M+H]+", elements = "CHNOP",
                                    profile = "qtof", formulaDatabase = NULL, fingerIDDatabase = "pubchem",
                                    noise = NULL, topMost = 100, logPath = file.path("log", "sirius"),
                                    maxProcAmount = getOption("patRoon.maxProcAmount"))
{
    anaInfo <- analysisInfo(fGroups)
    fTable <- featureTable(fGroups)
    ftind <- groupFeatIndex(fGroups)
    gTable <- groups(fGroups)
    gNames <- colnames(gTable)
    gCount <- length(fGroups)
    gInfo <- groupInfo(fGroups)
    pLists <- peakLists(MSPeakLists)

    cacheDB <- openCacheDB()
    setHash <- makeHash(fGroups, pLists, profile, adduct, maxMzDev, elements, formulaDatabase,
                        fingerIDDatabase, noise, topMost)
    cachedSet <- loadCacheSet("identifySirius", setHash, cacheDB)
    resultHashes <- vector("character", gCount)
    names(resultHashes) <- gNames

    printf("Processing %d feature groups with SIRIUS-CSI:FingerID...\n", gCount)

    cmdQueue <- sapply(gNames, function(grp)
    {
        ogind <- order(-gTable[[grp]])
        oanalyses <- anaInfo$analysis[ogind] # analyses ordered from highest to lowest intensity

        # filter out analyses without MS/MS
        oanalyses <- sapply(oanalyses, function(a) if (!is.null(pLists[[a]][[grp]]$MSMS)) a else "", USE.NAMES = FALSE)
        ogind <- ogind[oanalyses != ""]
        oanalyses <- oanalyses[oanalyses != ""]

        if (length(oanalyses) < 1)
            return(NULL)

        ana <- oanalyses[[1]] # take most sensitive analysis

        ftmz <- fTable[[ana]][["mz"]][ftind[[grp]][ogind[1]]]
        plmz <- getMZFromMSPeakList(ftmz, pLists[[ana]][[grp]]$MS)

        hash <- makeHash(plmz, pLists[[ana]][[grp]], profile, adduct, maxMzDev, elements,
                         formulaDatabase, fingerIDDatabase, noise, topMost)
        resultHashes[[grp]] <<- hash

        cmd <- getSiriusCommand(plmz, pLists[[ana]][[grp]]$MS, pLists[[ana]][[grp]]$MSMS, profile,
                                adduct, maxMzDev, elements, formulaDatabase, noise, TRUE,
                                fingerIDDatabase)
        db <- if (!is.null(fingerIDDatabase)) fingerIDDatabase else if (!is.null(formulaDatabase)) formulaDatabase else "pubchem"
        logf <- if (!is.null(logPath)) file.path(logPath, paste0("sirius-comp-", grp, ".txt")) else NULL

        return(c(list(hash = hash, adduct = adduct, cacheDB = cacheDB, MSMSSpec = pLists[[ana]][[grp]]$MSMS,
                      analysis = ana, database = db, topMost = topMost, stdoutFile = logf), cmd))
    }, simplify = FALSE)
    cmdQueue <- cmdQueue[!sapply(cmdQueue, is.null)]

    cachedResults <- sapply(cmdQueue, function(cmd)
    {
        cr <- NULL
        if (!is.null(cachedSet))
            cr <- cachedSet[[cmd$hash]]
        if (is.null(cr))
            cr <- loadCacheData("identifySirius", cmd$hash, cacheDB)
        return(cr)
    }, simplify = FALSE)
    cachedResults <- cachedResults[!sapply(cachedResults, is.null)]

    cmdQueue <- cmdQueue[setdiff(names(cmdQueue), names(cachedResults))] # remove cached results

    if (length(cmdQueue) > 0)
    {
        if (!is.null(logPath))
            mkdirp(logPath)
        ret <- executeMultiProcess(cmdQueue, processSiriusCompounds, maxProcAmount = maxProcAmount)
    }
    else
        ret <- list()

    if (length(cachedResults) > 0)
    {
        ret <- c(ret, cachedResults)
        ret <- ret[intersect(gNames, names(ret))] # re-order
    }

    # prune empty/NULL results
    ret <- ret[sapply(ret, function(r) !is.null(r) && nrow(r) > 0, USE.NAMES = FALSE)]

    ngrp <- length(ret)
    printf("Loaded %d compounds from %d features (%.2f%%).\n", sum(unlist(lapply(ret, nrow))),
           ngrp, ngrp * 100 / gCount)

    if (is.null(cachedSet))
        saveCacheSet("identifySirius", resultHashes[resultHashes != ""], setHash, cacheDB)

    closeCacheDB(cacheDB)

    return(compounds(compounds = ret, algorithm = "SIRIUS"))
}
