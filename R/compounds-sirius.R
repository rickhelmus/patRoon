#' @include main.R
#' @include compounds.R
#' @include utils-sirius.R
NULL

processSiriusCompounds <- function(cmd, exitStatus, retries)
{
    resultPath <- getSiriusResultPath(cmd$outPath, cmd$msFName, cmd$cmpName, cmd$isPre44)

    results <- data.table()
    scRanges <- list()

    summary <- file.path(resultPath, if (cmd$isPre44) "summary_csi_fingerid.csv" else "structure_candidates.csv")
    if (length(summary) != 0 && file.exists(summary)) # csi:fingerid got any results?
    {
        results <- fread(summary)
        results <- unifySirNames(results)
        
        # NOTE: so far SIRIUS only has one score
        if (nrow(results) > 0)
            scRanges <- list(score = range(results$score))

        if (!is.null(cmd$topMost))
        {
            if (nrow(results) > cmd$topMost)
                results <- results[seq_len(cmd$topMost)] # results should already be sorted on score
        }
        
        # manually sort results: SIRIUS seems to sometimes randomly order results with equal scores
        setorderv(results, c("score", "identifier"), order = c(-1, 1))

        # NOTE: fragment info is based on SIRIUS results, ie from formula
        # prediction and not by compounds! Hence, results are the same for all
        # isomers.
        fragFiles <- getSiriusFragFiles(resultPath, cmd$isPre44)
        for (ff in fragFiles)
        {
            precursor <- getFormulaFromSiriusFragFile(ff, cmd$isPre44)
            if (precursor %in% results$formula) # may not be there if filtered out or no compound was found
            {
                fragInfo <- fread(ff)
                fragInfo[, c("rel.intensity", "exactmass") := NULL]
                if (cmd$isPre44)
                    setnames(fragInfo, "explanation", "formula")
                fragInfo[, PLIndex := sapply(mz, function(omz) which.min(abs(omz - cmd$MSMSSpec$mz)))]

                # sirius neutralizes fragments, make them ion again
                fragInfo[, formula := calculateIonFormula(formula, cmd$adduct)]
                fragInfo[, neutral_loss := sapply(formula, subtractFormula,
                                                  formula1 = calculateIonFormula(precursor, cmd$adduct))]

                set(results, which(results$formula == precursor), "fragInfo", list(list(fragInfo)))
                set(results, which(results$formula == precursor), "explainedPeaks", nrow(fragInfo))
            }
        }

        if (is.null(results[["fragInfo"]]))
        {
            # warning(sprintf("no fragment info for %s", cmd$gName))
            results[, fragInfo := list(rep(list(data.table()), nrow(results)))]
            results[, explainedPeaks := 0]
        }

        results[, database := cmd$database][]
    }

    ret <- list(comptab = results, scRanges = scRanges)
    saveCacheData("compoundsSirius", ret, cmd$hash, cmd$cacheDB)
    return(ret)
}

processSIRIUSCompounds <- function(msFName, outPath, cmpName, MSMS, database, adduct, topMost, hash, isPre44, cacheDB)
{
    resFile <- if (isPre44) "summary_csi_fingerid.csv" else "structure_candidates.csv"
    resultPath <- getSiriusResultPath(outPath, msFName, cmpName, isPre44)
    summary <- file.path(resultPath, resFile)
    results <- scRanges <- list()
    
    if (length(summary) != 0 && file.exists(summary))
    {
        results <- fread(summary)
        results <- unifySirNames(results)
        
        # NOTE: so far SIRIUS only has one score
        if (nrow(results) > 0)
            scRanges <- list(score = range(results$score))
        
        # manually sort results: SIRIUS seems to sometimes randomly order results with equal scores
        setorderv(results, c("score", "identifier"), order = c(-1, 1))
        
        if (!is.null(topMost))
        {
            if (nrow(results) > topMost)
                results <- results[seq_len(topMost)] # results should already be sorted on score
        }
        
        # NOTE: fragment info is based on SIRIUS results, ie from formula
        # prediction and not by compounds! Hence, results are the same for
        # candidates with the same formula.
        fragFiles <- getSiriusFragFiles(resultPath, isPre44)
        for (ff in fragFiles)
        {
            precursor <- getFormulaFromSiriusFragFile(ff, isPre44)
            if (precursor %in% results$formula) # may not be there if filtered out or no compound was found
            {
                fragInfo <- fread(ff)
                fragInfo[, c("rel.intensity", "exactmass") := NULL]
                if (isPre44)
                    setnames(fragInfo, "explanation", "formula")
                fragInfo[, PLIndex := sapply(mz, function(omz) which.min(abs(omz - MSMS$mz)))]
                
                # sirius neutralizes fragments, make them ion again
                fragInfo[, formula := calculateIonFormula(formula, ..adduct)]
                fragInfo[, neutral_loss := sapply(formula, subtractFormula,
                                                  formula1 = calculateIonFormula(precursor, ..adduct))]
                
                set(results, which(results$formula == precursor), "fragInfo", list(list(fragInfo)))
                set(results, which(results$formula == precursor), "explainedPeaks", nrow(fragInfo))
            }
        }
        
        if (is.null(results[["fragInfo"]]))
        {
            # warning(sprintf("no fragment info for %s", cmd$gName))
            results[, fragInfo := list(rep(list(data.table()), nrow(results)))]
            results[, explainedPeaks := 0]
        }
        
        results[, database := database][]
    }
    
    ret <- list(comptab = results, scRanges = scRanges)
    saveCacheData("compoundsSIRIUS", ret, hash, cacheDB)
    return(ret)
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
#' @return \code{generateCompoundsSirius} returns a \code{\link{compounds}}
#'   object.
#' @param fingerIDDatabase Database specifically used for
#'   \command{CSI:FingerID}. If \code{NULL}, the value of the
#'   \code{formulaDatabase} parameter will be used or \code{"pubchem"} when that
#'   is also \code{NULL}. Sets the \option{--fingerid-db} option.
#' @param topMostFormulas Do not return more than this number of candidate
#'   formulae. Note that only compounds for these formulae will be searched.
#'   Sets the \option{--candidates} commandline option.
#'
#' @references \insertRef{Duhrkop2015}{patRoon} \cr\cr
#'   \insertRef{Duhrkop2015-2}{patRoon} \cr\cr \insertRef{Bcker2008}{patRoon}
#'
#' @rdname compound-generation
#' @export
generateCompoundsSirius <- function(fGroups, MSPeakLists, relMzDev = 5, adduct = "[M+H]+", elements = "CHNOP",
                                    profile = "qtof", formulaDatabase = NULL, fingerIDDatabase = "pubchem",
                                    noise = NULL, errorRetries = 2, topMost = 100, topMostFormulas = 5, extraOpts = NULL,
                                    logPath = file.path("log", "sirius"),
                                    maxProcAmount = getOption("patRoon.maxProcAmount"))
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertNumber(relMzDev, lower = 0, finite = TRUE, add = ac)
    aapply(checkmate::assertString, . ~ elements + profile + fingerIDDatabase, fixed = list(add = ac))
    aapply(checkmate::assertString, . ~ formulaDatabase + fingerIDDatabase, null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertNumber(noise, lower = 0, finite = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertCount(errorRetries, add = ac)
    checkmate::assertCount(topMost, positive = TRUE, add = ac)
    checkmate::assertCount(topMostFormulas, positive = TRUE, add = ac)
    checkmate::assertCharacter(extraOpts, null.ok = TRUE, add = ac)
    assertMultiProcArgs(logPath, maxProcAmount, add = ac)
    checkmate::reportAssertions(ac)

    adduct <- checkAndToAdduct(adduct)
    if (is.null(fingerIDDatabase))
        fingerIDDatabase <- if (!is.null(formulaDatabase)) formulaDatabase else "pubchem"

    anaInfo <- analysisInfo(fGroups)
    fTable <- featureTable(fGroups)
    ftind <- groupFeatIndex(fGroups)
    gTable <- groups(fGroups)
    gNames <- colnames(gTable)
    gCount <- length(fGroups)
    gInfo <- groupInfo(fGroups)

    cacheDB <- openCacheDBScope()
    setHash <- makeHash(fGroups, MSPeakLists, profile, adduct, relMzDev, elements, formulaDatabase,
                        fingerIDDatabase, noise, topMost, topMostFormulas, extraOpts)
    cachedSet <- loadCacheSet("compoundsSirius", setHash, cacheDB)
    resultHashes <- vector("character", gCount)
    names(resultHashes) <- gNames

    printf("Processing %d feature groups with SIRIUS-CSI:FingerID...\n", gCount)

    isPre44 <- isSIRIUSPre44()
    cmdQueue <- sapply(gNames, function(grp)
    {
        plist <- MSPeakLists[[grp]]
        if (is.null(plist[["MS"]]) || is.null(plist[["MSMS"]]))
            return(NULL)

        plmz <- plist$MS[precursor == TRUE, mz]
        if (length(plmz) == 0)
            return(NULL)

        hash <- makeHash(plmz, plist, profile, adduct, relMzDev, elements,
                         formulaDatabase, fingerIDDatabase, noise, topMost,
                         topMostFormulas, extraOpts)
        resultHashes[[grp]] <<- hash

        cmd <- getSiriusCommand(plmz, plist$MS, plist$MSMS, profile,
                                adduct, relMzDev, elements, formulaDatabase, noise, TRUE,
                                fingerIDDatabase, topMostFormulas, extraOpts, isPre44)
        logf <- if (!is.null(logPath)) file.path(logPath, paste0("sirius-comp-", grp, ".txt")) else NULL
        return(c(list(hash = hash, adduct = adduct, cacheDB = cacheDB, MSMSSpec = plist$MSMS,
                      database = fingerIDDatabase, topMost = topMost, logFile = logf, gName = grp), cmd))
    }, simplify = FALSE)
    cmdQueue <- cmdQueue[!sapply(cmdQueue, is.null)]

    cachedResults <- sapply(cmdQueue, function(cmd)
    {
        cr <- NULL
        if (!is.null(cachedSet))
            cr <- cachedSet[[cmd$hash]]
        if (is.null(cr))
            cr <- loadCacheData("compoundsSirius", cmd$hash, cacheDB)
        return(cr)
    }, simplify = FALSE)
    cachedResults <- cachedResults[!sapply(cachedResults, is.null)]

    cmdQueue <- cmdQueue[setdiff(names(cmdQueue), names(cachedResults))] # remove cached results

    if (length(cmdQueue) > 0)
    {
        if (!is.null(logPath))
            mkdirp(logPath)
        ret <- executeMultiProcess(cmdQueue, processSiriusCompounds, errorHandler = function(cmd, exitStatus, retries) {
            if (retries < errorRetries)
            {
                warning(sprintf("Restarting failed SIRIUS for %s - exit: %d (retry %d/%d)",
                                cmd$gName, exitStatus, retries+1, errorRetries))
                unlink(cmd$outPath, TRUE)
                return(TRUE)
            }
            stop(sprintf("Fatal: Failed to execute SIRIUS for %s - exit code: %d\nLog: %s", cmd$gName, exitStatus,
                         cmd$logFile))
        }, maxProcAmount = maxProcAmount)
    }
    else
        ret <- list()

    if (length(cachedResults) > 0)
    {
        ret <- c(ret, cachedResults)
        ret <- ret[intersect(gNames, names(ret))] # re-order
    }

    # prune empty/NULL results
    if (length(ret) > 0)
        ret <- ret[sapply(ret, function(r) !is.null(r$comptab) && nrow(r$comptab) > 0, USE.NAMES = FALSE)]

    ngrp <- length(ret)
    printf("Loaded %d compounds from %d features (%.2f%%).\n", sum(unlist(lapply(ret, function(r) nrow(r$comptab)))),
           ngrp, if (gCount == 0) 0 else ngrp * 100 / gCount)

    if (is.null(cachedSet))
        saveCacheSet("compoundsSirius", resultHashes[resultHashes != ""], setHash, cacheDB)

    return(compounds(compounds = lapply(ret, "[[", "comptab"), scoreTypes = "score",
                     scoreRanges = lapply(ret, "[[", "scRanges"),
                     algorithm = "sirius"))
}

generateCompoundsSIRIUS <- function(fGroups, MSPeakLists, relMzDev = 5, adduct = "[M+H]+", elements = "CHNOP",
                                    profile = "qtof", formulaDatabase = NULL, fingerIDDatabase = "pubchem",
                                    noise = NULL, errorRetries = 2, cores = NULL, topMost = 100, topMostFormulas = 5,
                                    extraOptsGeneral = NULL, extraOptsFormula = NULL, verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertNumber(relMzDev, lower = 0, finite = TRUE, add = ac)
    aapply(checkmate::assertString, . ~ elements + profile + fingerIDDatabase, fixed = list(add = ac))
    aapply(checkmate::assertString, . ~ formulaDatabase + fingerIDDatabase, null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertNumber(noise, lower = 0, finite = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertCount(errorRetries, add = ac)
    checkmate::assertCount(cores, positive = TRUE, null.ok = TRUE, add = ac)
    aapply(checkmate::assertCount, . ~ topMost + topMostFormulas, positive = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertCharacter, . ~ extraOptsGeneral + extraOptsFormula, null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)

    gNames <- names(fGroups)
    adduct <- checkAndToAdduct(adduct)
    if (is.null(fingerIDDatabase))
        fingerIDDatabase <- if (!is.null(formulaDatabase)) formulaDatabase else "pubchem"

    gCount <- length(fGroups)
    printf("Processing %d feature groups with SIRIUS-CSI:FingerID...\n", gCount)
    
    results <- doSIRIUS2(gNames, MSPeakLists, FALSE, profile, adduct, relMzDev, elements,
                         formulaDatabase, noise, cores, TRUE, fingerIDDatabase, topMost, extraOptsGeneral, extraOptsFormula,
                         verbose, "compoundsSIRIUS", processSIRIUSCompounds,
                         list(database = fingerIDDatabase, topMost = topMost))
    
    # prune empty/NULL results
    if (length(results) > 0)
        results <- results[sapply(results, function(r) !is.null(r$comptab) && nrow(r$comptab) > 0, USE.NAMES = FALSE)]
    
    if (verbose)
    {
        printf("-------\n")
        ngrp <- length(results)
        printf("Assigned %d compounds to %d feature groups (%.2f%%).\n", sum(unlist(lapply(results, function(r) nrow(r$comptab)))),
               ngrp, if (gCount == 0) 0 else ngrp * 100 / gCount)
    }

    return(compounds(compounds = lapply(results, "[[", "comptab"), scoreTypes = "score",
                     scoreRanges = lapply(results, "[[", "scRanges"),
                     algorithm = "sirius"))
}
