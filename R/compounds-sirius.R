#' @include main.R
#' @include compounds.R
#' @include utils-sirius.R
NULL

processSIRIUSCompounds <- function(msFName, outPath, cmpName, MSMS, database, adduct, topMost, hash, isPre44, cacheDB)
{
    resFile <- if (isPre44) "summary_csi_fingerid.csv" else "structure_candidates.tsv"
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
        # NOTE: do this before topMost filter step below to ensure proper filtering
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

#' @details \code{generateCompoundsSIRIUS} uses
#'   \href{https://bio.informatik.uni-jena.de/software/sirius/}{SIRIUS} in
#'   combination with \href{https://www.csi-fingerid.uni-jena.de/}{CSI:FingerID}
#'   for compound identification. Similar to
#'   \code{\link{generateFormulasSIRIUS}}, candidate formulae are generated with
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
#' @return \code{generateCompoundsSIRIUS} returns a \code{\link{compounds}}
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
generateCompoundsSIRIUS <- function(fGroups, MSPeakLists, relMzDev = 5, adduct = "[M+H]+", elements = "CHNOP",
                                    profile = "qtof", formulaDatabase = NULL, fingerIDDatabase = "pubchem",
                                    noise = NULL, errorRetries = 2, cores = NULL, topMost = 100, topMostFormulas = 5,
                                    extraOptsGeneral = NULL, extraOptsFormula = NULL, verbose = TRUE,
                                    batchSize = 0, logPath = file.path("log", "sirius_compounds"),
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
    checkmate::assertCount(cores, positive = TRUE, null.ok = TRUE, add = ac)
    aapply(checkmate::assertCount, . ~ topMost + topMostFormulas, positive = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertCharacter, . ~ extraOptsGeneral + extraOptsFormula, null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertFlag(verbose, add = ac)
    checkmate::assertCount(batchSize, add = ac)
    assertMultiProcArgs(logPath, maxProcAmount, add = ac)
    checkmate::reportAssertions(ac)

    gNames <- names(fGroups)
    adduct <- checkAndToAdduct(adduct)
    if (is.null(fingerIDDatabase))
        fingerIDDatabase <- if (!is.null(formulaDatabase)) formulaDatabase else "pubchem"

    gCount <- length(fGroups)
    printf("Processing %d feature groups with SIRIUS-CSI:FingerID...\n", gCount)
    
    results <- doSIRIUS(gNames, MSPeakLists, FALSE, profile, adduct, relMzDev, elements,
                        formulaDatabase, noise, cores, TRUE, fingerIDDatabase, topMost, extraOptsGeneral, extraOptsFormula,
                        verbose, "compoundsSIRIUS", processSIRIUSCompounds,
                        list(database = fingerIDDatabase, topMost = topMost),
                        batchSize, logPath, maxProcAmount)
    
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
