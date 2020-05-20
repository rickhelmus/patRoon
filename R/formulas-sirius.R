#' @include main.R
#' @include mspeaklists.R
#' @include formulas.R
#' @include utils-sirius.R
NULL

# callback for executeMultiProcess()
processSiriusFormulas <- function(cmd, exitStatus, retries)
{
    noResult <- forms <- data.table(neutral_formula = character(0), formula = character(0),
                                    adduct = character(0), score = numeric(0), MSMSScore = numeric(0),
                                    isoScore = numeric(0), byMSMS = logical(0),
                                    frag_neutral_formula = character(0), frag_formula = character(0),
                                    frag_mz = numeric(0), frag_formula_mz = numeric(0), frag_intensity = numeric(0),
                                    neutral_loss = character(0), explainedPeaks = integer(0), explainedIntensity = numeric(0))

    resultPath <- getSiriusResultPath(cmd$outPath, cmd$msFName, cmd$cmpName, cmd$isPre44)

    summary <- file.path(resultPath, if (cmd$isPre44) "summary_sirius.csv" else "formula_candidates.csv")
    if (!file.exists(summary))
        forms <- noResult
    else
    {
        forms <- fread(summary)
        fragFiles <- getSiriusFragFiles(resultPath, cmd$isPre44)

        if (nrow(forms) == 0 || length(fragFiles) == 0)
            forms <- noResult
        else
        {
            if (cmd$isPre44)
                setnames(forms, c("formula", "treeScore"), c("neutral_formula", "MSMSScore"))
            else
            {
                setnames(forms, c("molecularFormula", "TreeIsotope_Score", "Tree_Score", "Isotope_Score"),
                         c("neutral_formula", "score", "MSMSScore", "isoScore"))
                forms[, precursorFormula := NULL] # seems same as molecularFormula
            }
            
            setkey(forms, neutral_formula)

            frags <- rbindlist(lapply(fragFiles, function(ff)
            {
                fragInfo <- fread(ff)
                setnames(fragInfo,
                         c("mz", "intensity", "exactmass", if (cmd$isPre44) "explanation" else "formula"),
                         c("frag_mz", "frag_intensity", "frag_formula_mz", "frag_neutral_formula"))
                fragInfo[, rel.intensity := NULL]
                if (!cmd$isPre44)
                    fragInfo[, ionization := NULL]
                fragInfo[, neutral_formula := getFormulaFromSiriusFragFile(ff, cmd$isPre44)]
                return(fragInfo)
            }))
            setkey(frags, neutral_formula)

            forms <- forms[frags] # merge fragment info

            forms[, formula := calculateIonFormula(neutral_formula, cmd$adduct)]
            forms[, frag_formula := calculateIonFormula(frag_neutral_formula, cmd$adduct)]
            forms[, neutral_loss := as.character(Vectorize(subtractFormula)(formula, frag_formula))]
            forms[, byMSMS := TRUE]
            forms[, rank := NULL]

            # Precursor is always present in MS/MS spectrum: it's added by
            # SIRIUS if necessarily (with zero intensity). Remove it and use its
            # frag_formula_mz to get the formula_mz
            forms[, formula_mz := .SD[frag_formula == formula, frag_formula_mz], by = "formula"]
            forms <- forms[frag_intensity != 0 | formula != frag_formula]

            # set nice column order
            setcolorder(forms, c("neutral_formula", "formula", "adduct", "score", "MSMSScore", "isoScore", "byMSMS",
                                 "frag_neutral_formula", "frag_formula", "frag_mz", "frag_formula_mz", "frag_intensity", "neutral_loss",
                                 "explainedPeaks", "explainedIntensity"))

            forms <- rankFormulaTable(forms)
        }
    }

    saveCacheData("formulasSirius", forms, cmd$hash, cmd$cacheDB)

    return(forms)
}

processSIRIUSFormulas <- function(gNames, adduct, hashes, runData, isPre44, cacheDB)
{
    noResult <- forms <- data.table(neutral_formula = character(0), formula = character(0),
                                    adduct = character(0), score = numeric(0), MSMSScore = numeric(0),
                                    isoScore = numeric(0), byMSMS = logical(0),
                                    frag_neutral_formula = character(0), frag_formula = character(0),
                                    frag_mz = numeric(0), frag_formula_mz = numeric(0), frag_intensity = numeric(0),
                                    neutral_loss = character(0), explainedPeaks = integer(0), explainedIntensity = numeric(0))
    
    
    resFile <- if (isPre44) "summary_sirius.csv" else "formula_candidates.csv"
    
    ret <- mapply(runData$msFNames, hashes, SIMPLIFY = FALSE, FUN = function(msFName, hash)
    {
        resultPath <- getSiriusResultPath(runData$outPath, msFName, runData$cmpName, isPre44)
        summary <- file.path(resultPath, resFile)
        if (length(summary) == 0 || length(summary) == 0 || !file.exists(summary))
            forms <- noResult
        else
        {
            forms <- fread(summary)
            fragFiles <- getSiriusFragFiles(resultPath, isPre44)
            
            if (nrow(forms) == 0 || length(fragFiles) == 0)
                forms <- noResult
            else
            {
                if (isPre44)
                    setnames(forms, c("formula", "treeScore"), c("neutral_formula", "MSMSScore"))
                else
                {
                    setnames(forms, c("molecularFormula", "TreeIsotope_Score", "Tree_Score", "Isotope_Score"),
                             c("neutral_formula", "score", "MSMSScore", "isoScore"))
                    forms[, precursorFormula := NULL] # seems same as molecularFormula
                }
                
                setkey(forms, neutral_formula)
                
                frags <- rbindlist(lapply(fragFiles, function(ff)
                {
                    fragInfo <- fread(ff)
                    setnames(fragInfo,
                             c("mz", "intensity", "exactmass", if (isPre44) "explanation" else "formula"),
                             c("frag_mz", "frag_intensity", "frag_formula_mz", "frag_neutral_formula"))
                    fragInfo[, rel.intensity := NULL]
                    if (!isPre44)
                        fragInfo[, ionization := NULL]
                    fragInfo[, neutral_formula := getFormulaFromSiriusFragFile(ff, isPre44)]
                    return(fragInfo)
                }))
                setkey(frags, neutral_formula)
                
                forms <- forms[frags] # merge fragment info
                
                forms[, formula := calculateIonFormula(neutral_formula, ..adduct)]
                forms[, frag_formula := calculateIonFormula(frag_neutral_formula, ..adduct)]
                forms[, neutral_loss := as.character(Vectorize(subtractFormula)(formula, frag_formula))]
                forms[, byMSMS := TRUE]
                forms[, rank := NULL]
                
                # Precursor is always present in MS/MS spectrum: it's added by
                # SIRIUS if necessarily (with zero intensity). Remove it and use its
                # frag_formula_mz to get the formula_mz
                forms[, formula_mz := .SD[frag_formula == formula, frag_formula_mz], by = "formula"]
                forms <- forms[frag_intensity != 0 | formula != frag_formula]
                
                # set nice column order
                setcolorder(forms, c("neutral_formula", "formula", "adduct", "score", "MSMSScore", "isoScore", "byMSMS",
                                     "frag_neutral_formula", "frag_formula", "frag_mz", "frag_formula_mz", "frag_intensity", "neutral_loss",
                                     "explainedPeaks", "explainedIntensity"))
                
                forms <- rankFormulaTable(forms)
            }
        }
        
        saveCacheData("formulasSIRIUS", forms, hash, cacheDB)
        return(forms)
    })
    names(ret) <- gNames
    
    return(ret)
}

processSIRIUSFormulas2 <- function(msFName, outPath, cmpName, adduct, hash, isPre44, cacheDB, ...)
{
    noResult <- forms <- data.table(neutral_formula = character(0), formula = character(0),
                                    adduct = character(0), score = numeric(0), MSMSScore = numeric(0),
                                    isoScore = numeric(0), byMSMS = logical(0),
                                    frag_neutral_formula = character(0), frag_formula = character(0),
                                    frag_mz = numeric(0), frag_formula_mz = numeric(0), frag_intensity = numeric(0),
                                    neutral_loss = character(0), explainedPeaks = integer(0), explainedIntensity = numeric(0))
    
    
    resFile <- if (isPre44) "summary_sirius.csv" else "formula_candidates.csv"
    
    resultPath <- getSiriusResultPath(outPath, msFName, cmpName, isPre44)
    summary <- file.path(resultPath, resFile)
    if (length(summary) == 0 || length(summary) == 0 || !file.exists(summary))
        forms <- noResult
    else
    {
        forms <- fread(summary)
        fragFiles <- getSiriusFragFiles(resultPath, isPre44)
        
        if (nrow(forms) == 0 || length(fragFiles) == 0)
            forms <- noResult
        else
        {
            if (isPre44)
                setnames(forms, c("formula", "treeScore"), c("neutral_formula", "MSMSScore"))
            else
            {
                setnames(forms, c("molecularFormula", "TreeIsotope_Score", "Tree_Score", "Isotope_Score"),
                         c("neutral_formula", "score", "MSMSScore", "isoScore"))
                forms[, precursorFormula := NULL] # seems same as molecularFormula
            }
            
            setkey(forms, neutral_formula)
            
            frags <- rbindlist(lapply(fragFiles, function(ff)
            {
                fragInfo <- fread(ff)
                setnames(fragInfo,
                         c("mz", "intensity", "exactmass", if (isPre44) "explanation" else "formula"),
                         c("frag_mz", "frag_intensity", "frag_formula_mz", "frag_neutral_formula"))
                fragInfo[, rel.intensity := NULL]
                if (!isPre44)
                    fragInfo[, ionization := NULL]
                fragInfo[, neutral_formula := getFormulaFromSiriusFragFile(ff, isPre44)]
                return(fragInfo)
            }))
            setkey(frags, neutral_formula)
            
            forms <- forms[frags] # merge fragment info
            
            forms[, formula := calculateIonFormula(neutral_formula, ..adduct)]
            forms[, frag_formula := calculateIonFormula(frag_neutral_formula, ..adduct)]
            forms[, neutral_loss := as.character(Vectorize(subtractFormula)(formula, frag_formula))]
            forms[, byMSMS := TRUE]
            forms[, rank := NULL]
            
            # Precursor is always present in MS/MS spectrum: it's added by
            # SIRIUS if necessarily (with zero intensity). Remove it and use its
            # frag_formula_mz to get the formula_mz
            forms[, formula_mz := .SD[frag_formula == formula, frag_formula_mz], by = "formula"]
            forms <- forms[frag_intensity != 0 | formula != frag_formula]
            
            # set nice column order
            setcolorder(forms, c("neutral_formula", "formula", "adduct", "score", "MSMSScore", "isoScore", "byMSMS",
                                 "frag_neutral_formula", "frag_formula", "frag_mz", "frag_formula_mz", "frag_intensity", "neutral_loss",
                                 "explainedPeaks", "explainedIntensity"))
            
            forms <- rankFormulaTable(forms)
        }
    }
    saveCacheData("formulasSIRIUS", forms, hash, cacheDB)
    return(forms)
}

#' @details \code{generateFormulasSirius} uses
#'   \href{https://bio.informatik.uni-jena.de/software/sirius/}{SIRIUS} to
#'   generate chemical formulae. Similarity of measured and theoretical isotopic
#'   patterns will be used for scoring candidates. Note that \command{SIRIUS}
#'   requires availability of MS/MS data.
#'
#' @templateVar ident FALSE
#' @template sirius-args
#'
#' @references \insertRef{Duhrkop2015}{patRoon} \cr\cr
#'   \insertRef{Bcker2008}{patRoon}
#'
#' @rdname formula-generation
#' @export
generateFormulasSirius <- function(fGroups, MSPeakLists, relMzDev = 5, adduct = "[M+H]+", elements = "CHNOP",
                                   profile = "qtof", database = NULL, noise = NULL, topMost = 100,
                                   extraOpts = NULL, calculateFeatures = TRUE, featThreshold = 0.75,
                                   logPath = file.path("log", "sirius"), maxProcAmount = getOption("patRoon.maxProcAmount"))
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertNumber(relMzDev, lower = 0, finite = TRUE, add = ac)
    aapply(checkmate::assertString, . ~ elements + profile, fixed = list(add = ac))
    checkmate::assertString(database, null.ok = TRUE, add = ac)
    checkmate::assertNumber(noise, lower = 0, finite = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertCount(topMost, positive = TRUE, add = ac)
    checkmate::assertCharacter(extraOpts, null.ok = TRUE, add = ac)
    checkmate::assertFlag(calculateFeatures, add = ac)
    checkmate::assertNumber(featThreshold, lower = 0, finite = TRUE, null.ok = TRUE, add = ac)
    assertMultiProcArgs(logPath, maxProcAmount, add = ac)
    checkmate::reportAssertions(ac)

    adduct <- checkAndToAdduct(adduct)

    anaInfo <- analysisInfo(fGroups)
    fTable <- featureTable(fGroups)
    featIndex <- groupFeatIndex(fGroups)
    gTable <- groups(fGroups)
    gNames <- colnames(gTable)
    gInfo <- groupInfo(fGroups)
    gCount <- length(fGroups)

    cacheDB <- openCacheDBScope() # open manually so caching code doesn't need to on each R/W access
    baseHash <- makeHash(profile, adduct, relMzDev, elements, database, noise, topMost, extraOpts)
    setHash <- makeHash(fGroups, MSPeakLists, baseHash)
    cachedSet <- loadCacheSet("formulasSirius", setHash, cacheDB)
    formHashes <- character(0)
    
    isPre44 <- isSIRIUSPre44()

    doSIRIUS <- function(featMZs, groupPeakLists)
    {
        doFGroups <- names(featMZs)

        # skip fgroups without proper ms peak lists
        if (length(doFGroups) > 0)
            doFGroups <- doFGroups[sapply(doFGroups,
                                          function(grp) !is.null(groupPeakLists[[grp]][["MS"]]) &&
                                                        any(groupPeakLists[[grp]][["MS"]]$precursor) &&
                                                        !is.null(groupPeakLists[[grp]][["MSMS"]]))]

        hashes <- sapply(doFGroups, function(grp) makeHash(featMZs[grp], groupPeakLists[[grp]], baseHash))
        formHashes <<- c(formHashes, hashes)

        cachedResults <- pruneList(sapply(hashes, function(h)
        {
            forms <- NULL
            if (!is.null(cachedSet))
                forms <- cachedSet[[h]]
            if (is.null(forms))
                forms <- loadCacheData("formulasSirius", h, cacheDB)
            return(forms)
        }, simplify = FALSE))
        doFGroups <- setdiff(doFGroups, names(cachedResults))

        cmdQueue <- pruneList(sapply(doFGroups, function(grp)
        {
            plmz <- groupPeakLists[[grp]][["MS"]][precursor == TRUE, mz]

            cmd <- getSiriusCommand(plmz, groupPeakLists[[grp]][["MS"]], groupPeakLists[[grp]][["MSMS"]], profile,
                                    adduct, relMzDev, elements, database, noise, FALSE, NULL, topMost, extraOpts,
                                    isPre44)
            logf <- if (!is.null(logPath)) file.path(logPath, paste0("sirius-form-", grp, ".txt")) else NULL

            return(c(list(hash = hashes[grp], adduct = adduct, cacheDB = cacheDB, logFile = logf,
                          gName = grp), cmd))
        }, simplify = FALSE))

        if (length(cmdQueue) > 0)
        {
            if (!is.null(logPath))
                mkdirp(logPath)

            ret <- executeMultiProcess(cmdQueue, processSiriusFormulas, errorHandler = function(cmd, exitStatus, retries)
            {
                stop(sprintf("Fatal: Failed to execute SIRIUS for %s - exit code: %d - Log: %s", cmd$gName,
                             exitStatus, cmd$logFile))
            }, maxProcAmount = maxProcAmount)

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
            ret <- ret[intersect(gNames, names(ret))] # re-order
        }

        # prune after combining with cached results: these may also contain zero row results
        ret <- pruneList(ret, checkZeroRows = TRUE)

        printf("Loaded %d formulas for %d %s (%.2f%%).\n", countUniqueFormulas(ret),
               ngrp, if (calculateFeatures) "features" else "feature groups",
               if (gCount == 0) 0 else ngrp * 100 / gCount)

        return(ret)
    }

    if (calculateFeatures)
    {
        pLists <- peakLists(MSPeakLists)

        formTable <- sapply(seq_along(anaInfo$analysis), function(anai)
        {
            ana <- anaInfo$analysis[anai]

            if (is.null(pLists[[ana]]))
                return(NULL)

            ftinds <- sapply(gNames, function(grp) featIndex[[grp]][anai])
            ftinds <- ftinds[ftinds != 0] # prune missing
            featMZs <- sapply(ftinds, function(fti) fTable[[ana]][["mz"]][fti])

            printf("Loading all formulas for analysis '%s'...\n", ana)
            return(doSIRIUS(featMZs, pLists[[ana]]))
        }, simplify = FALSE)

        names(formTable) <- anaInfo$analysis
        formTable <- pruneList(formTable, TRUE)

        if (is.null(cachedSet))
            saveCacheSet("formulasGenForm", formHashes, setHash, cacheDB)

        if (length(formTable) > 0)
            groupFormulas <- generateGroupFormulasByConsensus(formTable, featThreshold, gNames)
        else
            groupFormulas <- list()
    }
    else
    {
        featMZs <- setNames(gInfo[, "mzs"], gNames)
        groupFormulas <- doSIRIUS(featMZs, averagedPeakLists(MSPeakLists))
        formTable <- list()
    }

    return(formulas(formulas = groupFormulas, featureFormulas = formTable, algorithm = "sirius"))

}

generateFormulasSIRIUS <- function(fGroups, MSPeakLists, relMzDev = 5, adduct = "[M+H]+", elements = "CHNOP",
                                   profile = "qtof", database = NULL, noise = NULL, cores = NULL, topMost = 100,
                                   extraOptsGeneral = NULL, extraOptsFormula = NULL, calculateFeatures = TRUE,
                                   featThreshold = 0.75, verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertNumber(relMzDev, lower = 0, finite = TRUE, add = ac)
    aapply(checkmate::assertString, . ~ elements + profile, fixed = list(add = ac))
    checkmate::assertString(database, null.ok = TRUE, add = ac)
    checkmate::assertNumber(noise, lower = 0, finite = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertCount(cores, positive = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertCount(topMost, positive = TRUE, add = ac)
    aapply(checkmate::assertCharacter, . ~ extraOptsGeneral + extraOptsFormula, null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertFlag(calculateFeatures, add = ac)
    checkmate::assertNumber(featThreshold, lower = 0, finite = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    adduct <- checkAndToAdduct(adduct)
    
    anaInfo <- analysisInfo(fGroups)
    fTable <- featureTable(fGroups)
    featIndex <- groupFeatIndex(fGroups)
    gTable <- groups(fGroups)
    gNames <- colnames(gTable)
    gInfo <- groupInfo(fGroups)
    gCount <- length(fGroups)
    
    cacheDB <- openCacheDBScope() # open manually so caching code doesn't need to on each R/W access
    baseHash <- makeHash(profile, adduct, relMzDev, elements, database, noise, topMost, extraOptsGeneral, extraOptsFormula)
    setHash <- makeHash(fGroups, MSPeakLists, baseHash)
    cachedSet <- loadCacheSet("formulasSIRIUS", setHash, cacheDB)
    formHashes <- character(0)
    
    isPre44 <- isSIRIUSPre44()
    
    if (calculateFeatures)
    {
        pLists <- peakLists(MSPeakLists)
        
        formTable <- sapply(seq_along(anaInfo$analysis), function(anai)
        {
            ana <- anaInfo$analysis[anai]
            
            if (is.null(pLists[[ana]]))
                return(NULL)
            
            ftinds <- sapply(gNames, function(grp) featIndex[[grp]][anai])
            ftinds <- ftinds[ftinds != 0] # prune missing
            featMZs <- sapply(ftinds, function(fti) fTable[[ana]][["mz"]][fti])
            
            printf("Loading all formulas for analysis '%s'...\n", ana)
            return(doSIRIUS(gNames, featMZs, pLists[[ana]], profile, adduct, relMzDev, elements,
                            database, noise, cores, FALSE, NULL, topMost, extraOptsGeneral, extraOptsFormula,
                            verbose, isPre44, "formulasSIRIUS", cacheDB, processSIRIUSFormulas))
        }, simplify = FALSE)
        
        names(formTable) <- anaInfo$analysis
        formTable <- pruneList(formTable, TRUE)
        
        if (length(formTable) > 0)
            groupFormulas <- generateGroupFormulasByConsensus(formTable, featThreshold, gNames)
        else
            groupFormulas <- list()
    }
    else
    {
        featMZs <- setNames(gInfo[, "mzs"], gNames)
        # groupFormulas <- doSIRIUS(featMZs, averagedPeakLists(MSPeakLists))
        groupFormulas <- doSIRIUS(gNames, featMZs, averagedPeakLists(MSPeakLists), profile, adduct, relMzDev, elements,
                                  database, noise, cores, FALSE, NULL, topMost, extraOptsGeneral, extraOptsFormula,
                                  verbose, isPre44, "formulasSIRIUS", cacheDB, processSIRIUSFormulas)
        formTable <- list()
    }

    # if (is.null(cachedSet))
    #     saveCacheSet("formulasSIRIUS", formHashes, setHash, cacheDB)
    
    return(formulas(formulas = groupFormulas, featureFormulas = formTable, algorithm = "sirius"))
}

generateFormulasSIRIUS2 <- function(fGroups, MSPeakLists, relMzDev = 5, adduct = "[M+H]+", elements = "CHNOP",
                                    profile = "qtof", database = NULL, noise = NULL, cores = NULL, topMost = 100,
                                    extraOptsGeneral = NULL, extraOptsFormula = NULL, calculateFeatures = TRUE,
                                    featThreshold = 0.75, verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertNumber(relMzDev, lower = 0, finite = TRUE, add = ac)
    aapply(checkmate::assertString, . ~ elements + profile, fixed = list(add = ac))
    checkmate::assertString(database, null.ok = TRUE, add = ac)
    checkmate::assertNumber(noise, lower = 0, finite = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertCount(cores, positive = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertCount(topMost, positive = TRUE, add = ac)
    aapply(checkmate::assertCharacter, . ~ extraOptsGeneral + extraOptsFormula, null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertFlag(calculateFeatures, add = ac)
    checkmate::assertNumber(featThreshold, lower = 0, finite = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    adduct <- checkAndToAdduct(adduct)
    
    gNames <- names(fGroups)
    
    formTable <- doSIRIUS2(gNames, MSPeakLists, calculateFeatures, profile, adduct, relMzDev, elements,
                           database, noise, cores, FALSE, NULL, topMost, extraOptsGeneral, extraOptsFormula,
                           verbose, "formulasSIRIUS", processSIRIUSFormulas2, NULL)
        
    if (calculateFeatures)
    {
        if (length(formTable) > 0)
        {
            formTable <- lapply(formTable, pruneList, checkZeroRows = TRUE)
            groupFormulas <- generateGroupFormulasByConsensus(formTable, featThreshold, gNames)
        }
        else
            groupFormulas <- list()
    }
    else
    {
        groupFormulas <- pruneList(formTable, checkZeroRows = TRUE)
        formTable <- list()
    }
    
    if (verbose)
    {
        printf("-------\n")
        if (calculateFeatures)
        {
            fCounts <- sapply(formTable, countUniqueFormulas)
            fTotCount <- sum(fCounts)
            printf("Formula statistics:\n")
            printf("%s: %d (%.1f%%)\n", names(formTable), fCounts, if (fTotCount == 0) 0 else fCounts * 100 / fTotCount)
            printf("Total: %d\n", fTotCount)
        }
        ngrp <- length(groupFormulas)
        gCount <- length(fGroups)
        printf("Assigned %d unique formulas to %d feature groups (%.2f%% coverage).\n", countUniqueFormulas(groupFormulas),
               ngrp, if (gCount == 0) 0 else ngrp * 100 / gCount)
    }
    
    return(formulas(formulas = groupFormulas, featureFormulas = formTable, algorithm = "sirius"))
}
