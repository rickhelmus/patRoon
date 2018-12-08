#' @include main.R
#' @include mspeaklists.R
#' @include formulas.R
#' @include utils-sirius.R
NULL

# callback for executeMultiProcess()
processSiriusFormulas <- function(cmd, exitStatus, retries)
{
    noResult <- forms <- data.table(neutral_formula = character(0), formula = character(0),
                                    adduct = character(0), rank = integer(0), score = numeric(0), treeScore = numeric(0),
                                    isoScore = numeric(0), byMSMS = logical(0), frag_formula = character(0),
                                    frag_mz = numeric(0), frag_formula_mz = numeric(0), frag_intensity = numeric(0),
                                    neutral_loss = character(0), explainedPeaks = integer(0), explainedIntensity = numeric(0))

    # format is resultno_specname_compoundname
    resultPath <- file.path(cmd$outPath, sprintf("1_%s_%s", basename(tools::file_path_sans_ext(cmd$msFName)), cmd$cmpName))

    summary <- file.path(resultPath, "summary_sirius.csv")
    if (!file.exists(summary))
        forms <- noResult
    else
    {
        forms <- fread(summary)
        fragFiles <- getSiriusFragFiles(resultPath)

        if (nrow(forms) == 0 || length(fragFiles) == 0)
            forms <- noResult
        else
        {
            setnames(forms, "formula", "neutral_formula")
            setkey(forms, neutral_formula)

            frags <- rbindlist(lapply(fragFiles, function(ff)
            {
                fragInfo <- fread(ff)
                setnames(fragInfo,
                         c("mz", "intensity", "exactmass", "explanation"),
                         c("frag_mz", "frag_intensity", "frag_formula_mz", "frag_formula"))
                fragInfo[, rel.intensity := NULL]
                fragInfo[, neutral_formula := getFormulaFromSiriusFragFile(ff)]
                return(fragInfo)
            }))
            setkey(frags, neutral_formula)

            forms <- forms[frags] # merge fragment info

            forms[, formula := calculateIonFormula(neutral_formula, cmd$adduct)]
            forms[, neutral_loss := as.character(Vectorize(subtractFormula)(formula, frag_formula))]
            forms[, byMSMS := TRUE]

            # set nice column order
            setcolorder(forms, c("neutral_formula", "formula", "adduct", "rank", "score", "treeScore", "isoScore", "byMSMS",
                                 "frag_formula", "frag_mz", "frag_formula_mz", "frag_intensity", "neutral_loss",
                                 "explainedPeaks", "explainedIntensity"))
            
            forms <- rankFormulaTable(forms)
        }
    }

    saveCacheData("formulasSirius", forms, cmd$hash, cmd$cacheDB)

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
generateFormulasSirius <- function(fGroups, MSPeakLists, maxMzDev = 5, adduct = "[M+H]+", elements = "CHNOP",
                                   profile = "qtof", database = NULL, noise = NULL,
                                   calculateFeatures = TRUE, formFeatThreshold = 0.75,
                                   logPath = file.path("log", "sirius"), maxProcAmount = getOption("patRoon.maxProcAmount"))
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertNumber(maxMzDev, lower = 0, finite = TRUE, add = ac)
    aapply(checkmate::assertString, . ~ adduct + elements + profile, fixed = list(add = ac))
    checkmate::assertString(database, null.ok = TRUE, add = ac)
    checkmate::assertNumber(noise, lower = 0, finite = TRUE, null.ok = TRUE, add = ac)
    assertMultiProcArgs(logPath, maxProcAmount, add = ac)
    checkmate::reportAssertions(ac)
    
    anaInfo <- analysisInfo(fGroups)
    fTable <- featureTable(fGroups)
    featIndex <- groupFeatIndex(fGroups)
    gTable <- groups(fGroups)
    gNames <- colnames(gTable)
    gInfo <- groupInfo(fGroups)
    gCount <- length(fGroups)

    cacheDB <- openCacheDBScope() # open manually so caching code doesn't need to on each R/W access
    baseHash <- makeHash(profile, adduct, maxMzDev, elements, database, noise)
    setHash <- makeHash(fGroups, MSPeakLists, baseHash)
    cachedSet <- loadCacheSet("formulasSirius", setHash, cacheDB)
    formHashes <- character(0)

    doSIRIUS <- function(featMZs, groupPeakLists)
    {
        doFGroups <- names(featMZs)
        
        # skip fgroups without proper ms peak lists
        if (length(doFGroups) > 0)
            doFGroups <- doFGroups[sapply(doFGroups,
                                          function(grp) !is.null(groupPeakLists[[grp]][["MS"]]) &&
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
            plmz <- getMZFromMSPeakList(featMZs[grp], groupPeakLists[[grp]][["MS"]])

            cmd <- getSiriusCommand(plmz, groupPeakLists[[grp]][["MS"]], groupPeakLists[[grp]][["MSMS"]], profile,
                                    adduct, maxMzDev, elements, database, noise, FALSE, NULL)
            logf <- if (!is.null(logPath)) file.path(logPath, paste0("sirius-form-", grp, ".txt")) else NULL
            logfe <- if (!is.null(logPath)) file.path(logPath, paste0("sirius-form-err-", grp, ".txt")) else NULL
            
            return(c(list(hash = hashes[grp], adduct = adduct, cacheDB = cacheDB, stdoutFile = logf,
                          stderrFile = logfe, gName = grp), cmd))
        }, simplify = FALSE))
        
        if (length(cmdQueue) > 0)
        {
            if (!is.null(logPath))
                mkdirp(logPath)
            
            ret <- executeMultiProcess(cmdQueue, processSiriusFormulas, errorHandler = function(cmd, exitStatus, retries)
            {
                stop(sprintf("Fatal: Failed to execute SIRIUS for %s - exit code: %d\nCommand: %s", cmd$gName, exitStatus,
                             paste(cmd$command, paste0(cmd$args, collapse = " "))))
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
        
        printf("Loaded %d formulas for %d %s (%.2f%%).\n", sum(unlist(lapply(ret, nrow))),
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
            groupFormulas <- generateGroupFormulasByConsensus(formTable, formFeatThreshold)
        else
            groupFormulas <- list()
    }
    else
    {
        featMZs <- setNames(gInfo[, "mzs"], gNames)
        groupFormulas <- doSIRIUS(featMZs, averagedPeakLists(MSPeakLists))
        formTable <- list()
    }
    
    return(formulas(formulas = formTable, groupFormulas = groupFormulas, algorithm = "SIRIUS"))
    
}
