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
    ftind <- groupFeatIndex(fGroups)
    gTable <- groups(fGroups)
    gNames <- colnames(gTable)
    gInfo <- groupInfo(fGroups)
    gcount <- ncol(ftind)
    pLists <- peakLists(MSPeakLists)

    formTable <- list()

    cacheDB <- openCacheDB() # open manually so caching code doesn't need to on each R/W access
    setHash <- makeHash(fGroups, pLists, profile, adduct, maxMzDev, elements, database, noise)
    cachedSet <- loadCacheSet("formulasSirius", setHash, cacheDB)
    formHashes <- character(0)

    for (anai in seq_along(anaInfo$analysis))
    {
        ana <- anaInfo$analysis[anai]
        gforms <- list()

        printf("Loading all formulas from analysis '%s'...\n", ana)

        cmdQueue <- sapply(gNames, function(grp)
        {
            if (ftind[[grp]][anai] == 0)
                return(NULL) # feature not present

            if (is.null(pLists[[ana]][[grp]]) || nrow(pLists[[ana]][[grp]]$MS) == 0)
                return(NULL) # MS or MSMS spectrum probably filtered away

            if (is.null(pLists[[ana]][[grp]]$MSMS))
                return(NULL) # no MSMS

            ftmz <- fTable[[ana]][["mz"]][ftind[[grp]][anai]]
            plmz <- getMZFromMSPeakList(ftmz, pLists[[ana]][[grp]]$MS)
            hash <- makeHash(ana, grp, ftmz, pLists[[ana]][[grp]], profile, adduct, maxMzDev, elements, database, noise)

            cmd <- getSiriusCommand(plmz, pLists[[ana]][[grp]]$MS, pLists[[ana]][[grp]]$MSMS, profile,
                                    adduct, maxMzDev, elements, database, noise, FALSE, NULL)
            logf <- if (!is.null(logPath)) file.path(logPath, paste0("sirius-form-", grp, ".txt")) else NULL
            logfe <- if (!is.null(logPath)) file.path(logPath, paste0("sirius-form-err-", grp, ".txt")) else NULL

            return(c(list(hash = hash, adduct = adduct, cacheDB = cacheDB, stdoutFile = logf,
                          stderrFile = logfe, gName = grp), cmd))
        }, simplify = FALSE)
        cmdQueue <- cmdQueue[!sapply(cmdQueue, is.null)]

        formHashes <- c(formHashes, sapply(cmdQueue, function(cmd) cmd$hash))

        cachedResults <- sapply(cmdQueue, function(cmd)
        {
            forms <- NULL
            if (!is.null(cachedSet))
                forms <- cachedSet[[cmd$hash]]
            if (is.null(forms))
                forms <- loadCacheData("formulasSirius", cmd$hash, cacheDB)
            return(forms)
        }, simplify = FALSE)
        cachedResults <- cachedResults[!sapply(cachedResults, is.null)]

        cmdQueue <- cmdQueue[setdiff(names(cmdQueue), names(cachedResults))] # remove cached results

        if (length(cmdQueue) > 0)
        {
            if (!is.null(logPath))
                mkdirp(logPath)

            formTable[[ana]] <- executeMultiProcess(cmdQueue, processSiriusFormulas, errorHandler = function(cmd, exitStatus, retries) {
                stop(sprintf("Fatal: Failed to execute SIRIUS for %s - exit code: %d\nCommand: %s", cmd$gName, exitStatus,
                             paste(cmd$command, paste0(cmd$args, collapse = " "))))
                }, maxProcAmount = maxProcAmount)
            ngrp <- length(formTable[[ana]])
        }
        else
        {
            formTable[[ana]] <- list()
            ngrp <- 0
        }

        if (length(cachedResults) > 0)
        {
            ngrp <- ngrp + length(cachedResults)
            formTable[[ana]] <- c(formTable[[ana]], cachedResults)
            formTable[[ana]] <- formTable[[ana]][intersect(gNames, names(formTable[[ana]]))] # re-order
        }

        printf("Loaded %d formulas from %d features (%.2f%%).\n", sum(unlist(lapply(formTable[[ana]], nrow))),
               ngrp, if (gcount == 0) 0 else ngrp * 100 / gcount)
    }

    if (is.null(cachedSet))
        saveCacheSet("formulasSirius", formHashes, setHash, cacheDB)

    closeCacheDB(cacheDB)

    formTable <- sapply(formTable, function(ft) ft[sapply(ft, nrow) > 0], simplify = FALSE)

    return(formulas(formulas = formTable, algorithm = "SIRIUS"))
}
