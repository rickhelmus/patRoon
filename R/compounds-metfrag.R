#' @include main.R
#' @include compounds.R
#' @include mspeaklists.R
#' @include formulas.R
NULL

unifyMFNames <- function(mfr)
{
    unNames <- c(NoExplPeaks = "explainedPeaks",
                 Score = "score",
                 MonoisotopicMass = "neutralMass",
                 SMILES = "SMILES",
                 InChIKey = "InChIKey",
                 InChIKey2 = "InChIKey2",
                 InChIKey1 = "InChIKey1",
                 InChI = "InChi",
                 Identifier = "identifier",
                 PubChemNumberPatents = "numberPatents",
                 fragInfo = "fragInfo", # keep it
                 PubChemNumberPubMedReferences = "pubMedReferences",
                 ChemSpiderNumberPubMedReferences = "pubMedReferences",
                 FragmenterScore = "fragScore",
                 MolecularFormula = "formula",
                 XlogP3 = "XlogP",
                 CHEMSPIDER_XLOGP = "XlogP",
                 TrivialName = "trivialName",
                 CompoundName = "trivialName",
                 CHEMSPIDER_ALOGP = "AlogP",
                 ChemSpiderNumberExternalReferences = "extReferenceCount",
                 ChemSpiderDataSourceCount = "dataSourceCount",
                 ChemSpiderReferenceCount = "referenceCount",
                 ChemSpiderRSCCount = "RSCCount",
                 OfflineMetFusionScore = "metFusionScore",
                 OfflineIndividualMoNAScore = "individualMoNAScore",
                 SmartsSubstructureInclusionScore = "smartsInclusionScore")

    unNames <- unNames[names(unNames) %in% names(mfr)] # filter out missing
    setnames(mfr, names(unNames), unNames)

    return(mfr[, unNames, with = FALSE]) # filter out any other columns
}

# MetFragCL gives bracketed fragment formulas including charge and weird (de)protonation adducts
# For comparison with other results these should be converted to a simple formula format
cleanFragFormulas <- function(forms)
{
    forms <- gsub("\\[|\\]", "", forms) # remove brackets
    forms <- gsub("[-\\+]$", "", forms) # remove trailing charge

    # add single counts to hydrogen adducts without count (e.g. -H becomes -1H)
    forms <- gsub("([-\\+])H", "\\11H", forms)

    adducts <- regmatches(forms, gregexpr("([-\\+][0-9]+)H", forms)) # get "-1H", "+2H" etc

    # Get 'regular' part of formula and update H count
    baseForms <- gsub("^([[:alnum:]]+)[-|\\+].*", "\\1", forms)

    return(sapply(seq_along(forms), function(fi)
    {
        addHCount <- sum(as.integer(gsub("H", "", adducts[[fi]])))
        flist <- splitFormulaToList(baseForms[[fi]])

        if (addHCount != 0)
        {
            if (!"H" %in% names(flist))
            {
                flist <- c(flist, addHCount)
                names(flist)[length(flist)] <- "H"
            }
            else
                flist[["H"]] <- flist[["H"]] + addHCount
        }
        return(formulaListToString(flist))
    }))
}

getMFFragmentInfo <- function(spec, mfResult)
{
    if (mfResult$NoExplPeaks == 0 || mfResult$FormulasOfExplPeaks == "NA")
        return(data.table(mz = numeric(0), formula = character(0), score = numeric(0), PLIndex = numeric(0)))

    # format of FormulasOfExplPeaks: list of strings with mz1:formula1;mz2:formula2;...
    fi <- unlist(strsplit(mfResult$FormulasOfExplPeaks, "[;:]")) # split into list with subsequent m/z / formula pairs

    ret <- data.table(mz = as.numeric(fi[c(TRUE, FALSE)]),
                      formula = cleanFragFormulas(fi[c(FALSE, TRUE)]),
                      score = as.numeric(unlist(strsplit(mfResult$FragmenterScore_Values, ";"))))
    ret[, PLIndex := sapply(mz, function(omz) which.min(abs(omz - spec$mz)))]
    ret[, intensity := spec$intensity[PLIndex]]

    return(ret)
}

initMetFragCLCommand <- function(mfSettings, spec, mfBin, logFile)
{
    paramFile <- tempfile("parameters", fileext = ".txt")
    paramCon <- file(paramFile, "w")

    writeParam <- function(param, val) cat(sprintf("%s = %s\n", param, paste0(as.character(val), collapse = ",")), file = paramCon)

    for (param in names(mfSettings))
        writeParam(param, mfSettings[[param]])

    outFile <- tempfile("results", fileext = ".csv")
    writeParam("MetFragCandidateWriter", "CSV")
    writeParam("SampleName", basename(tools::file_path_sans_ext(outFile)))
    writeParam("ResultsPath", dirname(outFile))

    specFile <- tempfile("spectrum", fileext = ".txt")
    write.table(spec[, c("mz", "intensity")], specFile, sep = "\t", row.names = FALSE, col.names = FALSE)
    writeParam("PeakListPath", specFile)

    close(paramCon)

    return(list(command = "java", args = c("-jar", mfBin, paramFile), stderrFile = logFile, outFile = outFile))
}

generateMetFragRunData <- function(fGroups, MSPeakLists, mfSettings, topMost, identifiers, method, addTrivialNames)
{
    gNames <- names(fGroups)
    gTable <- groups(fGroups)
    gInfo <- groupInfo(fGroups)
    anaInfo <- analysisInfo(fGroups)
    pLists <- peakLists(MSPeakLists)

    ret <- sapply(gNames, function(grp)
    {
        if (!is.null(identifiers) && is.null(identifiers[[grp]]))
            return(NULL)

        ogind <- order(-gTable[[grp]])
        oanalyses <- anaInfo$analysis[ogind] # analyses ordered from highest to lowest intensity

        # filter out analyses without MS/MS
        oanalyses <- sapply(oanalyses, function(a) if (!is.null(pLists[[a]][[grp]]$MSMS)) a else "", USE.NAMES = FALSE)
        oanalyses <- oanalyses[oanalyses != ""]

        if (length(oanalyses) < 1)
            return(NULL)

        spec <- pLists[[oanalyses[[1]]]][[grp]]$MSMS
        mfSettings$IonizedPrecursorMass <- gInfo[grp, "mzs"]
        mfSettings$ExperimentalRetentionTimeValue <- gInfo[grp, "rts"] / 60

        if (!is.null(identifiers))
            mfSettings$PrecursorCompoundIDs <- identifiers[[grp]]

        hash <- makeHash(method, mfSettings, spec, topMost, if (method == "R") addTrivialNames else FALSE)

        return(list(hash = hash, analysis = oanalyses[[1]], gName = grp, spec = spec, mfSettings = mfSettings))
    }, simplify = FALSE)

    return(ret[!sapply(ret, is.null)])
}

processMFResults <- function(metf, analysis, spec, db, topMost, lfile = "")
{
    if (nrow(metf) > 0)
    {
        if (!is.null(topMost) && nrow(metf) > topMost)
            metf <- metf[seq_len(topMost)]

        # make character: needed for getMFFragmentInfo()
        # note: this column is not present in empty tables so can't do this with colClasses
        metf[, FragmenterScore_Values := as.character(FragmenterScore_Values)]

        # fill in fragment info
        # NOTE: double wrap in list to nest table
        if (!is.null(lfile))
            cat(sprintf("\n%s - Done! Processing frags...\n", date()), file = lfile, append = TRUE)
        for (r in seq_len(nrow(metf)))
            set(metf, r, "fragInfo", list(list(getMFFragmentInfo(spec, metf[r]))))
        if (!is.null(lfile))
            cat(sprintf("\n%s - Done!\n", date()), file = lfile, append = TRUE)

        # unify column names & filter unnecessary columns
        metf <- unifyMFNames(metf)

        metf[, analysis := analysis]
        metf[, database := db]
    }

    return(metf)
}

#' @details \code{generateCompoundsMetfrag} uses the \pkg{metfRag} package for
#'   compound identification (see \url{http://c-ruttkies.github.io/MetFrag/}).
#'   Several online compound databases such as
#'   \href{https://pubchem.ncbi.nlm.nih.gov/}{PubChem} and
#'   \href{http://www.chemspider.com/}{ChemSpider} may be chosen for retrieval
#'   of candidate structures. In addition, many options exist to score and
#'   filter resulting data, and it is highly suggested to optimize these to
#'   improve results. While MS/MS data is not mandatory, it will usually greatly
#'   improve candidate scoring. The \command{MetFrag} options \code{PeakList},
#'   \code{IonizedPrecursorMass} and \code{ExperimentalRetentionTimeValue} (in
#'   minutes) fields are automatically set from feature data.
#'
#' @param method Which method should be used for MetFrag execution: \code{"CL"}
#'   for \command{MetFragCL} and \code{"R"} for \command{MetFragR}. The former
#'   might be faster.
#' @param timeout Maximum time (in seconds) before a metFrag query for a feature
#'   group is stopped. Also see \code{timeoutRetries} argument.
#' @param timeoutRetries Maximum number of retries after reaching a timeout
#'   before completely skipping the metFrag query for a feature group. Also see
#'   \code{timeout} argument.
#' @param dbRelMzDev Relative mass deviation (in ppm) for database search. Sets
#'   the \option{DatabaseSearchRelativeMassDeviation} option.
#' @param fragRelMzDev Relative mass deviation (in ppm) for fragment matching.
#'   Sets the \option{FragmentPeakMatchRelativeMassDeviation} option.
#' @param fragAbsMzDev Absolute mass deviation (in Da) for fragment matching.
#'   Sets the \option{FragmentPeakMatchAbsoluteMassDeviation} option.
#' @param isPositive Set to \code{TRUE} for data measured with positive
#'   ionization. Sets the \option{IsPositiveIonMode} option.
#' @param database Compound database to use. Valid values are: \code{"PubChem"},
#'   \code{"ExtendedPubChem"}, \code{"ChemSpider"} and \code{"KEGG"}. The
#'   \code{"ExtendedPubChem"} is a PubChem database which includes number of
#'   patents and references information, which an be used for further scoring
#'   (see \code{scoreTypes} parameter.) In order to use the \code{ChemSpider}
#'   database the \code{chemSpiderToken} should be set. Sets the
#'   \option{MetFragDatabaseType} option.
#' @param chemSpiderToken A character string with the
#'   \href{http://www.chemspider.com/AboutServices.aspx}{ChemSpider security
#'   token} that should be set when the ChemSpider database is used. Sets the
#'   \option{ChemSpiderToken} option.
#' @param scoreTypes A character vector defining the scoring methods
#'   (\emph{e.g.} \code{"FragmenterScore"}, \code{"OfflineMetFusionScore"},
#'   \code{"RetentionTimeScore"}). Some methods require further options to be
#'   set. Additional scoring methods become available when the
#'   \option{"ExtendedPubChem"} (\emph{i.e.} \option{"PubChemNumberPatents"} and
#'   \option{"PubChemNumberPubMedReferences"}) or \option{"ChemSpider"}
#'   (\emph{i.e.} \option{"ChemSpiderReferenceCount"},
#'   \option{"ChemSpiderNumberExternalReferences"},
#'   \option{"ChemSpiderRSCCount"}, \option{"ChemSpiderNumberPubMedReferences"}
#'   and \option{"ChemSpiderDataSourceCount"}) database type is used. For all
#'   scoring types and more information refer to the \verb{Candidate Scores}
#'   section on the
#'   \href{http://c-ruttkies.github.io/MetFrag/projects/metfragr/}{MetFragR
#'   homepage}. Sets the \option{MetFragScoreTypes} option.
#' @param scoreWeights Numeric vector containing weights of the used scoring
#'   types. Order is the same as set in \code{scoreTypes}. Values are recycled
#'   if necessary. Sets the \option{MetFragScoreWeights} option.
#' @param preProcessingFilters,postProcessingFilters A character vector defining
#'   pre/post filters applied before/after fragmentation and scoring
#'   (\emph{e.g.} \code{"UnconnectedCompoundFilter"}, \code{"IsotopeFilter"},
#'   \code{"ElementExclusionFilter"}). Some methods require further options to
#'   be set. For all filters and more information refer to the \verb{Candidate
#'   Filters} section on the
#'   \href{http://c-ruttkies.github.io/MetFrag/projects/metfragr/}{MetFragR
#'   homepage}. Sets the \option{MetFragPreProcessingCandidateFilter} and
#'   \code{MetFragPostProcessingCandidateFilter} options.
#' @param maxCandidatesToStop If more than this number of candidate structures
#'   are found then processing will be aborted and no results this feature group
#'   will be reported. Low values increase the chance of missing data, whereas
#'   too high values will use too much computer resources and signficantly
#'   slowdown the process. Sets the \option{MaxCandidateLimitToStop} option.
#' @param addTrivialNames If \code{TRUE} and the PubChem database is used then
#'   trivial names will be added after compound search.
#' @param identifiers A \code{list} containing for each feature group a
#'   character vector with database identifiers that should be used to find
#'   candidates for a feature group (the list should be named by feature group
#'   names). If \code{NULL} all relevant candidates will be retrieved from the
#'   specified database. An example usage scenario is to obtain the list of
#'   candidate identifiers from a \code{\link{compounds}} object obtained with
#'   \code{\link{generateCompoundsSirius}} using the \code{\link{identifiers}}
#'   method. This way, only those candidates will be searched by MetFrag that
#'   were generated by SIRIUS+CSI:FingerID. Sets the
#'   \option{PrecursorCompoundIDs} option.
#' @param extraOpts A named \code{list} containing further settings to be passed
#'   to \code{\link[metfRag]{run.metfrag}}. See the
#'   \href{http://c-ruttkies.github.io/MetFrag/projects/metfragr/}{MetFragR} and
#'   \href{http://c-ruttkies.github.io/MetFrag/projects/metfragcl/}{MetFrag CL}
#'   homepages for all available options.
#'
#' @references \insertRef{Ruttkies2016}{patRoon}
#'
#' @rdname compound-generation
#' @export
generateCompoundsMetfrag <- function(fGroups, MSPeakLists, method = "CL", logPath = file.path("log", "metfrag"),
                                     timeout = 300, timeoutRetries = 2, errorRetries = 2, topMost = 100,
                                     dbRelMzDev = 5, fragRelMzDev = 5, fragAbsMzDev = 0.002, isPositive, adduct,
                                     database = "PubChem", chemSpiderToken = "",
                                     scoreTypes = c("FragmenterScore", "OfflineMetFusionScore"), scoreWeights = 1.0,
                                     preProcessingFilters = c("UnconnectedCompoundFilter","IsotopeFilter"),
                                     postProcessingFilters = c("InChIKeyFilter"),
                                     maxCandidatesToStop = 2500, addTrivialNames = TRUE,
                                     identifiers = NULL, extraOpts = NULL,
                                     maxProcAmount = getOption("patRoon.maxProcAmount"))
{
    # UNDONE: does addTrivialNames actually work?
    
    if (method == "R")
        checkPackage("metfRag", "c-ruttkies/MetFragR")
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertChoice(method, c("CL", "R"), add = ac)
    assertMultiProcArgs(logPath, maxProcAmount, add = ac)
    aapply(checkmate::assertNumber, . ~ timeout + dbRelMzDev + fragRelMzDev + fragAbsMzDev,
           lower = 0, finite = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertCount, . ~ timeoutRetries + errorRetries, fixed = list(add = ac))
    aapply(checkmate::assertCount, . ~ topMost + maxCandidatesToStop, positive = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertFlag, . ~ isPositive + addTrivialNames, fixed = list(add = ac))
    checkmate::assertInt(adduct, add = ac)
    aapply(checkmate::assertString, . ~ database + chemSpiderToken, fixed = list(add = ac))
    aapply(checkmate::assertCharacter, . ~ scoreTypes + preProcessingFilters + postProcessingFilters,
           any.missing = FALSE, fixed = list(add = ac))
    checkmate::assertNumeric(scoreWeights, lower = 0, finite = TRUE, any.missing = FALSE, min.len = 1, add = ac)
    aapply(checkmate::assertList, . ~ identifiers + extraOpts, any.missing = FALSE,
           names = "unique", null.ok = TRUE, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    anaInfo <- analysisInfo(fGroups)
    ftind <- groupFeatIndex(fGroups)
    gTable <- groups(fGroups)
    gNames <- names(fGroups)
    gCount <- length(fGroups)
    gInfo <- groupInfo(fGroups)
    pLists <- peakLists(MSPeakLists)

    mfSettings <- list(DatabaseSearchRelativeMassDeviation = dbRelMzDev,
                       FragmentPeakMatchRelativeMassDeviation = fragRelMzDev,
                       FragmentPeakMatchAbsoluteMassDeviation = fragAbsMzDev,
                       PrecursorIonMode = adduct, IsPositiveIonMode = isPositive,
                       MetFragDatabaseType = database, MetFragScoreTypes = scoreTypes,
                       MetFragScoreWeights = rep(scoreWeights, length.out = length(scoreTypes)),
                       MetFragPreProcessingCandidateFilter = preProcessingFilters,
                       MetFragPostProcessingCandidateFilter = postProcessingFilters,
                       MaxCandidateLimitToStop = maxCandidatesToStop)

    if (!is.null(chemSpiderToken) && nzchar(chemSpiderToken))
        mfSettings$ChemSpiderToken <- chemSpiderToken

    if (!is.null(extraOpts))
    {
        # add (or replace) any extra options
        mfSettings <- modifyList(mfSettings, extraOpts)
    }

    cacheDB <- openCacheDB()
    setHash <- makeHash(fGroups, pLists, method, mfSettings, topMost, identifiers, addTrivialNames)
    cachedSet <- loadCacheSet("identifyMetFrag", setHash, cacheDB)
    resultHashes <- vector("character", length(gNames))
    names(resultHashes) <- gNames

    printf("Identifying %d feature groups with metFrag...\n", gCount)

    runData <- generateMetFragRunData(fGroups, MSPeakLists, mfSettings, topMost, identifiers, method, addTrivialNames)

    cachedResults <- sapply(runData, function(rd)
    {
        resultHashes[rd$gName] <<- rd$hash
        metf <- NULL
        if (!is.null(cachedSet))
            metf <- cachedSet[[rd$hash]]
        if (is.null(metf))
            metf <- loadCacheData("identifyMetFrag", rd$hash, cacheDB)
        return(metf)
    }, simplify = FALSE)
    cachedResults <- cachedResults[!sapply(cachedResults, is.null)]

    runData <- runData[setdiff(names(runData), names(cachedResults))] # remove cached results

    if (length(runData) > 0)
    {
        if (method == "CL")
        {
            if (!is.null(logPath))
                mkdirp(logPath)
            mfBin <- path.expand(getOption("patRoon.path.metFragCL"))
            if (is.null(mfBin) || !nzchar(mfBin) || !file.exists(mfBin))
                stop("Please set the 'metFragCL' option with a (correct) path to the metFrag CL jar file. Example: options(patRoon.path.metFragCL = \"C:/MetFrag2.4.2-CL.jar\")")

            if (!nzchar(Sys.which("java")))
                stop("Please make sure that java is installed and its location is correctly set in PATH.")

            cmdQueue <- lapply(runData, function(rd)
            {
                logf <- if (!is.null(logPath)) file.path(logPath, paste0("mfcl-", rd$gName, ".txt")) else NULL
                return(c(rd[c("hash", "analysis", "gName", "spec")],
                         initMetFragCLCommand(rd$mfSettings, rd$spec, mfBin, logf)))
            })

            ret <- executeMultiProcess(cmdQueue, finishHandler = function(cmd, exitStatus, retries)
            {
                if (!is.null(cmd$stderrFile))
                    cat(sprintf("\n%s - Done with MF! Reading results...\n", date()), file = cmd$stderrFile, append = TRUE)
                metf <- fread(cmd$outFile, colClasses = c(Identifier = "character"))
                if (!is.null(cmd$stderrFile))
                    cat(sprintf("\n%s - Done! Processing results...\n", date()), file = cmd$stderrFile, append = TRUE)
                metf <- processMFResults(metf, cmd$analysis, cmd$spec, database, topMost, cmd$stderrFile)
                if (!is.null(cmd$stderrFile))
                    cat(sprintf("\n%s - Done! Caching results...\n", date()), file = cmd$stderrFile, append = TRUE)
                saveCacheData("identifyMetFrag", metf, cmd$hash, cacheDB)
                if (!is.null(cmd$stderrFile))
                    cat(sprintf("\n%s - Done!\n", date()), file = cmd$stderrFile, append = TRUE)

                return(metf)
            }, timeoutHandler = function(cmd, retries)
            {
                if (retries >= timeoutRetries)
                {
                    warning(sprintf("Could not run MetFrag for %s: timeout", cmd$gName))
                    return(FALSE)
                }
                warning(sprintf("Restarting timed out MetFrag command for %s (retry %d/%d)",
                                cmd$gName, retries+1, errorRetries))
                return(TRUE)
            }, errorHandler = function(cmd, exitStatus, retries)
            {
                if (exitStatus <= 6) # some error thrown by MF
                {
                    if (retries >= errorRetries)
                    {
                        warning(sprintf("Could not run MetFrag for %s - exit code: %d", cmd$gName, exitStatus))
                        return(FALSE)
                    }
                    warning(sprintf("Restarting failed MetFrag command for %s - exit: %d (retry %d/%d)",
                                    cmd$gName, exitStatus, retries+1, errorRetries))
                    return(TRUE)
                }
                
                # some other error (e.g. java not present)
                stop(sprintf("Fatal: Failed to execute MetFragCL for %s - exit code: %d", cmd$gName, exitStatus))
            }, maxProcAmount = maxProcAmount, procTimeout = timeout, delayBetweenProc = 200)
        }
        else
        {
            prog <- txtProgressBar(0, gCount, style = 3)

            ret <- lapply(runData, function(rd)
            {
                rd$mfSettings$PeakList <- as.matrix(rd$spec[, c("mz", "intensity")])
                metf <- metfRag::run.metfrag(rd$mfSettings)
                jgc() # hopefully reduce some memory usage

                if (nrow(metf) > 0)
                {
                    if (addTrivialNames && rd$mfSettings$MetFragDatabaseType %in% c("PubChem", "ExtendedPubChem"))
                    {
                        # fetching trivial names may sometimes fail with connection error, just ignore this for now
                        tryCatch(metf <<- metfRag::add.trivialname.pubchem(metf), error = function(e) metf$TrivialName <<- NA)
                    }

                    metf <- unFactorDF(metf)
                    metf <- as.data.table(metf)

                    metf <- processMFResults(metf, rd$analysis, rd$spec, database, topMost)

                    # BUG: metfRag seems to give second duplicate results where only NoExplPeaks may differ and have an incorrect value.
                    # for now, just remove all duplicates and re-assign NoExplPeaks
                    metf <- metf[!duplicated(identifier)]
                    metf[, explainedPeaks := sapply(fragInfo, nrow)]
                }

                saveCacheData("identifyMetFrag", metf, rd$hash, cacheDB)

                setTxtProgressBar(prog, match(rd$gName, gNames))

                return(metf)
            })

            setTxtProgressBar(prog, gCount)
            close(prog)
        }
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
        ret <- ret[sapply(ret, function(r) !is.null(r) && nrow(r) > 0, USE.NAMES = FALSE)]

    ngrp <- length(ret)
    printf("Loaded %d compounds from %d features (%.2f%%).\n", sum(unlist(lapply(ret, nrow))),
           ngrp, ngrp * 100 / gCount)

    if (is.null(cachedSet))
        saveCacheSet("identifyMetFrag", resultHashes[resultHashes != ""], setHash, cacheDB)

    closeCacheDB(cacheDB)

    return(compounds(compounds = ret, algorithm = "MetFrag"))
}
