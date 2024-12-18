# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include features.R
#' @include main.R
#' @include utils-bruker.R
NULL

# NOTE: No coverage calculation for Bruker tools as they cannot be run on CI
# nocov start

#' @rdname features-class
#' @export
featuresBruker <- setClass("featuresBruker", contains = "features")

setMethod("initialize", "featuresBruker",
          function(.Object, ...) callNextMethod(.Object, algorithm = "bruker_da", ...))

#' Find features using Bruker DataAnalysis
#'
#' Uses the 'Find Molecular Features' (FMF) algorithm of Bruker DataAnalysis vendor software to find features.
#'
#' @templateVar algo Bruker
#' @templateVar do automatically find features
#' @templateVar generic findFeatures
#' @templateVar algoParam bruker
#' @template algo_generator
#'
#' @inheritParams findFeatures
#'
#' @details The resulting 'compounds' are transferred from DataAnalysis and stored as features.
#'
#'   This algorithm only works with Bruker data files (\code{.d} extension) and requires Bruker DataAnalysis
#'   and the \pkg{RDCOMClient} package to be installed. Furthermore, DataAnalysis combines multiple related masses in a
#'   feature (\emph{e.g.} isotopes, adducts) but does not report the actual (monoisotopic) mass of the feature.
#'   Therefore, it is simply assumed that the feature mass equals that of the highest intensity mass peak.
#'
#' @param doFMF Run the 'Find Molecular Features' algorithm before loading compounds. Valid options are: \code{"auto"}
#'   (run FMF automatically if current results indicate it is necessary) and \code{"force"} (run FMF \emph{always}, even
#'   if cached results exist). Note that checks done if \code{doFMF="auto"} are fairly simplistic, hence set
#'   \code{doFMF="force"} if feature data needs to be updated.
#' @param startRange,endRange Start/End retention range (seconds) from which to collect features. A 0 (zero) for
#'   \code{endRange} marks the end of the analysis.
#'
#' @template dasaveclose-args
#' @template DA-restart-note
#' 
#' @inherit findFeatures return
#'
#' @export
findFeaturesBruker <- function(analysisInfo, doFMF = "auto", startRange = 0, endRange = 0,
                               save = TRUE, close = save, verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, "bruker", add = ac)
    checkmate::assertChoice(doFMF, c("auto", "force"), add = ac)
    checkmate::assertNumber(startRange, lower = 0, finite = TRUE, add = ac)
    checkmate::assertNumber(endRange, lower = 0, finite = TRUE, add = ac)
    assertDACloseSaveArgs(close, save, add = ac)
    checkmate::reportAssertions(ac)

    ret <- featuresBruker(analysisInfo = analysisInfo)

    DA <- getDAApplication()
    hideDAInScope()

    ret@features = sapply(seq_len(nrow(analysisInfo)),
                          function(i) getDAFeatures(DA, analysisInfo$analysis[i], analysisInfo$path[i], doFMF,
                                                    startRange, endRange, close, save, verbose), simplify = FALSE)
    names(ret@features) <- analysisInfo$analysis

    return(ret)
}

getDAFeatures <- function(DA, analysis, path, doFMF, startRange, endRange, close, save, verbose)
{
    printf("Loading bruker features for analysis '%s'...\n", analysis)

    hash <- makeHash(analysis, path, startRange, endRange) # UNDONE: better hash?
    ret <- if (doFMF == "force") NULL else loadCacheData("featuresBruker", hash)

    if (is.null(ret))
    {
        ind <- getDAFileIndex(DA, analysis, path)
        if (ind == -1)
            stop(sprintf("Failed to open analysis %s from path %s!", analysis, path))

        cmpds <- DA[["Analyses"]][[ind]][["Compounds"]]
        ccount <- cmpds[["Count"]]

        execFMF <- doFMF == "force" || ccount == 0

        # ensure that all DA compounds are from FMF
        if (!execFMF)
            execFMF <- any(sapply(seq_len(ccount), function(ci) cmpds[[ci]][["SeparationType"]] != "MolFeature"))

        if (execFMF)
        {
            if (verbose)
                printf("Running Find Molecular Features (FMF)... ")

            DA[["Analyses"]][[ind]][["Compounds"]]$Clear()
            DA[["Analyses"]][[ind]]$ClearChromatogramRangeSelections()
            DA[["Analyses"]][[ind]]$ClearResults()
            if (endRange > 0)
                DA[["Analyses"]][[ind]]$AddChromatogramRangeSelection(startRange, endRange, 0, 0)
            DA[["Analyses"]][[ind]]$FindMolecularFeatures()
            ccount <- cmpds[["Count"]]

            if (verbose)
                cat("Done!\n")
        }
        else if (verbose)
            printf("Skipping FMF for analysis '%s'\n", analysis)

        if (verbose && ccount > 0)
        {
            printf("Loading %d features from DataAnalysis...\n", ccount)
            prog <- openProgBar(0, ccount)
        }

        dt <- data.table(ID=numeric(ccount), ret=numeric(ccount), mz=numeric(ccount), intensity=numeric(ccount), peak_score=numeric(ccount),
                         retmin=numeric(ccount), retmax=numeric(ccount), mzmin=numeric(ccount), mzmax=numeric(ccount))
        dtCount <- 0

        for (i in seq_len(ccount))
        {
            cmpd <- cmpds[[i]]
            stopifnot(cmpd[["SeparationType"]] == "MolFeature")
            {
                maxintmz <- getDAMaxIntMZAndFWHM(cmpd[[1]])

                # estimate min/max mz from fwhm (https://en.wikipedia.org/wiki/Full_width_at_half_maximum)
                s <- maxintmz$fwhm / 2.355
                mzmn <- maxintmz$mz - 2*s
                mzmx <- maxintmz$mz + 2*s

                dtCount <- dtCount + 1
                dt[dtCount, c("ID", "ret", "mz", "intensity", "area", "peak_score", "retmin", "retmax", "mzmin", "mzmax") :=
                       .(i, cmpd[["RetentionTime"]], maxintmz$mz, cmpd[[1]][["MaximumIntensity"]], cmpd[["Area"]],
                         cmpd[["Area"]] / cmpd[["Intensity"]], cmpd[["RetentionTimeStart"]], cmpd[["RetentionTimeEnd"]],
                         mzmn, mzmx)]
            }

            if (verbose && (i %% 10) == 0)
                setTxtProgressBar(prog, i) # update every 10 iterations
        }

        if (verbose && ccount > 0)
        {
            setTxtProgressBar(prog, ccount)
            close(prog)
        }

        ret <- dt[seq_len(dtCount)]
        saveCacheData("featuresBruker", ret, hash)

        closeSaveDAFile(DA, ind, close, save)
    }

    if (verbose)
        cat("... Done!\n")

    return(ret)
}

# nocov end
