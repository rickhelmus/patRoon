#' @include features.R
#' @include main.R
#' @include utils-bruker.R
NULL

#' @rdname feature-finding
#' @export
featuresBruker <- setClass("featuresBruker", contains = "features")

#' @details \code{findFeaturesBruker} uses the 'Find Molecular Features' (FMF)
#'   algorithm of Bruker DataAnalysis vendor software to find features. The
#'   resulting 'compounds' are then transferred from DataAnalysis and stored as
#'   features.
#'
#' @note \code{findFeaturesBruker} only works with Bruker data files (\code{.d}
#'   extension) and requires Bruker DataAnalysis and the
#'   \pkg{RDCOMClient} package to be installed.
#'
#' @param doFMF Run the 'Find Molecular Features' algorithm before loading
#'   compounds. Valid options are: \code{TRUE} (always), \code{FALSE} (never)
#'   and \code{"auto"} (automatically if no compounds present).
#' @param startRange,endRange Start/End retention range (seconds) from which to
#'   collect features. A 0 (zero) for \code{endRange} marks the end of the
#'   analysis.
#' @rdname feature-finding
#' @export
findFeaturesBruker <- function(analysisInfo, doFMF = "auto", startRange = 0, endRange = 0, verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    assertAnalysisInfo(analysisInfo, "d", add = ac)
    checkmate::assertChoice(as.character(doFMF), c("auto", FALSE, TRUE), add = ac)
    checkmate::assertNumber(startRange, lower = 0, finite = TRUE, add = ac)
    checkmate::assertNumber(endRange, lower = 0, finite = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    ret <- featuresBruker(analysisInfo = analysisInfo)

    DA <- getDAApplication()
    hideDAInScope()

    ret@features = sapply(seq_len(nrow(analysisInfo)),
                          function(i) getDAFeatures(DA, analysisInfo$analysis[i], analysisInfo$path[i], doFMF,
                                                    startRange, endRange, verbose), simplify = FALSE)
    names(ret@features) <- analysisInfo$analysis

    return(ret)
}

getDAFeatures <- function(DA, analysis, path, doFMF, startRange, endRange, verbose)
{
    printf("Loading bruker features for analysis '%s'...\n", analysis)

    hash <- makeHash(analysis, path, startRange, endRange) # UNDONE: better hash?
    ret <- loadCacheData("featuresBruker", hash)

    if (is.null(ret))
    {
        ind <- getDAFileIndex(DA, analysis, path)

        if (ind == -1)
            return()

        cmpds <- DA[["Analyses"]][[ind]][["Compounds"]]
        ccount <- cmpds[["Count"]]

        if (ccount < 1)
        {
            if (doFMF == "auto")
                doFMF = TRUE
            else if (doFMF == FALSE)
                return()
        }

        if (doFMF == "auto")
        {
            # see if we need to do FMF
            doFMF <- TRUE
            for (i in 1:ccount)
            {
                if (cmpds[[i]][["SeparationType"]] == "MolFeature")
                {
                    doFMF <- FALSE
                    break
                }
            }
        }

        if (doFMF == TRUE)
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

        if (verbose)
        {
            printf("Loading %d features from DataAnalysis...\n", ccount, analysis)
            prog <- txtProgressBar(0, ccount, style=3)
        }

        dt <- data.table(ID=numeric(ccount), ret=numeric(ccount), mz=numeric(ccount), intensity=numeric(ccount), peak_score=numeric(ccount),
                         retmin=numeric(ccount), retmax=numeric(ccount), mzmin=numeric(ccount), mzmax=numeric(ccount))
        dtCount <- 0

        for (i in seq_len(ccount))
        {
            cmpd <- cmpds[[i]]
            if (cmpd[["SeparationType"]] == "MolFeature")
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

        if (verbose)
        {
            setTxtProgressBar(prog, ccount)
            close(prog)
        }

        ret <- dt[1:dtCount]
        saveCacheData("featuresBruker", ret, hash)
    }

    if (verbose)
        cat("... Done!\n")

    return(ret)
}
