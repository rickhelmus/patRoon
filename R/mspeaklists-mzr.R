#' @include main.R
#' @include mspeaklists.R
#' @include utils-mzr.R
NULL

# use mzR to generate MS peaklists.
# limitations compared to DA: no bg subtraction, no isotope information

#' @details \code{generateMSPeakListsMzR} uses the \pkg{mzR} package to
#'   extract MS peak lists. For this analyses should be either in \file{.mzXML}
#'   or \file{.mzML} format. This function averages multiple spectra over a
#'   chromatgraphic peak to improve accuracy.
#'
#' @param precursorMzWindow The \emph{m/z} window (in Da) to find MS/MS spectra
#'   of a precursor. This is typically used for Data-Dependent like MS/MS data
#'   and should correspond to the isolation \emph{m/z} window (\emph{i.e.} +/-
#'   the precursor \emph{m/z}) that was used to collect the data. For
#'   Data-Independent MS/MS experiments, where precursor ions are not isolated
#'   prior to fragmentation (\emph{e.g.} bbCID, MSe, all-ion, ...) the value
#'   should be \code{NULL}.
#'
#' @rdname MSPeakLists-generation
#' @export
setMethod("generateMSPeakListsMzR", "featureGroups", function(fGroups, maxMSRtWindow = 5,
                                                              precursorMzWindow = 4, topMost = NULL,
                                                              avgFeatParams = getDefAvgPListParams(),
                                                              avgFGroupParams = getDefAvgPListParams())
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertNumber(maxMSRtWindow, lower = 1, finite = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertNumber(precursorMzWindow, lower = 0, finite = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertCount(topMost, positive = TRUE, null.ok = TRUE, add = ac)
    assertAvgPListParams(avgFeatParams, add = ac)
    assertAvgPListParams(avgFGroupParams, add = ac)
    checkmate::reportAssertions(ac)

    ftindex <- groupFeatIndex(fGroups)
    fTable <- featureTable(fGroups)
    anaInfo <- analysisInfo(fGroups)
    gTable <- groupTable(fGroups)
    gCount <- length(fGroups)
    gNames <- names(fGroups)
    anaCount <- nrow(anaInfo)

    if (gCount == 0)
        return(MSPeakLists(algorithm = "mzr"))

    cacheDB <- openCacheDBScope()
    setHash <- makeHash(fGroups, maxMSRtWindow, precursorMzWindow, topMost, avgFeatParams)
    cachedSet <- loadCacheSet("MSPeakListsMzR", setHash, cacheDB)
    resultHashes <- vector("character", anaCount * gCount)
    resultHashCount <- 0

    avgFeatParamsMS <- avgFeatParamsMSMS <-
        avgFeatParams[setdiff(names(avgFeatParams), c("pruneMissingPrecursorMS", "retainPrecursorMSMS"))]
    avgFeatParamsMS$retainPrecursor <- TRUE;
    avgFeatParamsMS$pruneMissingPrecursor <- avgFeatParams$pruneMissingPrecursorMS
    avgFeatParamsMSMS$pruneMissingPrecursor <- FALSE
    avgFeatParamsMSMS$retainPrecursor <- avgFeatParams$retainPrecursorMSMS

    # structure: [[analysis]][[fGroup]][[MSType]][[MSPeak]]
    plists <- list()
    metadata <- list()

    # if topMost is specified, make list (topAna) of topMost intense analyses
    if (!is.null(topMost) && topMost > 0)
        topAna <- sapply(seq_along(gTable), function(grpi) anaInfo$analysis[order(gTable[[grpi]], decreasing = TRUE)[seq_len(topMost)]])

    for (anai in seq_len(anaCount))
    {
        ana <- anaInfo$analysis[anai]
        fp <- getMzMLOrMzXMLAnalysisPath(ana, anaInfo$path[anai])
        spectra <- NULL

        baseHash <- makeHash(ana, maxMSRtWindow, precursorMzWindow, topMost, avgFeatParams)

        printf("Loading all MS peak lists for %d feature groups in analysis '%s'...\n", gCount, ana)
        prog <- openProgBar(0, gCount)

        for (grpi in seq_along(ftindex))
        {
            if (!is.null(topMost) && !ana %in% topAna[[grpi]])
                next # not intense enough

            grp <- gNames[grpi]

            fti <- ftindex[[grpi]][anai]
            if (fti == 0)
                next
            ft <- fTable[[ana]][fti, ]

            hash <- makeHash(baseHash, ft)
            resultHashCount <- resultHashCount + 1
            resultHashes[resultHashCount] <- hash

            results <- NULL
            if (!is.null(cachedSet))
                results <- cachedSet[[hash]]
            if (is.null(results))
                results <- loadCacheData("MSPeakListsMzR", hash, cacheDB)

            if (is.null(results))
            {
                results <- list(plists = list(), metatadata = list())

                rtRange <- c(ft$retmin, ft$retmax)
                if (!is.null(maxMSRtWindow) && diff(rtRange) > maxMSRtWindow*2)
                    rtRange <- c(max(rtRange[1], ft$ret - maxMSRtWindow), min(rtRange[2], ft$ret + maxMSRtWindow))

                if (is.null(spectra))
                    spectra <- loadSpectra(fp, verbose = FALSE)

                # NOTE: precursor is set here only for precursor assignment,
                # keeping precursorMzWindow unset for the header makes sure that
                # no spectra selection is made.
                results$metadata$MS <- getSpectraHeader(spectra, rtRange, 1, NULL, NULL)
                results$plists$MS <- do.call(averageSpectraMZR,
                                             c(list(spectra = spectra, hd = results$metadata$MS, precursor = ft$mz), avgFeatParamsMS))

                hdMSMS <- getSpectraHeader(spectra, rtRange, 2, ft$mz, precursorMzWindow)
                MSMS <- do.call(averageSpectraMZR, c(list(spectra = spectra, hd = hdMSMS, precursor = ft$mz),
                                                     avgFeatParamsMSMS))
                if (nrow(MSMS) > 0)
                {
                    results$plists$MSMS <- MSMS
                    results$metadata$MSMS <- hdMSMS
                }

                saveCacheData("MSPeakListsMzR", results, hash, cacheDB)
            }

            plists[[ana]][[grp]] <- results$plists
            metadata[[ana]][[grp]] <- results$metadata

            setTxtProgressBar(prog, grpi)
        }

        setTxtProgressBar(prog, gCount)
        close(prog)
    }

    if (is.null(cachedSet))
        saveCacheSet("MSPeakListsMzR", resultHashes[seq_len(resultHashCount)], setHash, cacheDB)

    return(MSPeakLists(peakLists = plists, metadata = metadata, avgPeakListArgs = avgFGroupParams,
                       origFGNames = gNames, algorithm = "mzr"))
})

setMethod("generateMSPeakListsMzR", "featureGroupsSet", function(fGroups, ...,
                                                                 avgSetParams = getDefAvgPListParams())
{
    generateMSPeakListsSet(fGroups, generateMSPeakListsMzR, ..., avgSetParams = avgSetParams)
})
