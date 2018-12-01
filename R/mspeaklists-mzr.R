#' @include main.R
#' @include mspeaklists.R
#' @include utils-mzr.R
NULL

# use mzR to generate MS peaklists.
# limitations compared to DA: no bg subtraction, no isotope information

#' @details \code{generateMSPeakListsMzR} uses the \pkg{\link{mzR}} package to
#'   extract MS peak lists. For this analyses should be either in \file{.mzXML}
#'   or \file{.mzML} format. This function can optionally average multiple
#'   spectra over a chromatgraphic peak to improve accuracy.
#'
#' @param avgMzWindow \emph{m/z} window (in Da) used for clustering \emph{m/z}
#'   values when spectra are averaged. Too small windows will prevent clustering
#'   \emph{m/z} values (thus erroneously creating 'extra' values), whereas too
#'   big windows may cluster unrelated \emph{m/z} values from different or even
#'   the same spectrum together.
#' @param avgTopPeaks Only retain a maximum number of \code{avgTopPeaks} MS
#'   peaks when generating averaged spectra. Lowering this number may exclude
#'   more irrelevant (noisy) MS peaks and decrease processing time, whereas
#'   higher values may avoid excluding lower intense MS peaks that may still be
#'   of interest.
#' @param avgMinIntensity MS peaks with intensities below this value will be
#'   completely excluded (overrides \code{avgTopPeaks}).
#' @param avgMassFun Function that is used to calculate average \emph{m/z}
#'   values.
#' @param avgMethod Method used for producing averaged MS spectra. Valid values
#'   are \code{"hclust"}, used for hierarchical clustering, and
#'   \code{"distance"}, to use the between peak distance. The latter method
#'   significantly reduces processing time and memory requirements, at the
#'   (potential) cost of slightly reducing accuracy.
#' @param precursorMzWindow The \emph{m/z} window (in Da) to find MS/MS spectra
#'   of a precursor. This is typically used for Data-Dependent like MS/MS data
#'   and should correspond to the isolation \emph{m/z} width that was used to
#'   collect the data. For Data-Independent MS/MS experiments, where precursor
#'   ions are not isolated prior to fragmentation (\emph{e.g.} bbCID, MSe,
#'   all-ion, ...) the value should be \code{NULL}.
#'
#' @references Averaging of mass spectra algorithms used by
#'   \code{generateMSPeakListsMzR} are based on the
#'   \href{https://github.com/zeehio/msProcess}{msProcess} R package (now
#'   archived on CRAN). \cr\cr
#'   \addCitations{mzR}{1} \cr\cr
#'   \addCitations{mzR}{2} \cr\cr
#'   \addCitations{mzR}{3} \cr\cr
#'   \addCitations{mzR}{4} \cr\cr
#'   \addCitations{mzR}{5}
#'
#' @rdname MSPeakLists-generation
#' @export
generateMSPeakListsMzR <- function(fGroups, maxRtMSWidth = 20, precursorMzWindow = 8, topMost = NULL,
                                   avgFeatParams = getDefAvgPListParams(),
                                   avgFGroupParams = getDefAvgPListParams())
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    checkmate::assertNumber(maxRtMSWidth, lower = 1, finite = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertNumber(precursorMzWindow, lower = 0, finite = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertCount(topMost, positive = TRUE, null.ok = TRUE, add = ac)
    assertAvgPListParams(avgFeatParams, add = ac)
    assertAvgPListParams(avgFGroupParams, add = ac)
    checkmate::reportAssertions(ac)
    
    ftindex <- groupFeatIndex(fGroups)
    fTable <- featureTable(fGroups)
    anaInfo <- analysisInfo(fGroups)
    gTable <- groups(fGroups)
    gCount <- length(fGroups)
    gNames <- names(fGroups)
    anaCount <- nrow(anaInfo)
    
    if (gCount == 0)
        return(MSPeakLists(algorithm = "mzR"))

    cacheDB <- openCacheDBScope()
    setHash <- makeHash(fGroups, maxRtMSWidth, precursorMzWindow, topMost, avgFeatParams)
    cachedSet <- loadCacheSet("MSPeakListsMzR", setHash, cacheDB)
    resultHashes <- vector("character", anaCount * gCount)
    resultHashCount <- 0

    # structure: [[analysis]][[fGroup]][[MSType]][[MSPeak]]
    plists <- list()

    # if topMost is specified, make list (topAna) of topMost intense analyses
    if (!is.null(topMost) && topMost > 0)
        topAna <- sapply(seq_along(gTable), function(grpi) anaInfo$analysis[order(gTable[[grpi]], decreasing = TRUE)[seq_len(topMost)]])

    for (anai in seq_len(anaCount))
    {
        ana <- anaInfo$analysis[anai]
        fp <- getMzMLOrMzXMLAnalysisPath(ana, anaInfo$path[anai])
        spectra <- NULL

        baseHash <- makeHash(ana, maxRtMSWidth, precursorMzWindow, topMost, avgFeatParams)

        printf("Loading all MS peak lists for %d feature groups in analysis '%s'...\n", gCount, ana)
        prog <- txtProgressBar(0, gCount, style = 3)

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
                results <- list()

                rtRange <- c(ft$retmin, ft$retmax)
                if (!is.null(maxRtMSWidth) && diff(rtRange) > maxRtMSWidth)
                    rtRange <- c(max(rtRange[1], ft$ret - maxRtMSWidth/2), min(rtRange[2], ft$ret + maxRtMSWidth/2))

                if (is.null(spectra))
                    spectra <- loadSpectra(fp, verbose = FALSE)

                results$MS <- do.call(averageSpectraMZR, c(list(spectra, rtRange), avgFeatParams))

                stopifnot(!is.null(ft$mz))
                MSMS <- do.call(averageSpectraMZR, c(list(spectra = spectra, rtRange = rtRange, MSLevel = 2,
                                                          precursor = ft$mz, precursorMzWindow = precursorMzWindow),
                                                     avgFeatParams))
                if (nrow(MSMS) > 0)
                    results$MSMS <- MSMS

                saveCacheData("MSPeakListsMzR", results, hash, cacheDB)
            }

            plists[[ana]][[grp]] <- results

            setTxtProgressBar(prog, grpi)
        }

        setTxtProgressBar(prog, gCount)
        close(prog)
    }

    if (is.null(cachedSet))
        saveCacheSet("MSPeakListsMzR", resultHashes[seq_len(resultHashCount)], setHash, cacheDB)
    
    return(MSPeakLists(peakLists = plists, avgPeakListArgs = avgFGroupParams, algorithm = "mzR"))
}
