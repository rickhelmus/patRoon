#' @include features.R
#' @include main.R
NULL

#' @rdname features-class
#' @export
featuresOpenMS <- setClass("featuresOpenMS", contains = "features")

setMethod("initialize", "featuresOpenMS",
          function(.Object, ...) callNextMethod(.Object, algorithm = "openms", ...))


#' @details \code{findFeaturesOpenMS} uses the
#'   \href{https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_FeatureFinderMetabo.html}{FeatureFinderMetabo}
#'    TOPP tool (see \url{http://www.openms.de}).
#'
#' @note The file format of analyses for \code{findFeaturesOpenMS} must be
#'   \file{mzML}. This functionality has been tested with OpenMS version >= 2.0.
#'   Please make sure it is installed and its binaries are added to the PATH
#'   environment variable or the \code{patRoon.path.OpenMS} option is set.
#'
#' @param noiseThrInt Noise intensity threshold. Sets
#'   \code{algorithm:common:noise_threshold_int} option.
#' @param chromSNR Minimum S/N of a mass trace. Sets
#'   \code{algorithm:common:chrom_peak_snr} option.
#' @param chromFWHM Expected chromatographic peak width (in seconds). Sets
#'   \code{algorithm:common:chrom_fwhm} option.
#' @param mzPPM Allowed mass deviation (ppm) for trace detection. Sets
#'   \code{algorithm:mtd:mass_error_ppm}.
#' @param reEstimateMTSD If \code{TRUE} then enables dynamic re-estimation of
#'   m/z variance during mass trace collection stage. Sets
#'   \code{algorithm:mtd:reestimate_mt_sd}.
#' @param traceTermCriterion,traceTermOutliers,minSampleRate Termination
#'   criterion for the extension of mass traces. See
#'   \href{https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_FeatureFinderMetabo.html}{FeatureFinderMetabo}.
#'    Sets the \code{algorithm:mtd:trace_termination_criterion},
#'   \code{algorithm:mtd:trace_termination_outliers} and
#'   \code{algorithm:mtd:min_sample_rate} options, respectively.
#' @param minTraceLength,maxTraceLength Minimum/Maximum length of mass trace
#'   (seconds). Set negative value for maxlength to disable maximum. Sets
#'   \code{algorithm:mtd:min_trace_length} and
#'   \code{algorithm:mtd:min_trace_length}, respectively.
#' @param widthFiltering,minFWHM,maxFWHM Enable filtering of unlikely peak
#'   widths. See
#'   \href{https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_FeatureFinderMetabo.html}{FeatureFinderMetabo}.
#'    Sets \code{algorithm:epd:width_filtering}, \code{algorithm:epd:min_fwhm}
#'   and \code{algorithm:epd:max_fwhm}, respectively.
#' @param traceSNRFiltering If \code{TRUE} then apply post-filtering by
#'   signal-to-noise ratio after smoothing. Sets the
#'   \code{algorithm:epd:masstrace_snr_filtering} option.
#' @param localRTRange,localMZRange Retention/MZ range where to look for
#'   coeluting/isotopic mass traces. Sets the
#'   \code{algorithm:ffm:local_rt_range} and \code{algorithm:ffm:local_mz_range}
#'   options, respectively.
#' @param isotopeFilteringModel Remove/score candidate assemblies based on
#'   isotope intensities. See
#'   \href{https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_FeatureFinderMetabo.html}{FeatureFinderMetabo}.
#'    Sets the \code{algorithm:ffm:isotope_filtering_model} option.
#' @param MZScoring13C Use the 13C isotope as the expected shift for isotope
#'   mass traces. See
#'   \href{https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_FeatureFinderMetabo.html}{FeatureFinderMetabo}.
#'    Sets \code{algorithm:ffm:mz_scoring_13C}.
#' @param useSmoothedInts If \code{TRUE} then use LOWESS intensities instead of
#'   raw intensities. Sets the \code{algorithm:ffm:use_smoothed_intensities}
#'   option.
#' @param extraOpts Named \code{list} containing extra options that will be
#'   passed to \command{FeatureFinderMetabo}. Any options specified here will
#'   override any of the above. Example:
#'   \code{extraOpts=list("-algorithm:common:noise_threshold_int"=1000)}
#'   (corresponds to setting \code{noiseThrInt=1000}). Set to \code{NULL} to
#'   ignore.
#' @param intSearchRTWindow Retention time window (in seconds, +/- feature
#'   retention time) that is used to find the closest data point to the
#'   retention time to obtain the intensity of a feature (this is needed since
#'   OpenMS does not provide this data).
#'
#' @template refs-openms
#'
#' @rdname feature-finding
#' @export
findFeaturesOpenMS <- function(analysisInfo, noiseThrInt = 1000, chromSNR = 3, chromFWHM = 5, mzPPM = 10, reEstimateMTSD = TRUE,
                               traceTermCriterion = "sample_rate", traceTermOutliers = 5, minSampleRate = 0.5,
                               minTraceLength = 3, maxTraceLength = -1, widthFiltering = "fixed", minFWHM = 3,
                               maxFWHM = 60, traceSNRFiltering = FALSE, localRTRange = 10, localMZRange = 6.5,
                               isotopeFilteringModel = "metabolites (5% RMS)", MZScoring13C = FALSE, useSmoothedInts = TRUE,
                               extraOpts = NULL, intSearchRTWindow = 3, verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, "mzML", add = ac)
    aapply(checkmate::assertNumber, . ~ noiseThrInt + chromSNR + chromFWHM + mzPPM + minSampleRate +
               minTraceLength + minFWHM + maxFWHM + localRTRange + localMZRange + intSearchRTWindow,
           lower = 0, finite = TRUE, fixed = list(add = ac))
    checkmate::assertChoice(traceTermCriterion, c("sample_rate", "outlier"))
    checkmate::assertCount(traceTermOutliers, add = ac)
    checkmate::assertNumber(maxTraceLength, finite = TRUE, add = ac)
    checkmate::assertChoice(widthFiltering, c("fixed", "off", "auto"))
    checkmate::assertChoice(isotopeFilteringModel, c("metabolites (2% RMS)", "metabolites (5% RMS)",
                                                     "peptides", "none"), add = ac)
    aapply(checkmate::assertFlag, . ~ reEstimateMTSD + traceSNRFiltering + MZScoring13C + useSmoothedInts, fixed = list(add = ac))
    checkmate::assertList(extraOpts, any.missing = FALSE, names = "unique", null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    if (verbose)
        printf("Finding features with OpenMS for %d analyses ...\n", nrow(analysisInfo))

    params <- list(noiseThrInt, chromSNR, chromFWHM, mzPPM, reEstimateMTSD,
                   traceTermCriterion, traceTermOutliers, minSampleRate,
                   minTraceLength, maxTraceLength, widthFiltering, minFWHM,
                   maxFWHM, traceSNRFiltering, localRTRange, localMZRange,
                   isotopeFilteringModel, MZScoring13C, useSmoothedInts, extraOpts)
    paramsHash <- makeHash(params)

    cmdQueue <- lapply(seq_len(nrow(analysisInfo)), function(anai)
    {
        dfile <- getMzMLAnalysisPath(analysisInfo$analysis[anai], analysisInfo$path[anai])
        hash <- makeHash(makeFileHash(dfile), paramsHash)
        logf <- paste0("ffm-", analysisInfo$analysis[anai], ".txt")
        return(list(hash = hash, dataFile = dfile, logFile = logf))
    })
    names(cmdQueue) <- analysisInfo$analysis

    fList <- list()
    if (length(cmdQueue) > 0)
    {
        fList <- executeMultiProcess(cmdQueue, function(cmd, ...)
        {
            fts <- patRoon:::importFeatureXML(cmd$featFile)
            return(fts)
        }, prepareHandler = function(cmd)
        {
            ffile <- tempfile(fileext = ".featureXML")
            cmd <- do.call(patRoon:::getOpenMSFFCommand, c(list(cmd$dataFile, ffile), params))
            return(c(list(featFile = ffile), cmd))
        }, showProgress = verbose, logSubDir = "openms", cacheName = "featuresOpenMS")
        
        # load intensities afterwards: we want to use the cache if possible,
        # which wouldn't be possible if future MP is used.
        fList <- sapply(names(fList), function(ana)
        {
            if (verbose)
                printf("Loading peak intensities for '%s'...\n", ana)
            
            fts <- loadIntensities(cmdQueue[[ana]]$dataFile, fList[[ana]], intSearchRTWindow)
            
            # BUG: OpenMS sporadically reports features with 0 intensity
            # (noticed this with a feature of only two datapoints in hull).
            fts <- fts[intensity > 0]
            
            return(fts)
        }, simplify = FALSE)
    }

    if (verbose)
    {
        printf("Done!\n")
        printFeatStats(fList)
    }

    return(featuresOpenMS(analysisInfo = analysisInfo, features = fList))
}

getOpenMSFFCommand <- function(datafile, out, noiseThrInt, chromSNR, chromFWHM, mzPPM, reEstimateMTSD,
                               traceTermCriterion, traceTermOutliers, minSampleRate,
                               minTraceLength, maxTraceLength, widthFiltering, minFWHM,
                               maxFWHM, traceSNRFiltering, localRTRange, localMZRange,
                               isotopeFilteringModel, MZScoring13C, useSmoothedInts, extraOpts)
{
    boolToChr <- function(b) if (b) "true" else "false"

    settings <- list("-algorithm:common:noise_threshold_int" = noiseThrInt,
                     "-algorithm:common:chrom_peak_snr" = chromSNR,
                     "-algorithm:common:chrom_fwhm" = chromFWHM,
                     "-algorithm:mtd:mass_error_ppm" = mzPPM,
                     "-algorithm:mtd:reestimate_mt_sd" = boolToChr(reEstimateMTSD),
                     "-algorithm:mtd:trace_termination_criterion" = traceTermCriterion,
                     "-algorithm:mtd:trace_termination_outliers" = traceTermOutliers,
                     "-algorithm:mtd:min_sample_rate" = minSampleRate,
                     "-algorithm:mtd:min_trace_length" = minTraceLength,
                     "-algorithm:mtd:max_trace_length" = maxTraceLength,
                     "-algorithm:epd:width_filtering" = widthFiltering,
                     "-algorithm:epd:min_fwhm" = minFWHM,
                     "-algorithm:epd:max_fwhm" = maxFWHM,
                     "-algorithm:ffm:local_rt_range" = localRTRange,
                     "-algorithm:ffm:local_mz_range" = localMZRange,
                     "-algorithm:ffm:isotope_filtering_model" = isotopeFilteringModel,
                     "-algorithm:ffm:mz_scoring_13C" = boolToChr(MZScoring13C),
                     "-algorithm:ffm:use_smoothed_intensities" = boolToChr(useSmoothedInts),
                     "-algorithm:ffm:report_convex_hulls" = "true")

    # figure out if we're running OpenMS version >= 2.5
    oldFFM <- TRUE
    FFMHelp <- suppressWarnings(executeCommand(getCommandWithOptPath("FeatureFinderMetabo", "OpenMS"), stdout = TRUE, stderr = TRUE))
    FFMHelp <- FFMHelp[grepl("Version:", FFMHelp, fixed = TRUE)]
    if (length(FFMHelp) == 1) # should be fine, but fallback to old version just in case...
    {
        FFMVer <- unlist(regmatches(FFMHelp, regexec("[0-9\\.]+", FFMHelp)))
        oldFFM <- utils::compareVersion(FFMVer, "2.5") == -1
    }
    
    if (oldFFM)
        settings <- c(settings, "-algorithm:epd:masstrace_snr_filtering" = boolToChr(traceSNRFiltering))
    else if (traceSNRFiltering) # changed from boolean value to just be there or not...
        settings <- c(settings, "-algorithm:epd:masstrace_snr_filtering" = "")
    
    if (!is.null(extraOpts))
        settings <- modifyList(settings, extraOpts)

    return(list(command = getCommandWithOptPath("FeatureFinderMetabo", "OpenMS"),
                args = c(OpenMSArgListToOpts(settings), "-in", datafile, "-out", out)))
}

importFeatureXML <- function(ffile)
{
    return(as.data.table(parseFeatureXMLFile(ffile)))
}

# OpenMS doesn't support peak intensities. Estimate them from retention times
loadIntensities <- function(dfile, features, rtWindow)
{
    spectra <- loadSpectra(dfile, verbose = FALSE)
    features <- copy(features) # HACK: avoid sR crash caused by data.table

    if (nrow(features) == 0)
        features[, intensity := 0]
    else
        features[, intensity := loadEICIntensities(spectra, features, rtWindow)]

    return(features)
}
