#' @include features.R
#' @include main.R
NULL

#' @rdname feature-finding
#' @export
featuresOpenMS <- setClass("featuresOpenMS", contains = "features")

#' @details \code{findFeaturesOpenMS} uses the
#'   \href{http://ftp.mi.fu-berlin.de/pub/OpenMS/release2.1.0-documentation/html/TOPP_FeatureFinderMetabo.html}{FeatureFinderMetabo}
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
#' @param comFWHM Expected chromatographic peak width (in seconds). Sets
#'   \code{algorithm:common:chrom_fwhm} option.
#' @param mzPPM Allowed mass deviation (ppm) for trace detection. Sets
#'   \code{algorithm:mtd:mass_error_ppm}.
#' @param reEstimateMTSD If \code{TRUE} then enables dynamic re-estimation of
#'   m/z variance during mass trace collection stage. Sets
#'   \code{algorithm:mtd:reestimate_mt_sd}.
#' @param traceTermCriterion,traceTermOutliers,minSampleRate Termination
#'   criterion for the extension of mass traces. See
#'   \href{http://ftp.mi.fu-berlin.de/pub/OpenMS/release2.1.0-documentation/html/TOPP_FeatureFinderMetabo.html}{FeatureFinderMetabo}.
#'    Sets the \code{algorithm:mtd:trace_termination_criterion},
#'   \code{algorithm:mtd:trace_termination_outliers} and
#'   \code{algorithm:mtd:min_sample_rate} options, respectively.
#' @param minTraceLength,maxTraceLength Minimum/Maximum length of mass trace
#'   (seconds). Set negative value for maxlength to disable maximum. Sets
#'   \code{algorithm:mtd:min_trace_length} and
#'   \code{algorithm:mtd:min_trace_length}, respectively.
#' @param widthFiltering,minFWHM,maxFWHM Enable filtering of unlikely peak
#'   widths. See
#'   \href{http://ftp.mi.fu-berlin.de/pub/OpenMS/release2.1.0-documentation/html/TOPP_FeatureFinderMetabo.html}{FeatureFinderMetabo}.
#'    Sets \code{algorithm:epd:width_filtering}, \code{algorithm:epd:min_fwhm}
#'   and \code{algorithm:epd:max_fwhm}, respectively.
#' @param traceSNRFiltering If \code{TRUE} then apply post-filtering by
#'   signal-to-noise ratio after smoothing. Sets the
#'   \code{algorithm:epd:masstrace_snr_filtering} option.
#' @param localRTRange,localMZRange Retention/MZ range where to look for
#'   coeluting/isotopic mass traces. Sets the
#'   \code{algorithm:ffm:local_rt_range} and \code{algorithm:ffm:local_mz_range}
#'   options, respectively.
#' @param isotopeFilteringModel,MZScoring13C Remove/score candidate assemblies
#'   based on isotope intensities. See
#'   \href{http://ftp.mi.fu-berlin.de/pub/OpenMS/release2.1.0-documentation/html/TOPP_FeatureFinderMetabo.html}{FeatureFinderMetabo}.
#'    Sets the \code{algorithm:ffm:isotope_filtering_model} and
#'   \code{algorithm:ffm:mz_scoring_13C} options, respectivly.
#' @param useSmoothedInts If \code{TRUE} then use LOWESS intensities instead of
#'   raw intensities. Sets the \code{algorithm:ffm:use_smoothed_intensities}
#'   option.
#' @param extraOpts Named character \code{vector} containing extra options that
#'   will be passed to \code{FeatureFinderMetabo}. Any options specified here
#'   will override any of the above.
#' @param intSearchRTWindow Retention time window (in seconds) that is used to
#'   find the closest data point to the retention time to obtain the intensity
#'   of a feature (this is needed since OpenMS does not provide this data).
#' @template multiProc-args
#'
#' @template refs-openms
#'
#' @rdname feature-finding
#' @export
findFeaturesOpenMS <- function(analysisInfo, noiseThrInt = 1000, chromSNR = 3, comFWHM = 5, mzPPM = 10, reEstimateMTSD = TRUE,
                               traceTermCriterion = "sample_rate", traceTermOutliers = 5, minSampleRate = 0.5,
                               minTraceLength = 3, maxTraceLength = -1, widthFiltering = "fixed", minFWHM = 3,
                               maxFWHM = 60, traceSNRFiltering = FALSE, localRTRange = 10, localMZRange = 6.5,
                               isotopeFilteringModel = "metabolites (5% RMS)", MZScoring13C = FALSE, useSmoothedInts = FALSE,
                               extraOpts = NULL, intSearchRTWindow = 3,
                               logPath = file.path("log", "openms"), maxProcAmount = getOption("patRoon.maxProcAmount"))
{
    ac <- checkmate::makeAssertCollection()
    assertAnalysisInfo(analysisInfo, "mzML", add = ac)
    aapply(checkmate::assertNumber, . ~ thr + minfwhm + maxfwhm + minlength + mzppm, lower = 0,
           finite = TRUE, fixed = list(add = ac))
    checkmate::assertNumber(maxlength, finite = TRUE, add = ac)
    checkmate::assertList(extraOpts, any.missing = FALSE, names = "unique", null.ok = TRUE, add = ac)
    assertMultiProcArgs(logPath, maxProcAmount, add = ac)
    checkmate::reportAssertions(ac)
    
    printf("Finding features with OpenMS for %d analyses ...\n", nrow(analysisInfo))

    cmdQueue <- lapply(seq_len(nrow(analysisInfo)), function(anai)
    {
        dfile <- getMzMLAnalysisPath(analysisInfo$analysis[anai], analysisInfo$path[anai])
        ffile <- tempfile(analysisInfo$analysis[anai], fileext = ".featureXML")
        hash <- makeHash(makeFileHash(dfile), thr, comfwhm, minfwhm, maxfwhm, minlength,
                         maxlength, mzppm, extraOpts)
        cmd <- getOpenMSFFCommand(dfile, ffile, thr, comfwhm, minfwhm, maxfwhm, minlength, maxlength, mzppm, extraOpts)

        logf <- if (!is.null(logPath)) file.path(logPath, paste0("ffm-", analysisInfo$analysis[anai], ".txt")) else NULL

        return(c(list(hash = hash, dataFile = dfile, featFile = ffile, stdoutFile = logf), cmd))
    })
    names(cmdQueue) <- analysisInfo$analysis

    cachedResults <- sapply(cmdQueue, function(cmd) loadCacheData("featuresOpenMS", cmd$hash), simplify = FALSE)
    cachedResults <- cachedResults[!sapply(cachedResults, is.null)]
    cmdQueue <- cmdQueue[setdiff(names(cmdQueue), names(cachedResults))] # remove cached results

    fList <- list()
    if (length(cmdQueue) > 0)
    {
        if (!is.null(logPath))
            mkdirp(logPath)

        fList <- executeMultiProcess(cmdQueue, function(cmd, ...)
        {
            # fts <- importFeatureXML(cmd$featFile)
            fts <- importFeatureXMLCpp(cmd$featFile)
            fts <- loadIntensitiesCPP(cmd$dataFile, fts, intSearchRTWindow)

            # BUG: OpenMS sporadically reports features with 0 intensity
            # (noticed this with a feature of only two datapoints in hull).
            fts <- fts[intensity > 0]

            saveCacheData("featuresOpenMS", fts, cmd$hash)

            return(fts)
        }, maxProcAmount = maxProcAmount)
    }

    if (length(cachedResults) > 0)
    {
        fList <- c(fList, cachedResults)
        fList <- fList[intersect(analysisInfo$analysis, names(fList))] # re-order
    }

    fCounts <- sapply(fList, nrow)
    fTotCount <- sum(fCounts)
    printf("Done! Feature statistics:\n")
    printf("%s: %d (%.1f%%)\n", analysisInfo$analysis, fCounts, if (fTotCount == 0) 0 else fCounts * 100 / fTotCount)
    printf("Total: %d\n", fTotCount)

    return(featuresOpenMS(analysisInfo = analysisInfo, features = fList))
}

getOpenMSFFCommand <- function(datafile, out, thr, comfwhm, minfwhm, maxfwhm, minlength, maxlength, mzppm, extraOpts)
{
    settings <- c("-algorithm:common:noise_threshold_int" = thr,
                  "-algorithm:common:chrom_fwhm" = comfwhm,
                  "-algorithm:mtd:mass_error_ppm" = mzppm,
                  "-algorithm:mtd:trace_termination_criterion" = "sample_rate",
                  "-algorithm:mtd:min_trace_length" = minlength,
                  "-algorithm:mtd:max_trace_length" = maxlength, # positive value makes it much slower
                  "-algorithm:epd:width_filtering" = "fixed",
                  "-algorithm:epd:min_fwhm" = minfwhm,
                  "-algorithm:epd:max_fwhm" = maxfwhm,
                  "-algorithm:ffm:report_convex_hulls" = "true")

    if (!is.null(extraOpts))
        settings <- c(extraOpts, settings)

    # convert to unnamed character vector where previous names are followed by set values
    settings <- as.vector(sapply(names(settings), function(s) c(s, settings[[s]])))

    return(list(command = getCommandWithOptPath("FeatureFinderMetabo", "OpenMS"),
                args = c(settings, "-in", datafile, "-out", out)))
}

importFeatureXML <- function(ffile)
{
    doc <- XML::xmlTreeParse(ffile)
    docrt <- XML::xmlRoot(doc)

    fcount <- as.numeric(XML::xmlAttrs(docrt[["featureList"]])[["count"]])
    ret <- data.table(ID=character(fcount), ret=numeric(fcount), mz=numeric(fcount), area=numeric(fcount),
                      retmin=numeric(fcount), retmax=numeric(fcount), mzmin=numeric(fcount), mzmax=numeric(fcount))
    ind <- 1

    XML::xmlSApply(docrt[["featureList"]], function(xn)
    {
        ret[[ind, "ID"]] <<- XML::xmlAttrs(xn)[["id"]]
        ret[[ind, "area"]] <<- as.numeric(XML::xmlValue(xn[["intensity"]]))

        lapply(XML::xmlElementsByTagName(xn, "position"), function(xsn)
        {
            if (XML::xmlAttrs(xsn)[["dim"]] == "0")
                ret[[ind, "ret"]] <<- as.numeric(XML::xmlValue(xsn))
            else
                ret[[ind, "mz"]] <<- as.numeric(XML::xmlValue(xsn))
        })

        # assume first hull describes main ion
        hull <- XML::xmlSApply(xn[["convexhull"]], XML::xmlAttrs)
        # odd indices (1, 3, ...): retention times, others (2, 4, ...): mz's
        ret[[ind, "retmin"]] <<- min(as.numeric(hull[c(T, F)]))
        ret[[ind, "retmax"]] <<- max(as.numeric(hull[c(T, F)]))
        ret[[ind, "mzmin"]] <<- min(as.numeric(hull[c(F, T)]))
        ret[[ind, "mzmax"]] <<- max(as.numeric(hull[c(F, T)]))

        ind <<- ind + 1
    })

    return(ret)
}

importFeatureXMLCpp <- function(ffile)
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
    {
        # find closest data point to retention time in +/- 5 sec window
        for (f in seq_len(nrow(features)))
        {
            ft <- features[f]
            eic <- getEIC(spectra, ft$ret + c(-rtWindow, rtWindow), c(ft$mzmin, ft$mzmax))
            set(features, f, "intensity", eic$intensity[which.min(abs(ft$ret - eic$time))])
        }
    }
    
    return(features)
}

loadIntensitiesCPP <- function(dfile, features, rtWindow)
{
    spectra <- loadSpectra(dfile, verbose = FALSE)
    features <- copy(features) # HACK: avoid sR crash caused by data.table
    
    if (nrow(features) == 0)
        features[, intensity := 0]
    else
        features[, intensity := loadEICIntensities(spectra, features, rtWindow)]

    return(features)
}
