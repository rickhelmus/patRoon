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
#' @param thr Noise intensity threshold. Sets
#'   \code{algorithm:common:noise_threshold_int} option.
#' @param comfwhm,minfwhm,maxfwhm Expected/Minimal/Maximal chromatographic peak
#'   width (seconds). Sets \code{algorithm:common:chrom_fwhm},
#'   \code{algorithm:epd:min_fwhm} and \code{algorithm:epd:max_fwhm},
#'   respectively.
#' @param minlength,maxlength Minimum/Maximum length of mass trace (seconds).
#'   Set negative value for maxlength to disable maximum. Sets
#'   \code{algorithm:mtd:min_trace_length} and
#'   \code{algorithm:mtd:min_trace_length}, respectively.
#' @param mzppm Allowed mass deviation (ppm) for trace detection. Sets
#'   \code{algorithm:mtd:mass_error_ppm}.
#' @param extraOpts Named character \code{vector} containing extra options that
#'   will be passed to \code{FeatureFinderMetabo}. Any options specified here
#'   will override any of the above.
#'
#' @template multiProc-args
#'
#' @rdname feature-finding
#' @export
findFeaturesOpenMS <- function(analysisInfo, thr = 1000, comfwhm = 5, minfwhm = 3, maxfwhm = 60, minlength = 3,
                               maxlength = -1, mzppm = 10, extraOpts = NULL,
                               logPath = file.path("log", "openms"), maxProcAmount = getOption("patRoon.maxProcAmount"))
{
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
            fts <- importFeatureXML(cmd$featFile)
            fts <- loadIntensities(cmd$dataFile, fts)

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
    printf("%s: %d (%.1f%%)\n", analysisInfo$analysis, fCounts, fCounts * 100 / fTotCount)
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

# test with xml2 package: nicer code but slower, so don't use for now
importFeatureXML2 <- function(ffile)
{
    xml <- xml2::read_xml(ffile)

    featPath <- xml_find_all(xml, "/featureMap/featureList/feature")
    featIds <- xml_text(xml_find_all(featPath, "@id"))
    rets <- xml_double(xml_find_all(featPath, "position[@dim = 0]"))
    mzs <- xml_double(xml_find_all(featPath, "position[@dim = 1]"))
    areas <- xml_double(xml_find_all(featPath, "intensity"))

    # this will get ALL retention/mz values for each hulls --> use hull sizes to
    # relate back actual hulls
    hullPath <- xml_find_all(featPath, "convexhull[@nr = 0]")
    hullSizes <- xml_length(hullPath)
    hullTab <- data.table(rt = xml_double(xml_find_all(hullPath, "pt/@x")),
                          mz = xml_double(xml_find_all(hullPath, "pt/@y")),
                          hull = unlist(sapply(seq_along(hullSizes), function(h) rep(h, hullSizes[h]))))

    rtMin <- hullTab[, min(rt), by = "hull"]
    rtMax <- hullTab[, max(rt), by = "hull"]
    mzMin <- hullTab[, min(mz), by = "hull"]
    mzMax <- hullTab[, max(mz), by = "hull"]

    return(data.table(ID = featIds, ret = rets, mz = mzs, area = areas,
                      retmin = rtMin, retmax = rtMax, mzmin = mzMin, mzmax = mzMax))
}

# OpenMS doesn't support peak intensities. Estimate them from retention times
loadIntensities <- function(dfile, features)
{
    spectra <- loadSpectra(dfile, verbose = FALSE)
    features <- copy(features) # HACK: avoid sR crash caused by data.table

    # take max intensity within +/- 5 sec window from retention time
    for (f in seq_len(nrow(features)))
    {
        ft <- features[f]
        eic <- getEIC(spectra, ft$ret + c(-5, 5), c(ft$mzmin, ft$mzmax))
        set(features, f, "intensity", max(eic$intensity))
    }

    return(features)
}
