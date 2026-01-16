# SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' MS data conversion
#'
#' Conversion of MS analysis files between several open and closed data formats.
#'
#' @param algorithm Either \code{"pwiz"} (ProteoWizard), \code{"openms"}, \code{"bruker"} (Bruker DataAnalysis) ,
#'   \code{"im_collapse"} or \code{"timsconvert"}.
#' @param direction A \code{character} specifying the direction of conversion. Either \code{"input"} or \code{"output"}.
#' @param inFiles,outFiles A \code{character} vector with input and output files, respectively. Lengths and order should
#'   be the same.
#' @param type,typeFrom,typeTo The type of the input or output files. See \code{getMSConversionTypes} for the supported
#'   types.
#' @param formatFrom,formatTo The input or output format. See \code{getMSConversionFormats} for the supported formats.
#' @param centroid Set to \code{TRUE} to perform centroiding.
#'
#'   For \code{convertMSFilesPWiz}: the value may be \code{"vendor"} to perform centroiding with the vendor algorithm or
#'   \code{"cwt"} to use ProteoWizard's wavelet algorithm.
#' @param IMS How to handle IMS data.
#'
#'   For \code{convertMSFilesPWiz}: if \code{TRUE} then IMS data is exported and spectra for each IMS frame are combined
#'   into a single spectrum (using the \command{--combineIonMobilitySpectra} option), which is the format supported by
#'   \pkg{patRoon}. Set to \code{NA} to collapse the IMS data by scan summing, which mimics 'regular' HRMS data. Set to
#'   \code{FALSE} for non-IMS data. \strong{NOTE}: do not set \code{IMS=FALSE} if the data has IMS data. This will
#'   result in very large files where MS spectra are not combined by frame, which \strong{cannot} be properly read by
#'   \pkg{patRoon}.
#'
#'   For \code{convertMSFilesTIMSCONVERT}: set to \code{TRUE} to keep IMS data or \code{FALSE} to exclude IMS data to
#'   mimic 'regular' LC-MS data.
#' @param extraOpts A \code{character} vector specifying any extra command line parameters passed to \command{msconvert}
#'   or \command{FileConverter}. Set to \code{NULL} to ignore. For options: see
#'   \href{https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_FileConverter.html}{FileConverter}
#'   and \href{http://proteowizard.sourceforge.net/tools/msconvert.html}{msconvert}.
#' @param overwrite Should existing destination file be overwritten (\code{TRUE}) or not (\code{FALSE})?
#' @param \dots For \code{convertMSFilesIMSCollapse}: further arguments passed to
#'   \code{\link[mzR:writeMSData]{mzR::writeMSData}}.
#'
#'   For \code{convertMSFilePaths} and \code{convertMSFiles}: further arguments passed to algorithm specific conversion
#'   functions.
#'
#' @templateVar what \code{convertMSFilesPWiz}, \code{convertMSFilesOpenMS} and \code{convertMSFilesTIMSCONVERT}
#' @template uses-multiProc
#'
#' @references \insertRef{Rst2016}{patRoon} \cr\cr
#'
#'   \insertRef{Chambers2012}{patRoon} \cr\cr
#'
#'   \insertRef{Luu2022}{patRoon} \cr\cr
#'
#'   \addCitations{mzR}
#'
#' @name MSConversion
NULL

#' @details \code{getMSConversionTypes} returns a \code{character} with all supported input or output conversion types
#'   for an algorithm.
#' @rdname MSConversion
#' @export
getMSConversionTypes <- function(algorithm, direction)
{
    assertMSConvAlgo(algorithm)
    checkmate::assertChoice(direction, c("input", "output"))
 
    if (direction == "input")
    {
        return(switch(algorithm,
                      pwiz = getMSFileTypes(),
                      openms = c("centroid", "profile"),
                      bruker = "raw",
                      im_collapse = c("raw", "ims"),
                      timsconvert = "raw"))
    }

    return(switch(algorithm,
                  pwiz = c("centroid", "profile", "ims"),
                  openms = c("centroid"),
                  bruker = c("centroid", "profile"),
                  im_collapse = "centroid",
                  timsconvert = c("centroid", "profile", "ims")))
}

#' @details \code{getMSConversionFormats} returns a \code{character} with all supported input or output conversion
#'   formats for an algorithm, optionally filtered by the given \code{type}.
#' @rdname MSConversion
#' @export
getMSConversionFormats <- function(algorithm, direction, type = NULL)
{
    assertMSConvAlgo(algorithm)
    checkmate::assertChoice(direction, c("input", "output"))
    checkmate::assertChoice(type, getMSConversionTypes(algorithm, direction), null.ok = TRUE)

    ret <- if (direction == "input")
    {
        if (algorithm == "pwiz")
            getMSFileFormats()
        else if (algorithm == "openms")
            c("mzML", "mzXML")
        else if (algorithm == "bruker")
            "bruker"
        else if (algorithm == "im_collapse")
            c("bruker_ims", "mzML")
        else if (algorithm == "timsconvert")
            "bruker_ims"
    }
    else # output
    {
        if (algorithm == "timsconvert")
            "mzML"
        else
            c("mzML", "mzXML")
    }

    if (!is.null(type))
        ret <- intersect(ret, getMSFileFormats(type))
    
    return(ret)
}

#' @details \code{convertMSFilesPWiz} converts and pre-treats HRMS data with the \command{msconvert} tool from
#'   \href{http://proteowizard.sourceforge.net/}{ProteoWizard}.
#'
#' @param minIntensity The minimum intensity of the mass peaks to be kept. Applying an intensity threshold is especially
#'   beneficial to reduce export file size when there are a lot of zero or very low intensity mass peaks. \strong{NOTE}
#'   this currently does \emph{not} work well with IMS data.
#' @param filters A \code{character} vector specifying one or more filters to \command{msconvert}. The elements of the
#'   specified vector are directly passed to the \code{--filter} option (see
#'   \href{http://proteowizard.sourceforge.net/tools/filters.html}{here})
#' @param PWizBatchSize The number of analyses to process by a single call to \command{msconvert}. Usually a value of
#'   one is most efficient. Set to zero to run all analyses all at once from a single call.
#'
#' @rdname MSConversion
#' @export
convertMSFilesPWiz <- function(inFiles, outFiles, formatTo = "mzML", centroid = TRUE, IMS = FALSE, minIntensity = 0,
                               filters = NULL, extraOpts = NULL, PWizBatchSize = 1)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertCharacter(inFiles, min.chars = 1, min.len = 1, add = ac)
    checkmate::assertCharacter(outFiles, min.chars = 1, len = length(inFiles), add = ac)
    checkmate::assertChoice(formatTo, getMSConversionFormats("pwiz", "output"), add = ac)
    checkmate::assert(checkmate::checkFlag(centroid),
                      checkmate::checkChoice(centroid, c("vendor", "cwt")),
                      .var.name = "centroid", add = ac)
    checkmate::assertFlag(IMS, na.ok = TRUE, add = ac)
    checkmate::assertNumber(minIntensity, lower = 0, finite = TRUE, add = ac)
    checkmate::assertCharacter(filters, min.chars = 1, null.ok = TRUE, add = ac)
    checkmate::assertCharacter(extraOpts, min.chars = 1, null.ok = TRUE, add = ac)
    checkmate::assertCount(PWizBatchSize, positive = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (!isFALSE(IMS) && minIntensity > 0)
        warning("The minIntensity argument currently results in incorrect file conversion of IMS data!", call. = FALSE)
    
    if (centroid != FALSE)
        filters <- c(paste("peakPicking", if (is.character(centroid)) centroid else ""), filters)
    
    if (is.na(IMS))
        filters <- c(filters, "scanSumming sumMs1=1")
    
    if (minIntensity > 0)
        filters <- c(filters, sprintf("threshold absolute %f most-intense", minIntensity))

    mainArgs <- paste0("--", formatTo)
    if (!isFALSE(IMS))
        mainArgs <- c(mainArgs, "--combineIonMobilitySpectra")
    if (!is.null(filters))
        mainArgs <- c(mainArgs, sapply(filters, function(f) c("--filter", f)))
    if (!is.null(extraOpts))
        mainArgs <- c(mainArgs, extraOpts)

    pwpath <- findPWizPath()
    if (is.null(pwpath) || !file.exists(file.path(pwpath, paste0("msconvert", if (Sys.info()[["sysname"]] == "Windows") ".exe" else ""))))
        stop("Could not find ProteoWizard. You may set its location in the patRoon.path.pwiz option. See ?patRoon for more details.")
    msc <- file.path(pwpath, "msconvert")

    if (PWizBatchSize != 1 && length(inFiles) > 1)
    {
        outDir <- dirname(outFiles)
        if (!allSame(outDir)) # UNDONE?
            stop("If PWizBatchSize>1 then all output files must go to the same directory.")
        outDir <- outDir[1]
        
        if (PWizBatchSize == 0)
            batches <- list(seq_along(inFiles))
        else
            batches <- splitInBatches(seq_along(inFiles), PWizBatchSize)
    
        cmdQueue <- lapply(seq_along(batches), function(bi)
        {
            input <- tempfile("msconvert")
            cat(inFiles[batches[[bi]]], sep = "\n", file = input)
            logf <- paste0("pwiz-batch_", bi, ".txt")
            # UNDONE: unlike PWizBatchSize==1 we don't (can't) set output file names here, is this a problem?
            return(list(logFile = logf, command = msc, args = c("-f", input,  "-o", outDir, mainArgs)))
        })
    }
    else
    {
        cmdQueue <- lapply(seq_along(inFiles), function(fi)
        {
            basef <- basename(tools::file_path_sans_ext(inFiles[fi]))
            logf <- paste0("pwiz-", basef, ".txt")
            return(list(logFile = logf, command = msc,
                        args = c(inFiles[fi], "--outfile", outFiles[fi],
                                 "-o", dirname(outFiles[fi]), mainArgs)))
        })
    }

    executeMultiProcess(cmdQueue, function(cmd) {}, logSubDir = "convert")

    invisible(NULL)
}

#' @details \code{convertMSFilesOpenMS} converts HRMS data with the \command{FileConvert} tool of
#'   \href{http://www.openms.de/}{OpenMS}.
#' @rdname MSConversion
#' @export
convertMSFilesOpenMS <- function(inFiles, outFiles, formatTo = "mzML", extraOpts = NULL)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertCharacter(inFiles, min.chars = 1, min.len = 1, add = ac)
    checkmate::assertCharacter(outFiles, min.chars = 1, len = length(inFiles), add = ac)
    checkmate::assertChoice(formatTo, getMSConversionFormats("openms", "output"), add = ac)
    checkmate::assertCharacter(extraOpts, min.chars = 1, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    mainArgs <- c()
    if (!is.null(extraOpts))
        mainArgs <- c(mainArgs, extraOpts)

    msc <- getExtDepPath("openms", "FileConverter")
    cmdQueue <- lapply(seq_along(inFiles), function(fi)
    {
        basef <- basename(tools::file_path_sans_ext(inFiles[fi]))
        logf <- paste0("openms-", basef, ".txt")
        return(list(logFile = logf, command = msc,
                    args = c("-in", inFiles[fi], "-out", outFiles[fi], mainArgs)))
    })

    executeMultiProcess(cmdQueue, function(cmd) {}, logSubDir = "convert")

    invisible(NULL)
}

#' @details \code{convertMSFilesBruker} converts and pre-treats Bruker HRMS data with Bruker DataAnalysis. Note that
#'   TIMS data currently is not supported.
#' @rdname MSConversion
#' @export
convertMSFilesBruker <- function(inFiles, outFiles, formatTo = "mzML", centroid = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertCharacter(inFiles, min.chars = 1, min.len = 1, add = ac)
    checkmate::assertCharacter(outFiles, min.chars = 1, len = length(inFiles), add = ac)
    checkmate::assertChoice(formatTo, getMSConversionFormats("bruker", "output"), add = ac)
    checkmate::assertFlag(centroid, add = ac)
    checkmate::reportAssertions(ac)
    
    # expConstant <- if (to == "mzXML") DAConstants$daMzXML else if (to == "mzData") DAConstants$daMzData else DAConstants$daMzML
    expConstant <- if (formatTo == "mzXML") DAConstants$daMzXML else DAConstants$daMzML
    expSpecConstant <- if (centroid) DAConstants$daLine else DAConstants$daProfile

    DA <- getDAApplication()
    hideDAInScope()

    fCount <- length(inFiles)
    prog <- openProgBar(0, fCount)

    for (i in seq_len(fCount))
    {
        ind <- getDAFileIndex(DA, inFiles[i], NULL)
        if (ind == -1)
            warning(paste("Failed to open file in DataAnalysis:", inFiles[i]))
        else
            DA[["Analyses"]][[ind]]$Export(outFiles[i], expConstant, expSpecConstant)

        setTxtProgressBar(prog, i)
    }

    setTxtProgressBar(prog, fCount)

    invisible(NULL)
}

#' @details \code{convertMSFilesIMSCollapse} is used to convert IMS data to data that mimics 'regular' HRMS data by
#'   collapsing the IMS dimension. The raw data interface of \pkg{patRoon} first sums up all spectra within each IMS
#'   frame, performs centroiding and finally exports the resulting data with the
#'   \code{\link[mzR:writeMSData]{mzR::writeMSData}} function. Several thresholds can be set to speed up the conversion
#'   process and reduce noise, but care should be taken that no mass peaks of interest are lost.
#'
#' @param mzRange,mobilityRange A two sized vector specifying the m/z and mobility range to be exported, respectively.
#'   Set to \code{NULL} to export the full range.
#' @param smoothWindow,halfWindow,maxGap Centroiding parameters: see \code{\link{getDefAvgPListParams}} for details.
#' @param clusterMethod,mzWindow The clustering method and window (see \link[=cluster-params]{clustering parameters})
#'   used to find and combine MS/MS spectra of precursors with close \emph{m/z}.
#' @param minIntensityIMS The minimum intensity for MS peaks in raw data.
#' @param includeMSMS Set to \code{TRUE} to include MS/MS spectra in the output. For IMS workflows where IMS data is
#'   only collapsed to produce compatible data files for feature detection, MS/MS data are not needed and can be
#'   excluded to reduce computational times and file sizes. Setting \code{includeMSMS=TRUE} is primarily intended to
#'   perform 'classical LC-MS workflows' with IMS data.
#'
#' @templateVar what \code{convertMSFilesIMSCollapse}
#' @template uses-msdata
#'
#' @rdname MSConversion
#' @export
convertMSFilesIMSCollapse <- function(inFiles, outFiles, typeFrom, formatTo = "mzML", mzRange = NULL, mobilityRange = NULL,
                                      smoothWindow = 0, halfWindow = 2, maxGap = defaultLim("mz", "medium"),
                                      clusterMethod = "distance_mean", mzWindow = defaultLim("mz", "medium"),
                                      minIntensityIMS = 0, includeMSMS = FALSE, ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertCharacter(inFiles, min.chars = 1, min.len = 1, add = ac)
    checkmate::assertCharacter(outFiles, min.chars = 1, len = length(inFiles), add = ac)
    checkmate::assertChoice(typeFrom, c("raw", "ims"), add = ac)
    checkmate::assertChoice(formatTo, getMSConversionFormats("im_collapse", "output"), add = ac)
    aapply(checkmate::assertCount, . ~ smoothWindow + halfWindow, positive = c(FALSE, TRUE), fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ maxGap + mzWindow + minIntensityIMS, lower = 0, finite = TRUE, fixed = list(add = ac))
    aapply(assertRange, . ~ mzRange + mobilityRange, null.ok = TRUE, fixed = list(add = ac))
    assertClusterMethod(clusterMethod, add = ac)
    checkmate::assertFlag(includeMSMS, add = ac)
    checkmate::reportAssertions(ac)
        
    makeHeader <- function(collapsedSpectra, MSLevel, meta)
    {
        # UNDONE: assume there is no polarity switching (currently not possible with IMS instruments anyway)
        # UNDONE: support polarities for MSTK
        pol <- if (meta$polarity[1] == 1) 1 else 0
        times <- if (MSLevel == 1) meta$time else meta$time[match(collapsedSpectra$framesMS2, meta$scan)]
        spectra <- if (MSLevel == 1) collapsedSpectra$MS1 else collapsedSpectra$MS2
        
        # NOTE: there may be empty spectra due to filtering or measurement errors
        nonEmptySpectra <- sapply(spectra, nrow) > 0
        specMZRange <- range(unlist(lapply(spectra[nonEmptySpectra], function(sp) range(sp[, "mz"]))))
        
        # NOTE: scan nums are filled in later
        data.frame(
            msLevel = MSLevel,
            polarity = pol,
            peaksCount = sapply(spectra, nrow), # UNDONE: correct?
            totIonCurrent = sapply(spectra, function(sp) if (nrow(sp) == 0) 0 else sum(sp[, "intensity"])),
            retentionTime = times,
            basePeakMZ = sapply(spectra, function(sp) if (nrow(sp) == 0) 0 else sp[which.max(sp[, "intensity"]), "mz"]),
            basePeakIntensity = sapply(spectra, function(sp) if (nrow(sp) == 0) 0 else max(sp[, "intensity"])),
            collisionEnergy = NA_real_,
            ionisationEnergy = 0,
            lowMZ = 0,
            highMZ = 0,
            precursorScanNum = NA_integer_,
            precursorMZ = if (MSLevel == 1) NA_real_ else collapsedSpectra$precursorMZs,
            precursorCharge = NA_integer_,
            precursorIntensity = NA_real_,
            mergedScan = NA_integer_,
            mergedResultScanNum = NA_integer_,
            mergedResultStartScanNum = NA_integer_,
            mergedResultEndScanNum = NA_integer_,
            injectionTime = 0,
            filterString = NA_character_,
            centroided = TRUE,
            ionMobilityDriftTime = NA_real_,
            isolationWindowTargetMZ = if (MSLevel == 1) NA_real_ else collapsedSpectra$precursorMZs,
            isolationWindowLowerOffset = if (MSLevel == 1) NA_real_ else collapsedSpectra$isolationStarts,
            isolationWindowUpperOffset = if (MSLevel == 1) NA_real_ else collapsedSpectra$isolationEnds,
            # UNDONE: below are technically not scan limits, but might be good enough?
            scanWindowLowerLimit = specMZRange[1],
            scanWindowUpperLimit = specMZRange[2],
            frame = if (MSLevel == 1) meta$scan else collapsedSpectra$framesMS2 # temporary, will be used and removed later
        )
    }
    
    # HACK: convert inFiles + outFiles to analysis information so we can use applyMSData
    inFilesAna <- basename(tools::file_path_sans_ext(inFiles))
    outFilesAna <- basename(tools::file_path_sans_ext(inFiles))
    
    if (!all.equal(inFilesAna, outFilesAna))
        stop("Input and output files must have the same base name", call. = FALSE)
    
    anaInfo <- data.table(analysis = inFilesAna, path_centroid = dirname(outFiles))
    anaInfo[, (paste0("path_", typeFrom)) := dirname(inFiles)]
    mkdirp(anaInfo$path_centroid)
    
    printf("Collapsing all %d analyses ...\n", nrow(anaInfo))
    applyMSData(anaInfo, anaInfo$path_centroid, needTypes = "ims", showProgress = TRUE, func = function(ana, path, backend, outd)
    {
        outp <- file.path(outd, paste0(ana, ".", formatTo))
        
        openMSReadBackend(backend, path)
        collapsedSpectra <- collapseIMSFrames(backend, NULLToZero(mzRange[1]), NULLToZero(mzRange[2]),
                                              NULLToZero(mobilityRange[1]), NULLToZero(mobilityRange[2]), smoothWindow,
                                              halfWindow, maxGap, clusterMethod, mzWindow, minIntensityIMS, includeMSMS)
        
        
        
        headerMS1 <- makeHeader(collapsedSpectra, 1, getMSMetadata(backend, 1))
        
        header <- if (!includeMSMS)
            headerMS1
        else
        {
            meta <- getMSMetadata(backend, 2)
            if (nrow(meta) > 0)
            {
                headerMS2 <- makeHeader(collapsedSpectra, 2, meta)
                rbind(headerMS1, headerMS2)
            }
            else
                headerMS1
        }
        
        ord <- order(header$frame, header$precursorMZ)
        header <- header[ord, ]
        header$seqNum <- header$acquisitionNum <- seq_len(nrow(header))
        header$spectrumId <- sprintf("merged=%d frame=%d", header$seqNum, header$frame)
        header <- header[, names(header) != "frame"]

        allSpectra <- c(collapsedSpectra$MS1, collapsedSpectra$MS2)
        allSpectra <- allSpectra[ord]

        mzR::writeMSData(allSpectra, outp, header, outformat = formatTo, ...)
    })
    
    invisible(NULL)
}

#' @details \code{convertMSFilesTIMSCONVERT} converts and pre-treats TIMS data with
#'   \href{https://gtluu.github.io/timsconvert/}{TIMSCONVERT}. The \code{\link{installTIMSCONVERT}} function can be used
#'   to automatically install \command{TIMSCONVERT}.
#'
#' @param centroidRaw Only applicable if \code{IMS=FALSE}. Sets the \code{mode} parameter of \command{TIMSCONVERT}:
#'   \code{raw} if \code{centroidRaw=TRUE} or \code{centroid} if \code{centroidRaw=FALSE}. See
#'   \url{https://gtluu.github.io/timsconvert/local.html#notes-on-mode-parameter} for more details.
#' @param virtualenv The virtual Python environment in which \command{TIMSCONVERT} is installed. This is passed to
#'   \code{\link[reticulate:use_virtualenv]{reticulate::use_virtualenv}}, which will ensure that the
#'   \command{TIMSCONVERT} command line utility can be found by \pkg{patRoon}. Set to \code{NULL} to skip this step.
#'
#' @rdname MSConversion
#' @export
convertMSFilesTIMSCONVERT <- function(inFiles, outFiles, formatTo = "mzML", centroid = TRUE, centroidRaw = FALSE,
                                      IMS = FALSE, extraOpts = NULL, virtualenv = "patRoon-TIMSCONVERT")
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertCharacter(inFiles, min.chars = 1, min.len = 1, add = ac)
    checkmate::assertCharacter(outFiles, min.chars = 1, len = length(inFiles), add = ac)
    checkmate::assertChoice(formatTo, getMSConversionFormats("timsconvert", "output"), add = ac)
    aapply(checkmate::assertFlag, . ~ centroid + centroidRaw + IMS, fixed = list(add = ac))
    checkmate::assertCharacter(extraOpts, min.chars = 1, null.ok = TRUE, add = ac)
    checkmate::assertCharacter(virtualenv, min.chars = 1, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    # NOTE: activate virtualenv will setup the right PATH
    if (!is.null(virtualenv))
    {
        reticulate::use_virtualenv(virtualenv)
        # UNDONE: the path is only added when loading the virtualenv on Windows?
        if (Sys.info()[["sysname"]] != "Windows")
            withr::local_path(file.path(reticulate::py_config()$virtualenv, "bin"))
    }
    
    mainArgs <- character()
    if (IMS)
    {
        # UNDONE: docs say to use "raw" and not "centroid", but in reality only "centroid" seems to work
        mainArgs <- c(mainArgs, "--mode", "centroid")
    }
    else if (centroid && centroidRaw)
        mainArgs <- c(mainArgs, "--mode", "raw")
    else if (centroid)
        mainArgs <- c(mainArgs, "--mode", "centroid")
    else
        mainArgs <- c(mainArgs, "--mode", "profile")
    
    if (!IMS)
        mainArgs <- c(mainArgs, "--exclude_mobility")
    
    if (!is.null(extraOpts))
        mainArgs <- c(mainArgs, extraOpts)
    
    inFiles <- normalizePath(inFiles); outFiles <- normalizePath(outFiles, mustWork = FALSE)
    cmdQueue <- Map(inFiles, outFiles, f = function(inF, outF)
    {
        if (basename(tools::file_path_sans_ext(inF)) != basename(tools::file_path_sans_ext(outF)))
            stop(sprintf("Input and output files must have the same base name: in '%s' - out '%s'", inF, outF),
                 call. = FALSE)
        logf <- paste0("timsconvert-", basename(tools::file_path_sans_ext(inF)), ".txt")
        return(list(logFile = logf, command = "timsconvert", args = c("--input", inF, "--outdir", dirname(outF), mainArgs)))
    })
    
    executeMultiProcess(cmdQueue, function(cmd) {}, logSubDir = "convert")
    
    invisible(NULL)
}

#' @details \code{convertMSFilePaths} is a wrapper function that simplifies the use of algorithm specific MS conversion
#'   functions, such as \code{convertMSFilesPWiz}, and \code{convertMSFilesTIMSCONVERT}.
#'
#' @param files,dirs The \code{files} argument should be a \code{character} vector with input files. If \code{files}
#'   contains directories and \code{dirs=TRUE} then files from these directories are also considered.
#' @param outPath A character vector specifying directories that should be used for the output. Will be re-cycled if
#'   necessary. If \code{NULL}, output directories will be kept the same as the input directories.
#'
#' @rdname MSConversion
#' @export
convertMSFilePaths <- function(files, formatFrom, formatTo = "mzML", outPath = NULL, dirs = TRUE, overwrite = FALSE,
                               algorithm = "pwiz", ...)
{
    ac <- checkmate::makeAssertCollection()
    assertConvertMSFilesArgs(formatFrom, formatTo, overwrite, algorithm, add = ac)
    checkmate::assertCharacter(files, min.len = 1, min.chars = 1, add = ac)
    checkmate::assertCharacter(outPath, min.chars = 1, min.len = 1, null.ok = TRUE, add = ac)
    assertCanCreateDirs(outPath, add = ac)
    checkmate::assertFlag(dirs, add = ac)
    checkmate::reportAssertions(ac)
    
    if (dirs)
    {
        if (formatFrom == formatTo)
            warning("Input and output formats are the same", call. = FALSE)
        
        dirs <- files[file.info(files, extra_cols = FALSE)$isdir]
        dirs <- dirs[sapply(dirs, function(d)
        {
            # NOTE: with some formats the analysis files are directories --> remove these
            !any(sapply(formatFrom, verifyFileForFormat, path = d))
        })]
        
        dirFiles <- listMSFiles(dirs, formatFrom)
        files <- union(dirFiles, setdiff(files, dirs))
    }
    
    if (is.null(outPath))
        outPath <- dirname(files)
    
    mkdirp(outPath)
    
    # NOTE: use normalizePath() here to convert to backslashes on Windows: needed by msconvert
    outPath <- normalizePath(rep(outPath, length.out = length(files)), mustWork = TRUE)
    files <- normalizePath(files, mustWork = FALSE) # no mustWork, file existence will be checked later
    
    basef <- basename(tools::file_path_sans_ext(files))
    output <- normalizePath(file.path(outPath, paste0(basef, ".", formatTo)), mustWork = FALSE)
    
    keepFiles <- sapply(seq_along(files), function(fi)
    {
        if (!file.exists(files[fi]))
            printf("Skipping non-existing input analysis %s\n", files[fi])
        else if (!overwrite && file.exists(output[fi]))
            printf("Skipping existing output analysis %s\n", output[fi])
        else
            return(TRUE)
        return(FALSE)
    }, USE.NAMES = FALSE)
    
    if (is.logical(keepFiles) && any(keepFiles))
    {
        files <- files[keepFiles]
        output <- output[keepFiles]
        
        if (algorithm == "pwiz")
            convertMSFilesPWiz(files, output, formatTo, ...)
        else if (algorithm == "openms")
            convertMSFilesOpenMS(files, output, formatTo, ...)
        else if (algorithm == "bruker")
            convertMSFilesBruker(files, output, formatTo, ...)
        else if (algorithm == "im_collapse")
            convertMSFilesIMSCollapse(files, output, formatTo = formatTo, ...)
        else # algorithm == "timsconvert"
            convertMSFilesTIMSCONVERT(files, output, formatTo = formatTo, ...)
    }
}

#' @details \code{convertMSFiles} is a wrapper function that simplifies the use of \code{convertMSFilePaths}.
#'
#' @param anaInfo An \link[=analysis-information]{analysis info table} that is used to retrieve the input files. The
#'   paths set by \code{path_centroid}, \code{path_profile} and \code{path_ims} are used to determine the output
#'   directories. This function automatically determines if and how centroiding and IMS conversions should be applied.
#'
#' @param centroidVendor Only for \code{algorithm="pwiz"}: whether centroiding should be performed with vendor
#'   algorithms.
#'
#' @rdname MSConversion
#' @export
convertMSFiles <- function(anaInfo, typeFrom = "raw", typeTo = "centroid", formatFrom, formatTo = "mzML",
                           overwrite = FALSE, algorithm = "pwiz", centroidVendor = TRUE, ...)
{
    ac <- checkmate::makeAssertCollection()
    assertConvertMSFilesArgs(formatFrom, formatTo, overwrite, algorithm, add = ac)
    checkmate::assertChoice(typeFrom, getMSConversionTypes(algorithm, "input"), add = ac)
    checkmate::assertChoice(typeTo, getMSConversionTypes(algorithm, "output"), add = ac)
    checkmate::assertFlag(centroidVendor, add = ac)
    checkmate::reportAssertions(ac)

    anaInfo <- assertAndPrepareAnaInfo(anaInfo, typeFrom, formatFrom)
    
    outPath <- getPathsFromAnaInfo(anaInfo, typeTo)
    if (is.null(outPath) || anyNA(outPath) || any(!nzchar(outPath)))
        stop("Please properly configure the analysis info for the output type ", typeTo, call. = FALSE)
    
    files <- getMSFilesFromAnaInfo(anaInfo, typeFrom, formatFrom)

    args <- list(files = files, formatFrom = formatFrom, formatTo = formatTo, outPath = outPath, dirs = FALSE,
                 overwrite = overwrite, algorithm = algorithm, ...)

    if (algorithm == "pwiz")
    {
        centroid <- if (typeTo == "profile")
            FALSE
        else if (centroidVendor)
            "vendor"
        else
            TRUE
        
        IMS <- formatFrom %in% c("agilent_ims", "bruker_ims") || typeFrom == "ims"
        if (IMS)
        {
            if (typeTo == "profile")
                stop("Converting IMS data to profile data is not supported.", call. = FALSE)
            if (typeTo == "centroid")
                IMS <- NA
            else if (typeTo == "ims")
                centroid <- FALSE # centroiding is ignored for Bruker data and messes up Agilent data
        }
        
        args <- c(args, list(centroid = centroid, IMS = IMS))
    }
    else if (algorithm == "im_collapse")
        args <- c(args, list(typeFrom = typeFrom))
    else if (algorithm == "timsconvert")
        args <- c(args, list(centroid = typeTo == "centroid", IMS = typeTo == "ims"))
    
    do.call(convertMSFilePaths, args)
}
