# SPDX-FileCopyrightText: 2016 - 2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' MS data conversion
#'
#' Conversion of MS analysis files between several open and closed data formats.
#'
#' @param algorithm Either \code{"pwiz"} (implemented by \command{msConvert} of
#'   ProteoWizard), \code{"openms"} (implemented by \command{FileConverter} of
#'   OpenMS) or \code{"bruker"} (implemented by DataAnalysis).
#'
#' @templateVar what \code{convertMSFiles} (except if \code{algorithm="bruker"})
#' @template uses-multiProc
#'
#'
#' @name convertMSFiles
NULL

#' @export
getMSConversionTypes <- function(algorithm, direction)
{
    assertMSConvAlgo(algorithm)
    checkmate::assertChoice(direction, c("input", "output"))
 
    if (direction == "input")
    {
        return(switch(algo,
                      pwiz = getMSFileTypes(),
                      openms = c("centroid", "profile"),
                      bruker = "raw",
                      im_collapse = c("raw", "ims"),
                      timsconvert = "raw"))
    }

    return(switch(algo,
                  pwiz = c("centroid", "profile", "ims"),
                  openms = c("centroid"),
                  bruker = c("centroid", "profile"),
                  im_collapse = "centroid",
                  timsconvert = c("centroid", "profile", "ims")))
}

#' @details \code{getMSInConversionFormats} returns a \code{character} with all supported
#'   input formats (see below).
#' @param vendor If \code{TRUE} only vendor formats are returned.
#' @rdname convertMSFiles
#' @export
getMSConversionFormats <- function(algorithm, direction, type = NULL)
{
    assertMSConvAlgo(algorithm)
    checkmate::assertChoice(direction, c("input", "output"))
    checkmate::assertChoice(type, getMSConversionTypes(algorithm, direction), null.ok = TRUE)

    ret <- if (direction == "input")
    {
        if (algorithm == "openms")
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

#' @export
convertMSFilesPWiz <- function(inFiles, outFiles, formatTo = "mzML", centroid = TRUE, IMS = FALSE, minIntensity = 5,
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

#' @export
convertMSFilesIMSCollapse <- function(inFiles, outFiles, typeFrom, formatTo = "mzML", mzRange = NULL, mobilityRange = NULL,
                                      clMethod = "distance", mzWindow = defaultLim("mz", "medium"),
                                      minAbundanceRel = 0, minAbundanceAbs = 0, topMost = NULL,
                                      minIntensityIMS = NULL, minIntensityPre = NULL, includeMSMS = FALSE, ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertCharacter(inFiles, min.chars = 1, min.len = 1, add = ac)
    checkmate::assertCharacter(outFiles, min.chars = 1, len = length(inFiles), add = ac)
    checkmate::assertChoice(typeFrom, c("raw", "ims"), add = ac)
    checkmate::assertChoice(formatTo, getMSConversionFormats("im_collapse", "output"), add = ac)
    aapply(assertRange, . ~ mzRange + mobilityRange, null.ok = TRUE, fixed = list(add = ac))
    assertClusterMethod(clMethod, add = ac)
    aapply(checkmate::assertNumber, . ~ mzWindow + minAbundanceRel + minAbundanceAbs + minIntensityIMS + minIntensityPre,
           lower = 0, finite = TRUE, null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertCount(topMost, positive = TRUE, null.ok = TRUE, add = ac)
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
            isolationWindowLowerOffset = if (MSLevel == 1) NA_real_ else collapsedSpectra$precursorMZs - collapsedSpectra$isolationStarts,
            isolationWindowUpperOffset = if (MSLevel == 1) NA_real_ else collapsedSpectra$isolationEnds - collapsedSpectra$precursorMZs,
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
    applyMSData(anaInfo, anaInfo$path_centroid, needIMS = TRUE, showProgress = TRUE, func = function(ana, path, backend, outd)
    {
        outp <- file.path(outd, paste0(ana, ".", formatTo))
        
        openMSReadBackend(backend, path)
        collapsedSpectra <- collapseIMSFrames(backend, NULLToZero(mzRange[1]), NULLToZero(mzRange[2]),
                                              NULLToZero(mobilityRange[1]), NULLToZero(mobilityRange[2]), clMethod,
                                              mzWindow, minAbundanceRel, minAbundanceAbs, NULLToZero(topMost),
                                              NULLToZero(minIntensityIMS), NULLToZero(minIntensityPre), includeMSMS)
        
        
        
        headerMS1 <- makeHeader(collapsedSpectra, 1, getMSMetadata(backend, 1))
        
        if (!includeMSMS)
            header <- headerMS1
        else
        {
            headerMS2 <- makeHeader(collapsedSpectra, 2, getMSMetadata(backend, 2))
            header <- rbind(headerMS1, headerMS2)
        }
        
        ord <- order(header$frame, header$precursorMZ)
        header <- header[ord, ]
        header$seqNum <- header$acquisitionNum <- seq_len(nrow(header))
        header$spectrumId <- sprintf("merged=%d frame=%d", header$seqNum, header$frame)
        header <- header[, names(header) != "frame"]

        allSpectra <- c(collapsedSpectra$MS1, collapsedSpectra$MS2)
        allSpectra <- allSpectra[ord]

        mzR::writeMSData(allSpectra, outp, header, outformat = formatTo, ...)

        doProgress()
    })
    
    invisible(NULL)
}

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
        reticulate::use_virtualenv(virtualenv)
    
    mainArgs <- character()
    if (IMS)
        mainArgs <- c(mainArgs, "--mode", "raw")
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

#' @details \code{convertMSFiles} converts the data format of an analysis to
#'   another. It uses tools from
#'   \href{http://proteowizard.sourceforge.net/}{ProteoWizard}
#'   (\command{msConvert} command), \href{http://www.openms.de/}{OpenMS}
#'   (\command{FileConverter} command) or Bruker DataAnalysis to perform the
#'   conversion. Supported input and output formats include \file{mzXML},
#'   \file{.mzML} and several vendor formats, depending on which algorithm is
#'   used.
#'
#' @param files,dirs The \code{files} argument should be a \code{character}
#'   vector with input files. If \code{files} contains directories and
#'   \code{dirs=TRUE} then files from these directories are also considered. An
#'   alternative method to specify input files is by the \code{anaInfo}
#'   argument. If the latter is specified \code{files} may be \code{NULL}.
#' @param outPath A character vector specifying directories that should be used
#'   for the output. Will be re-cycled if necessary. If \code{NULL}, output
#'   directories will be kept the same as the input directories.
#' @param anaInfo An \link[=analysis-information]{analysis info table} used to
#'   retrieve input files. Either this argument or \code{files} (or both) should
#'   be set (\emph{i.e.} not \code{NULL}).
#' @param from Input format (see below). These are used to find analyses when
#'   \code{dirs=TRUE} or \code{anaInfo} is set.
#' @param to Output format: \code{"mzXML"} or \code{"mzML"}.
#' @param overWrite Should existing destination file be overwritten
#'   (\code{TRUE}) or not (\code{FALSE})?
#' @param centroid Set to \code{TRUE} to enable centroiding (not supported if
#'   \code{algorithm="openms"}). In addition, when \code{algorithm="pwiz"} the
#'   value may be \code{"vendor"} to perform centroiding with the vendor
#'   algorithm or \code{"cwt"} to use ProteoWizard's wavelet algorithm.
#' @param filters When \code{algorithm="pwiz"}: a \code{character} vector
#'   specifying one or more filters. The elements of the specified vector are
#'   directly passed to the \code{--filter} option (see
#'   \href{http://proteowizard.sourceforge.net/tools/filters.html}{here})
#' @param extraOpts A \code{character} vector specifying any extra commandline
#'   parameters passed to \command{msConvert} or \command{FileConverter}. Set to
#'   \code{NULL} to ignore. For options: see
#'   \href{https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_FileConverter.html}{FileConverter}
#'    and
#'   \href{http://proteowizard.sourceforge.net/tools/msconvert.html}{msConvert}.
#' @param PWizBatchSize When \code{algorithm="pwiz"}: the number of analyses to
#'   process by a single call to \command{msConvert}. Usually a value of one is
#'   most efficient. Set to zero to run all analyses all at once from a single
#'   call.
#'
#' @section Conversion formats: Possible output formats (\code{to} argument) are
#'   \code{mzXML} and \code{mzML}.
#'
#'   Possible input formats (\code{from} argument) depend on the algorithm that
#'   was chosen and may include:
#'
#'   \itemize{
#'
#'   \item \code{thermo}: Thermo \file{.RAW} files (only
#'   \code{algorithm="pwiz"}).
#'
#'   \item \code{bruker}: Bruker \file{.d}, \file{.yep}, \file{.baf} and
#'   \file{.fid} files (only \code{algorithm="pwiz"} or
#'   \code{algorithm="bruker"}).
#'
#'   \item \code{agilent}: Agilent \file{.d} files (only
#'   \code{algorithm="pwiz"}).
#'
#'   \item \code{ab}: AB Sciex \file{.wiff} files (only
#'   \code{algorithm="pwiz"}).
#'
#'   \item \code{waters} Waters \file{.RAW} files (only
#'   \code{algorithm="pwiz"}).
#'
#'   \item \code{mzXML}/\code{mzML}: Open format \file{.mzXML}/\file{.mzML}
#'   files (only \code{algorithm="pwiz"} or \code{algorithm="openms"}).
#'
#'   }
#'
#'   Note that the actual supported file formats of ProteoWizard depend on how
#'   it was installed (see
#'   \href{http://proteowizard.sourceforge.net/formats/index.html}{here}).
#'
#' @examples \dontrun{
#' # Use FileConverter of OpenMS to convert between open mzXML/mzML format
#' convertMSFiles("standard-1.mzXML", to = "mzML", algorithm = "openms")
#'
#' # Convert all Thermo .RAW files in the analyses/raw directory to mzML and
#' # store the files in analyses/mzml. During conversion files are centroided by
#' # the peakPicking filter and only MS 1 data is kept.
#' convertMSFiles("analyses/raw", "analyses/mzml", dirs = TRUE, from = "thermo",
#'                centroid = "vendor", filters = "msLevel 1")
#' }
#'
#' @references \insertRef{Rst2016}{patRoon} \cr\cr
#'   \insertRef{Chambers2012}{patRoon}
#'
#' @rdname convertMSFiles
#' @export
convertMSFilePaths <- function(files, formatFrom, formatTo = "mzML", outPath = NULL, dirs = TRUE, overWrite = FALSE,
                               algorithm = "pwiz", ...)
{
    ac <- checkmate::makeAssertCollection()
    assertConvertMSFilesArgs(formatFrom, formatTo, overWrite, algorithm, add = ac)
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
        else if (!overWrite && file.exists(output[fi]))
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

#' @export
convertMSFilesAnaInfo <- function(anaInfo, typeFrom = "raw", typeTo = "centroid", formatFrom, formatTo = "mzML",
                                  overWrite = FALSE, algorithm = "pwiz", centroidVendor = TRUE, ...)
{
    ac <- checkmate::makeAssertCollection()
    assertConvertMSFilesArgs(formatFrom, formatTo, overWrite, algorithm, add = ac)
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
                 overWrite = overWrite, algorithm = algorithm, ...)

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
