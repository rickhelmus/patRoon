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

#' @details \code{getMSFileConversionFormats} returns a \code{character} with all supported
#'   input formats (see below).
#' @param vendor If \code{TRUE} only vendor formats are returned.
#' @rdname convertMSFiles
#' @export
getMSFileConversionFormats <- function(algorithm = "pwiz", vendor = FALSE)
{
    checkmate::assertChoice(algorithm, c("pwiz", "openms", "bruker"))
    checkmate::assertFlag(vendor)

    if (algorithm == "pwiz")
        ret <- names(MSFileExtensions())
    else if (algorithm == "openms")
        ret <- c("mzXML", "mzML")
    else # algorithm == "bruker"
        ret <- "bruker"

    if (vendor)
        ret <- setdiff(ret, c("mzXML", "mzML"))

    return(ret)
}

convertMSFilesPWiz <- function(inFiles, outFiles, to, centroid, IMS, filters, extraOpts, PWizBatchSize)
{
    if (centroid != FALSE)
    {
        if (is.null(filters))
            filters <- character()
        filters <- c(paste("peakPicking", if (is.character(centroid)) centroid else ""), filters)
    }
    
    if (is.na(IMS))
        filters <- c(filters, "scanSumming sumMs1=1")

    mainArgs <- paste0("--", to)
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

convertMSFilesOpenMS <- function(inFiles, outFiles, to, extraOpts)
{
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

convertMSFilesBruker <- function(inFiles, outFiles, to, centroid)
{
    # expConstant <- if (to == "mzXML") DAConstants$daMzXML else if (to == "mzData") DAConstants$daMzData else DAConstants$daMzML
    expConstant <- if (to == "mzXML") DAConstants$daMzXML else DAConstants$daMzML
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

collapseIMSFiles <- function(anaInfo, mzRange = NULL, mobilityRange = NULL, clMethod = "bin", mzWindow = 0.005,
                             minAbundance = 0, topMost = NULL, minIntensityIMS = NULL, minIntensityPre = NULL,
                             overWrite = FALSE)
{
    # UNDONE: checkmate
    # UNDONE: which default clust method?
    
    anaInfo <- assertAndPrepareAnaInfo(anaInfo)
    
    outDirs <- getPathsFromAnaInfo(anaInfo, "centroid")
    if (is.null(outDirs) || anyNA(outDirs) || any(!nzchar(outDirs)))
        stop("Please set a valid centroid path for all analyses in the analysis information.", call. = FALSE)
    mkdirp(outDirs)
    
    printf("Collapsing all %d analyses ...\n", nrow(anaInfo))
    applyMSData(anaInfo, outDirs, needIMS = TRUE, showProgress = FALSE, func = function(ana, path, backend, outd)
    {
        outp <- file.path(outd, paste0(ana, ".mzML"))
        
        printf("%s --> %s ... ", path, outp)
        
        if (!overWrite && file.exists(outp))
        {
            printf("skipped: already exists\n")
            return(NULL)
        }
        
        openMSReadBackend(backend, path)
        
        meta <- getMSMetadata(backend, 1)
        spectra <- collapseIMSFrames(backend, NULLToZero(mzRange[1]), NULLToZero(mzRange[2]), NULLToZero(mobilityRange[1]),
                                     NULLToZero(mobilityRange[2]), clMethod, mzWindow, minAbundance, NULLToZero(topMost),
                                     NULLToZero(minIntensityIMS), NULLToZero(minIntensityPre))

        specMZRange <- range(unlist(lapply(spectra, function(sp) range(sp[, "mz"]))))
        
        header <- data.frame(seqNum = seq_along(spectra),
                             acquisitionNum = meta$scan,
                             msLevel = 1,
                             polarity = fifelse(meta$scan == 1, 1, 0),
                             peaksCount = sapply(spectra, nrow), # UNDONE: correct?
                             totIonCurrent = sapply(spectra, function(sp) sum(sp[, "intensity"])),
                             retentionTime = meta$time,
                             basePeakMZ = sapply(spectra, function(sp) sp[which.max(sp[, "intensity"]), "mz"]),
                             basePeakIntensity = sapply(spectra, function(sp) max(sp[, "intensity"])),
                             collisionEnergy = NA_real_,
                             ionisationEnergy = 0,
                             lowMZ = 0,
                             highMZ = 0,
                             precursorScanNum = NA_integer_,
                             precursorMZ = NA_real_,
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
                             isolationWindowLowerOffset = NA_real_,
                             isolationWindowUpperOffset = NA_real_,
                             # UNDONE: below are technically not scan limits, but might be good enough?
                             scanWindowLowerLimit = specMZRange[1],
                             scanWindowUpperLimit = specMZRange[2])
        header$spectrumId <- paste0("scan=", meta$scan)
        
        mzR::writeMSData(spectra, outp, header)

        printf("done!\n")
    })
    
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
                               algorithm = "pwiz", centroid = algorithm != "openms", IMS = FALSE, filters = NULL,
                               extraOpts = NULL, PWizBatchSize = 1)
{
    ac <- checkmate::makeAssertCollection()
    assertConvertMSFilesArgs(formatFrom, formatTo, overWrite, algorithm, filters, extraOpts, PWizBatchSize, add = ac)
    checkmate::assertCharacter(files, min.len = 1, min.chars = 1, add = ac)
    checkmate::assertCharacter(outPath, min.chars = 1, min.len = 1, null.ok = TRUE, add = ac)
    assertCanCreateDirs(outPath, add = ac)
    checkmate::assertFlag(dirs, add = ac)
    checkmate::assert(checkmate::checkFlag(centroid),
                      checkmate::checkChoice(centroid, c("vendor", "cwt")),
                      .var.name = "centroid", add = ac)
    checkmate::assertFlag(IMS, na.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (centroid != FALSE && algorithm == "openms")
        stop("Centroiding with OpenMS is currently not supported.", call. = FALSE)
    else if ((centroid == "vendor" || centroid == "cwt") && algorithm != "pwiz")
        stop("Vendor/cwt centroiding is only supported when algorithm=\"pwiz\"", call. = FALSE)
    if (!isFALSE(IMS) && algorithm != "pwiz")
        stop("Handling of IMS data is only supported when algorithm=\"pwiz\"", call. = FALSE)
    
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
    output <- normalizePath(file.path(outPath, paste0(basef, ".", formatTo)),
                            mustWork = FALSE)
    
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
            convertMSFilesPWiz(files, output, formatTo, centroid, IMS, filters, extraOpts, PWizBatchSize)
        else if (algorithm == "openms")
            convertMSFilesOpenMS(files, output, formatTo, extraOpts)
        else # bruker
            convertMSFilesBruker(files, output, formatTo, centroid)
    }
}

#' @export
convertMSFilesAnaInfo <- function(anaInfo, typeFrom = "raw", typeTo = "centroid", formatFrom, formatTo = "mzML",
                                  overWrite = FALSE, algorithm = "pwiz", centroidVendor = TRUE, filters = NULL,
                                  extraOpts = NULL, PWizBatchSize = 1)
{
    validFromTypes <- switch(algorithm,
                             pwiz = getMSFileTypes(),
                             openms = c("centroid", "profile"),
                             bruker = "raw")
    validToTypes <- switch(algorithm,
                           pwiz = c("ims", "profile", "centroid"),
                           openms = c("centroid", "profile"),
                           bruker = c("centroid", "profile"))
    
    ac <- checkmate::makeAssertCollection()
    assertConvertMSFilesArgs(formatFrom, formatTo, overWrite, algorithm, filters, extraOpts, PWizBatchSize, add = ac)
    checkmate::assertChoice(typeFrom, validFromTypes, add = ac)
    checkmate::assertChoice(typeTo, validToTypes, add = ac)
    checkmate::assertFlag(centroidVendor, add = ac)
    checkmate::reportAssertions(ac)

    anaInfo <- assertAndPrepareAnaInfo(anaInfo, typeFrom, formatFrom)
    
    outPath <- getPathsFromAnaInfo(anaInfo, typeTo)
    if (is.null(outPath) || anyNA(outPath) || any(!nzchar(outPath)))
        stop("Please properly configure the analysis info for the output type ", typeTo, call. = FALSE)
    
    files <- getMSFilesFromAnaInfo(anaInfo, typeFrom, formatFrom)

    centroid <- if (typeTo == "profile")
        FALSE
    else if (centroidVendor)
        "vendor"
    else
        TRUE
    
    IMS <- formatFrom %in% c("agilent_ims", "bruker_ims")
    if (IMS)
    {
        if (typeTo == "profile")
            stop("Converting IMS data to profile data is not supported.", call. = FALSE)
        if (typeTo == "centroid")
            IMS <- NA
    }
    
    convertMSFilePaths(files, outPath, dirs = FALSE, formatFrom = formatFrom, formatTo = formatTo,
                       overWrite = overWrite, algorithm = algorithm, centroid = centroid, IMS = IMS, filters = filters,
                       extraOpts = extraOpts, PWizBatchSize = PWizBatchSize)
}
