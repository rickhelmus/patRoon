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

MSFileExtensions <- function()
{
    list(thermo = "raw",
         bruker = c("d", "yep", "baf", "fid"),
         agilent = "d",
         ab = "wiff",
         waters = "raw",
         mzXML = "mzXML",
         mzML = "mzML")
}

MSFileFormatIsDir <- function(format, ext)
{
    # UNDONE: is agilent .d also a directory?
    return((format == "bruker" && ext == "d") || (format == "waters" && ext == "raw"))
}

#' @details \code{MSFileFormats} returns a \code{character} with all supported
#'   input formats (see below).
#' @param vendor If \code{TRUE} only vendor formats are returned.
#' @rdname convertMSFiles
#' @export
MSFileFormats <- function(algorithm = "pwiz", vendor = FALSE)
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

filterMSFileDirs <- function(files, from)
{
    if (length(files) == 0)
        return(files)
    
    allFromExts <- MSFileExtensions()[from]
    keep <- sapply(files, function(file)
    {
        fExt <- tools::file_ext(file)
        
        fromExts <- pruneList(lapply(allFromExts, function(f) f[f %in% fExt]), checkEmptyElements = TRUE)
        if (length(fromExts) == 0)
            return(FALSE)
        
        fromCheck <- names(fromExts)
        shouldBeDir <- mapply(fromCheck, fromExts, SIMPLIFY = TRUE,
                              FUN = function(format, exts) sapply(exts, MSFileFormatIsDir, format = format))
        
        if (!allSame(shouldBeDir))
            return(TRUE) # can be either
        
        isDir <- file.info(file, extra_cols = FALSE)$isdir
        if (all(shouldBeDir))
            return(isDir)
        return(!isDir)
    })

    return(files[keep])    
}

listMSFiles <- function(dirs, from)
{
    allExts <- MSFileExtensions()
    allExts <- unique(unlist(allExts[from]))

    files <- list.files(dirs, full.names = TRUE, pattern = paste0("*\\.", allExts, "$", collapse = "|"),
                        ignore.case = TRUE)

    return(filterMSFileDirs(files, from))
}

convertMSFilesPWiz <- function(inFiles, outFiles, to, centroid, filters, extraOpts, PWizBatchSize)
{
    if (centroid != FALSE)
    {
        if (is.null(filters))
            filters <- character()
        filters <- c(filters, paste("peakPicking", if (is.character(centroid)) centroid else ""))
    }

    mainArgs <- paste0("--", to)
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

    msc <- getCommandWithOptPath("FileConverter", "OpenMS")
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
#'   \href{https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_FileConverter.html}{FileConverter}
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
convertMSFiles <- function(files = NULL, outPath = NULL, dirs = TRUE,
                           anaInfo = NULL, from = NULL, to = "mzML",
                           overWrite = FALSE, algorithm = "pwiz",
                           centroid = algorithm != "openms",
                           filters = NULL, extraOpts = NULL, PWizBatchSize = 1)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertCharacter(files, min.len = 1, min.chars = 1, null.ok = !is.null(anaInfo), add = ac)
    checkmate::assertCharacter(outPath, min.chars = 1, min.len = 1, null.ok = TRUE, add = ac)
    assertCanCreateDirs(outPath, add = ac)
    checkmate::assertFlag(dirs, add = ac)
    checkmate::assertChoice(to, c("mzXML", "mzML"), add = ac) # UNDONE: enough for now?
    checkmate::assertFlag(overWrite, add = ac)
    checkmate::assertChoice(algorithm, c("pwiz", "openms", "bruker"), add = ac)
    checkmate::assert(checkmate::checkFlag(centroid),
                      checkmate::checkChoice(centroid, c("vendor", "cwt")),
                      .var.name = centroid)
    checkmate::assertCharacter(filters, min.chars = 1, null.ok = TRUE, add = ac)
    checkmate::assertCharacter(extraOpts, null.ok = TRUE, add = ac)
    checkmate::assertCount(PWizBatchSize, add = ac)
    checkmate::reportAssertions(ac)

    if (centroid != FALSE && algorithm == "openms")
        stop("Centroiding with OpenMS is currently not supported.")
    else if ((centroid == "vendor" || centroid == "cwt") && algorithm != "pwiz")
        stop("Vendor/cwt centroiding is only supported when algorithm=\"pwiz\"")

    if (dirs || !is.null(anaInfo)) # from arg needs to be used?
    {
        if (algorithm == "pwiz")
            from <- checkmate::matchArg(from, c("thermo", "bruker", "agilent", "ab", "waters", "mzXML", "mzML"),
                                        several.ok = FALSE, add = ac)
        else if (algorithm == "openms")
            from <- checkmate::matchArg(from, c("mzXML", "mzML"), several.ok = FALSE, add = ac)
        else # bruker
            from <- checkmate::matchArg(from, "bruker", add = ac)

        if (from == to)
            stop("Input and output formats are the same")
    }

    anaInfo <- assertAndPrepareAnaInfo(anaInfo, from, null.ok = !is.null(files), add = ac)

    if (!is.null(files))
    {
        if (dirs)
        {
            dirs <- files[file.info(files, extra_cols = FALSE)$isdir]

            if (MSFileFormatIsDir(from))
                dirs <- dirs[!grepl(sprintf("(\\.%s)$", MSFileExtensions()[[from]]), files)] # filter out analyses "files" (are actually directories)

            dirFiles <- listMSFiles(dirs, from)
            files <- union(dirFiles, setdiff(files, dirs))
        }
    }
    else
        files <- character()

    if (!is.null(anaInfo))
    {
        ext <- MSFileExtensions()[[from]]
        afiles <- unlist(Map(anaInfo$path, anaInfo$analysis, f = function(p, a) paste0(file.path(p, a), ".", ext)))
        afiles <- afiles[file.exists(afiles)]
        afiles <- filterMSFileDirs(afiles, from)
        files <- c(files, afiles)
    }

    if (is.null(outPath))
        outPath <- dirname(files)

    mkdirp(outPath)

    # NOTE: use normalizePath() here to convert to backslashes on Windows: needed by msconvert
    outPath <- normalizePath(rep(outPath, length.out = length(files)), mustWork = TRUE)
    files <- normalizePath(files, mustWork = FALSE) # no mustWork, file existence will be checked later

    basef <- basename(tools::file_path_sans_ext(files))
    output <- normalizePath(file.path(outPath, paste0(basef, ".", to)),
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
            convertMSFilesPWiz(files, output, to, centroid, filters, extraOpts, PWizBatchSize)
        else if (algorithm == "openms")
            convertMSFilesOpenMS(files, output, to, extraOpts)
        else # bruker
            convertMSFilesBruker(files, output, to, centroid)
    }
}
