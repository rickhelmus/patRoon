#' MS data conversion
#'
#' Conversion of MS analysis files between several open and closed data formats.
#'
#' @param algorithm Either \code{"pwiz"} (implemented by \command{msConvert} of
#'   ProteoWizard), \code{"openms"} (implemented by \command{FileConverter} of
#'   OpenMS) or \code{"bruker"} (implemented by DataAnalysis).
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
    isDir <- file.info(files, extra_cols = FALSE)$isdir
    if ("bruker" %in% from)
    {
        isD <- grepl("(\\.d)$", files)
        # filter out directories unless they end with .d
        files <- files[isD | !isDir]

        # filter out any non directory files that end with .d
        files <- files[!isD | isDir]
    }
    else
        files <- files[!isDir]
}

listMSFiles <- function(dirs, from)
{
    allExts <- MSFileExtensions()
    allExts <- unique(unlist(allExts[from]))

    files <- list.files(dirs, full.names = TRUE, pattern = paste0("*\\.", allExts, "$", collapse = "|"),
                        ignore.case = TRUE)

    return(filterMSFileDirs(files, from))
}

convertMSFilesPWiz <- function(inFiles, outFiles, to, centroid, filters, extraOpts,
                               logPath, maxProcAmount)
{
    if (centroid != FALSE)
    {
        if (is.null(filters))
            filters <- character()
        filters <- c(filters, if (centroid == "vendor") "peakPicking vendor" else "peakPicking")
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

    cmdQueue <- lapply(seq_along(inFiles), function(fi)
    {
        basef <- basename(tools::file_path_sans_ext(inFiles[fi]))
        logf <- if (!is.null(logPath)) file.path(logPath, paste0("pwiz-", basef, ".txt")) else NULL
        return(list(logFile = logf, command = msc,
                    args = c(inFiles[fi], "--outfile", outFiles[fi],
                             "-o", dirname(outFiles[fi]), mainArgs)))
    })

    executeMultiProcess(cmdQueue, function(cmd) {}, maxProcAmount = maxProcAmount)

    invisible(NULL)
}

convertMSFilesOpenMS <- function(inFiles, outFiles, to, extraOpts, logPath, maxProcAmount)
{
    mainArgs <- c()
    if (!is.null(extraOpts))
        mainArgs <- c(mainArgs, extraOpts)

    msc <- getCommandWithOptPath("FileConverter", "OpenMS")
    cmdQueue <- lapply(seq_along(inFiles), function(fi)
    {
        basef <- basename(tools::file_path_sans_ext(inFiles[fi]))
        logf <- if (!is.null(logPath)) file.path(logPath, paste0("openms-", basef, ".txt")) else NULL
        return(list(logFile = logf, command = msc,
                    args = c("-in", inFiles[fi], "-out", outFiles[fi], mainArgs)))
    })

    executeMultiProcess(cmdQueue, function(cmd) {}, maxProcAmount = maxProcAmount)

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
#' @param from One or more input formats (see below). These are used to find
#'   analyses when \code{dirs=TRUE} or \code{anaInfo} is set.
#' @param to Output format: \code{"mzXML"} or \code{"mzML"}.
#' @param overWrite Should existing destination file be overwritten
#'   (\code{TRUE}) or not (\code{FALSE})?
#' @param centroid Set to \code{TRUE} to enable centroiding (not supported if
#'   \code{algorithm="openms"}). In addition, when \code{algorithm="pwiz"} the
#'   value may be \code{"vendor"} to perform centroiding with the vendor
#'   algorithm.
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
#'
#' @template multiProc-args
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
                           anaInfo = NULL, from = MSFileFormats(algorithm, vendor = algorithm != "openms"), to = "mzML",
                           overWrite = FALSE, algorithm = "pwiz",
                           centroid = algorithm != "openms",
                           filters = NULL, extraOpts = NULL,
                           logPath = file.path("log", "convert"),
                           maxProcAmount = getOption("patRoon.maxProcAmount"))
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
                      checkmate::checkSetEqual(centroid, "vendor"),
                      .var.name = centroid)
    checkmate::assertCharacter(filters, min.chars = 1, null.ok = TRUE, add = ac)
    checkmate::assertCharacter(extraOpts, null.ok = TRUE, add = ac)
    assertMultiProcArgs(logPath, maxProcAmount, add = ac)
    checkmate::reportAssertions(ac)

    if (centroid != FALSE && algorithm == "openms")
        stop("Centroiding with OpenMS is currently not supported.")
    else if (centroid == "vendor" && algorithm != "pwiz")
        stop("Vendor centroiding is only supported when algorithm=\"pwiz\"")

    if (dirs || !is.null(anaInfo)) # from arg needs to be used?
    {
        if (algorithm == "pwiz")
            from <- checkmate::matchArg(from, c("thermo", "bruker", "agilent", "ab", "waters", "mzXML", "mzML"),
                                        several.ok = TRUE, add = ac)
        else if (algorithm == "openms")
            from <- checkmate::matchArg(from, c("mzXML", "mzML"), several.ok = TRUE, add = ac)
        else # bruker
            from <- checkmate::matchArg(from, "bruker", add = ac)

        ofrom <- from
        from <- setdiff(from, to)
        if (length(from) < length(ofrom))
            warning(paste("Skipping input formats that are also specified as output: ",
                          paste0(setdiff(ofrom, from), collapse = ", ")))
        if (length(from) == 0)
            stop("No (valid) input formats specified.")
    }

    anaInfo <- assertAndPrepareAnaInfo(anaInfo, from, null.ok = !is.null(files), add = ac)

    if (!is.null(files))
    {
        if (dirs)
        {
            dirs <- files[file.info(files, extra_cols = FALSE)$isdir]

            # UNDONE: is agilent .d also a directory?
            if ("bruker" %in% from)
                dirs <- dirs[!grepl("(\\.d)$", files)] # filter out .d analyses "files" (are actually directories)

            dirFiles <- listMSFiles(dirs, from)
            files <- union(dirFiles, setdiff(files, dirs))
        }
    }
    else
        files <- character()

    if (!is.null(anaInfo))
    {
        fExts <- unique(unlist(MSFileExtensions()[from]))
        afiles <- unlist(Map(anaInfo$path, anaInfo$analysis, f = function(p, a) paste0(file.path(p, a), ".", fExts)))
        afiles <- afiles[file.exists(afiles)]
        afiles <- filterMSFileDirs(afiles, from)
        files <- c(files, afiles)
    }

    if (is.null(outPath))
        outPath <- dirname(files)

    mkdirp(outPath)
    if (!is.null(logPath))
        mkdirp(logPath)

    # NOTE: use normalizePath() here to convert to backslashes on Windows: needed by msconvert
    outPath <- normalizePath(rep(outPath, length.out = length(files)), mustWork = TRUE)
    files <- normalizePath(files, mustWork = TRUE)

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

    files <- files[keepFiles]
    output <- output[keepFiles]

    if (length(files) > 0)
    {
        if (algorithm == "pwiz")
            convertMSFilesPWiz(files, output, to, centroid, filters, extraOpts, logPath, maxProcAmount)
        else if (algorithm == "openms")
            convertMSFilesOpenMS(files, output, to, extraOpts, logPath, maxProcAmount)
        else # bruker
            convertMSFilesBruker(files, output, to, centroid)
    }
}
