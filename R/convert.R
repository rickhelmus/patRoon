convertMSFilesPWiz <- function(inFiles, outFiles, to, filters, extraOpts,
                               logPath, maxProcAmount)
{
    mainArgs <- paste0("--", to)
    if (!is.null(filters))
        mainArgs <- c(mainArgs, sapply(filters, function(f) c("--filter", f)))
    if (!is.null(extraOpts))
        mainArgs <- c(mainArgs, extraOpts)
    
    msc <- getCommandWithOptPath("msconvert", "pwiz")
    cmdQueue <- lapply(seq_along(inFiles), function(fi)
    {
        basef <- basename(tools::file_path_sans_ext(inFiles[fi]))
        
        logf <- if (!is.null(logPath)) file.path(logPath, paste0("pwiz-", basef, ".txt")) else NULL
        logfe <- if (!is.null(logPath)) file.path(logPath, paste0("pwiz-", basef, "-err.txt")) else NULL
        return(list(stdoutFile = logf, stderrFile = logfe, command = msc,
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
        logfe <- if (!is.null(logPath)) file.path(logPath, paste0("openms-", basef, "-err.txt")) else NULL
        
        return(list(stdoutFile = logf, stderrFile = logfe, command = msc,
                    args = c("-in", inFiles[fi], "-out", outFiles[fi], mainArgs)))
    })
    
    executeMultiProcess(cmdQueue, function(cmd) {}, maxProcAmount = maxProcAmount)
    
    invisible(NULL)
}

#' Conversion of analyses between several open and closed data formats
#'
#' This function converts the data format of an analysis to another. It either
#' uses tools from \href{http://proteowizard.sourceforge.net/}{ProteoWizard}
#' (\command{msConvert} command) or \href{http://www.openms.de/}{OpenMS}
#' (\command{FileConverter} command) to perform the conversion. The supported
#' input and output formats include \file{mzXML} and \file{.mzML}. Furthermore,
#' when ProteoWizard is used for conversion, most major (closed) vendor formats
#' are supported for input files.
#'
#' @param files A \code{character} vector with input files. Alternatively, if
#'   \code{dirs=TRUE}, then a \code{character} vector with one or more
#'   directories from which input files are automatically selected.
#' @param outPath A character vector specifying directories that should be used
#'   for the output. Will be re-cycled if necessary. If \code{NULL}, output
#'   directories will be kept the same as the input directories.
#' @param dirs If \code{TRUE} the \code{files} argument specifies directories
#'   from which input files should be selected.
#' @param from One or more input formats (see below). These are used to find
#'   analyses when \code{dirs=TRUE}.
#' @param to Output format: \code{"mzXML"} or \code{"mzML"}.
#' @param overWrite Should existing destination file be overwritten
#'   (\code{TRUE}) or not (\code{FALSE})?
#' @param algorithm Either \code{"pwiz"} (uses \command{msConvert} of
#'   ProteoWizard) or \code{"openms"} (uses \command{FileConverter} of OpenMS).
#' @param filters When \code{algorithm="pwiz"}: a \code{character} vector
#'   specifying one or more filters. Can be used for peak picking to obtain
#'   centroided data (see examples). The elements of the specified vector are
#'   directly passed to the
#'   \code{--filter} option (see
#'   \href{http://proteowizard.sourceforge.net/tools/filters.html}{here})
#' @param extraOpts A \code{character} vector specifying any extra commandline
#'   parameters passed to \command{msConvert} or \command{FileConverter}. Set to
#'   \code{NULL} to ignore. For options: see
#'   \href{http://ftp.mi.fu-berlin.de/pub/OpenMS/release-documentation/html/TOPP_FileConverter.html}{FileConverter}
#'   and
#'   \href{http://proteowizard.sourceforge.net/tools/msconvert.html}{msConvert}.
#'
#' @template multiProc-args
#'
#' @section Conversion formats: Input and output formats include \code{mzXML}
#'   and \code{mzML}. In addition, when ProteoWizard is used, the following
#'   vendor input formats may be chosen (specified with the \code{from} argument):
#'   \itemize{
#'   \item \code{thermo}: Thermo \file{.RAW} files.
#'   \item \code{bruker}: Bruker \file{.d}, \file{.yep}, \file{.baf} and \file{.fid} files.
#'   \item \code{agilent}: Agilent \file{.d} files.
#'   \item \code{ab}: AB Sciex \file{.wiff} files.
#'   \item \code{waters} Waters \file{.RAW} files.
#'   }
#'   Note that the actual supported file formats depend on the ProteoWizard
#'   installation (see
#'   \href{http://proteowizard.sourceforge.net/formats/index.html}{here}).
#'
#' @examples \donttest{
#' # Use FileConverter of OpenMS to convert between open mzXML/mzML format
#' convertMSFiles("standard-1.mzXML", to = "mzML", algorithm = "openms")
#' 
#' # Convert all Thermo .RAW files in the analyses/raw directory to mzML and
#' # store the files in analyses/mzml. During conversion files are centroided by
#' # the peakPicking filter and only MS 1 data is kept.
#' convertMSFiles("analyses/raw", "analyses/mzml", dirs = TRUE, from = "thermo",
#'                filters = c("peakPicking vendor", "msLevel 1"))
#' }
#'
#' @references \insertRef{Rst2016}{patRoon} \cr\cr
#'   \insertRef{Chambers2012}{patRoon}
#'
#' @export
convertMSFiles <- function(files, outPath = NULL, dirs = FALSE,
                           from = c("thermo", "bruker", "agilent", "ab", "waters", "mzXML", "mzML"),
                           to = "mzML", overWrite = FALSE, algorithm = "pwiz",
                           filters = NULL, extraOpts = NULL,
                           logPath = file.path("log", "convert"),
                           maxProcAmount = getOption("patRoon.maxProcAmount"))
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertCharacter(files, min.len = 1, min.chars = 1, add = ac)
    checkmate::assertCharacter(outPath, min.chars = 1, min.len = 1, null.ok = TRUE, add = ac)
    assertCanCreateDirs(outPath, add = ac)
    checkmate::assertLogical(dirs, add = ac)
    checkmate::assertChoice(to, c("mzXML", "mzML"), add = ac) # UNDONE: enough for now?
    checkmate::assertFlag(overWrite, add = ac)
    checkmate::assertChoice(algorithm, c("pwiz", "openms"), add = ac)
    checkmate::assertCharacter(filters, min.chars = 1, null.ok = TRUE, add = ac)
    checkmate::assertList(extraOpts, null.ok = TRUE, add = ac)
    assertMultiProcArgs(logPath, maxProcAmount, add = ac)
    checkmate::reportAssertions(ac)

    if (algorithm == "pwiz")
        from <- checkmate::matchArg(from, c("thermo", "bruker", "agilent", "ab", "waters", "mzXML", "mzML"),
                                    several.ok = TRUE, add = ac)
    else # OpenMS
        from <- checkmate::matchArg(from, c("mzXML", "mzML"), several.ok = TRUE, add = ac)
        
    if (dirs)
    {
        allExts <- list(thermo = ".raw",
                        bruker = c(".d", "yep", "baf", "fid"),
                        agilent = ".d",
                        ab = ".wiff",
                        waters = ".raw",
                        mzXML = ".mzXML",
                        mzML = ".mzML")
        fExt <- unique(unlist(sapply(from, function(f) allExts[[f]], USE.NAMES = FALSE)))
        files <- unique(unlist(sapply(fExt,
                                      function(e) list.files(files,
                                                             pattern = paste0("*\\", e, "$", collapse = " "),
                                                             ignore.case = TRUE))))
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
            convertMSFilesPWiz(files, output, to, filters, extraOpts, logPath, maxProcAmount)
        else # OpenMS
            convertMSFilesOpenMS(files, output, to, extraOpts, logPath, maxProcAmount)
    }
}
