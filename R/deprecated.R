# nocov start

#' Deprecated and renamed functions.
#'
#' Please do not use these functions anymore since they may be removed in the
#' future.
#'
#' @param \dots Passed to successor function.
#'
#' @name patRoon-deprecated
#' @keywords internal
NULL

#' @details \code{reportMD} performs HTML reporting, please use
#'   \code{\link{reportHTML}} instead.
#' @rdname patRoon-deprecated
#' @export
#' @keywords internal
reportMD <- function(...)
{
    .Deprecated("reportHTML")
    reportHTML(...)
}

#' @details \code{exportDAFiles} will export a set of analyses either in
#'   \file{.mzXML} or \file{.mzML} formats.
#'
#' @param format The output format of exported files. Should be either
#'   \code{"mzXML"}, \code{"mzML"} or \code{"mzData"}.
#' @param exportLine Export line spectra (\code{TRUE}) or profile spectra
#'   (\code{FALSE}). Usually line spectra are preferred, since profile spectra
#'   use signficantly more disk space and increase required memory during
#'   processing.
#' @param outPath Character vector of output paths for exported analyses. Will
#'   be recycled if necessary.
#' @param overWrite If \code{TRUE} existing files will be overwritten.
#'
#' @rdname patRoon-deprecated
#' @keywords internal
#' @export
exportDAFiles <- function(anaInfo, format = "mzML", exportLine = TRUE, outPath = anaInfo$path, overWrite = FALSE)
{
    .Deprecated("convertMSFiles")

    ac <- checkmate::makeAssertCollection()
    anaInfo <- assertAndPrepareAnaInfo(anaInfo, "bruker", add = ac)
    checkmate::assertChoice(format, c("mzML", "mzXML", "mzData"), add = ac)
    checkmate::assertFlag(exportLine, add = ac)
    checkmate::assertCharacter(outPath, min.chars = 1, min.len = 1, add = ac)
    checkmate::assertDirectoryExists(outPath, "w", add = ac)
    checkmate::assertFlag(overWrite, add = ac)
    checkmate::reportAssertions(ac)

    outPath <- rep(outPath, length.out = length(anaInfo$path))

    expConstant <- if (format == "mzXML") DAConstants$daMzXML else if (format == "mzData") DAConstants$daMzData else DAConstants$daMzML
    expSpecConstant <- if (exportLine) DAConstants$daLine else DAConstants$daProfile

    DA <- getDAApplication()
    hideDAInScope()

    for (i in seq_len(nrow(anaInfo)))
    {
        printf("Exporting analysis '%s' (%d/%d)... ", anaInfo$analysis[i], i, nrow(anaInfo))

        outf <- file.path(outPath[i], paste0(anaInfo$analysis[i], ".", format))
        if (!overWrite && file.exists(outf))
        {
            cat("Skipped: already exists.\n")
            next
        }

        ind <- getDAFileIndex(DA, anaInfo$analysis[i], anaInfo$path[i])
        if (ind == -1)
        {
            cat("Failed!!")
            next
        }

        DA[["Analyses"]][[ind]]$Export(outf, expConstant, expSpecConstant)

        cat("Done!\n")
    }

    invisible(NULL)
}

#' @details Please use \code{\link{plotChroms}} instead.
#' @rdname patRoon-deprecated
#' @export
#' @keywords internal
plotEIC <- function(obj, ...)
{
    .Deprecated("plotChroms")
    plotChroms(obj, ...)
}

#' @details Please use \code{\link{groupTable}} instead.
#' @rdname patRoon-deprecated
#' @export
#' @keywords internal
groups <- function(object, ...)
{
    .Deprecated("groupTable")
    groupTable(object, ...)
}

#' @details Please use \code{\link{plotSpectrum}} instead.
#' @rdname patRoon-deprecated
#' @export
#' @keywords internal
plotSpec <- function(obj, ...)
{
    .Deprecated("plotSpectrum")
    plotSpectrum(obj, ...)
}

# nocov end
