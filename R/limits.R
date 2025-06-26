#' Default limits and tolerances
#'
#' Get and configure default limits and tolerances used by \pkg{patRoon}.
#'
#' Tolerances for retention times, \emph{m/z} values and other numerical limits and tolerances are widely used in
#' \pkg{patRoon} to process data. Their defaults used to be hardcoded directly as defaults to function arguments. Since
#' version 3.0 the defaults are centralized in a \file{limits YAML} file. This simplifies their configuration and makes
#' it easier to switch the defaults between different HRMS instruments.
#'
#' The limits configuration file is taken from one of the following locations (in order): \itemize{
#'
#'   \item the path specified in the \code{patRoon.path.limits} option (if specified)
#'
#'   \item the \file{limits.yml} file in the current working directory (if present)
#'
#'   \item the default \file{limits.yml} file embedded in \pkg{patRoon}
#'
#' }
#'
#' @section limits YAML file format: The limits are configured in a simple \file{YAML} file. A brief summary of the
#'   format is given here.
#'
#'   The \strong{general} section has the \code{IMS} field which specifies the type of instrument used. Currently this
#'   is either \code{bruker} or \code{agilent}.
#'
#'   The next sections describe absolute and relative (suffixed by \code{_rel}) limits and tolerances for retention,
#'   \emph{m/z}, mobility and \acronym{CCS} data. These are divided into several tolerance levels: \code{very_narrow},
#'   \code{narrow}, \code{medium} and \code{wide}. In general, the \code{very narrow} tolerances are used to compare
#'   data that should be equivalent, but may be slightly different due to \emph{e.g.} small rounding errors. The
#'   \code{narrow} tolerances are generally sufficient when only small deviations are expected, \emph{e.g.} comparing
#'   different \emph{m/z} values from well calibrated HRMS instruments. The \code{medium} tolerances are generally used
#'   when a slightly larger tolerance is needed, \emph{e.g.} when \emph{m/z} values from raw (non-averaged) spectra. The
#'   \code{wide} values are mainly used for plot limits, \emph{e.g.} to provide a reasonable zoom-out. The sections only
#'   define values for the tolerance levels actually used in \pkg{patRoon}.
#'
#'   The \code{mobility_bruker} and \code{mobility_agilent} sections specify the default mobility limits for Bruker and
#'   Agilent systems, respectively. Which are used is set by the \code{IMS} variable in the \strong{general} section.
#'
#'   To see which limits are used in which functions, please refer to the \verb{Usage} section of the respective
#'   functions, specifically how the \code{defaultLim} function is used to assign function argument defaults.
#'
#'   \strong{NOTE}: the choice between using the \code{narrow} and \code{medium} tolerance is not always clear, and
#'   there is some inconsistency in its use throughout \pkg{patRoon} (primarily due to legacy code or keeping defaults
#'   from external algorithms).
#'
#' @eval paste0("@@section Default limits YAML file: \\preformatted{", patRoon:::readAllFile(system.file("misc",
#'   "limits.yml", package = "patRoon")), "}")
#'
#' @note Most of the defaults were derived with Bruker TOF HRMS instrumentation in mind, but should be reasonable with
#'   minor adjustments for others.
#'
#' @name limits
NULL

getLimitsFileClosure <- function()
{
    # NOTE: under certain circumstances, the system.file() call may be slow (noticed during tests on macOS), so we cache
    # the result here
    pkgCachedFile <- system.file("misc", "limits.yml", package = "patRoon")
    function()
    {
        for (p in c(getOption("patRoon.path.limits", ""), file.path(getwd(), "limits.yml"), pkgCachedFile))
        {
            if (!nzchar(p))
                next
            if (file.exists(p))
                return(p)
        }
        stop("Cannot find limits configuration file!")
    }
}

getLimitsFile <- getLimitsFileClosure()    

#' @details \code{defaultLim} returns the limits for a specific category and tolerance level.
#' @param category,level The category and level of the limit to be returned. See the detail sections below. For mobility
#'   related limits, the \code{category="mobility"} (instead of instrument specific categories) should be used.
#' @rdname limits
#' @export
defaultLim <- function(category, level) doDefaultLim(category, level)

defaultLimClosure <- function()
{
    cachedLimits <- cachedFileHash <- NULL
    function(category, level)
    {
        limFile <- getLimitsFile()
        limHash <- makeFileHash(limFile)
        if (is.null(cachedLimits) || is.null(cachedFileHash) || cachedFileHash != limHash)
        {
            cachedLimits <<- readYAML(limFile)
            cachedLimits <<- assertAndPrepareLimits(cachedLimits)
            cachedFileHash <<- makeFileHash(limFile)
        }
        
        checkmate::assertChoice(category, setdiff(names(cachedLimits), "general"))
        checkmate::assertChoice(level, names(cachedLimits[[category]]))
        
        return(cachedLimits[[category]][[level]])
    }
}

doDefaultLim <- defaultLimClosure()

#' @details \code{genLimitsFile} generates a new \file{limits.yml} configuration file in the specified path. The file is
#'   created with the defaults embedded in \pkg{patRoon} (see details below). Generating a custom limits is primarily
#'   useful for non-Bruker IMS workflows.
#'
#' @param outPath The output directory where the new \file{limits.yml} file will be created.
#' @param IMS The IMS instrument type. This sets the \code{IMS} variable in the \strong{general} section. Valid values
#'   are \code{"bruker"} and \code{"agilent"}.
#'
#' @rdname limits
#' @export
genLimitsFile <- function(outPath = normalizePath("."), IMS = "bruker")
{
    checkmate::assertString(outPath, min.chars = 1)
    checkmate::assertChoice(IMS, c("bruker", "agilent"))
    
    fp <- file.path(outPath, "limits.yml")
    checkmate::assertPathForOutput(fp, overwrite = TRUE)
    
    limits <- readYAML(system.file("misc", "limits.yml", package = "patRoon"))
    limits$general$IMS <- IMS
    writeYAML(limits, fp)

    return(invisible(NULL))    
}
