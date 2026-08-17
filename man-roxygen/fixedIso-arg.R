#' @param fixedIsolationWidth Configures how MS/MS spectra are selected for a feature: \itemize{
#'
#'   \item \code{NA}: select all spectra unconditionally, \emph{i.e.} ignoring any precursor information of spectra.
#'
#'   \item \code{FALSE}: select the spectra that were recorded for the feature m/z, using the instrumental precursor
#'   isolation window (\emph{e.g.} quadrupole m/z width) as the selection tolerance.
#'
#'   \item A numeric value: select the spectra that were recorded for the feature m/z within this tolerance window (+/-
#'   precursor m/z).
#'
#'   \item A two-sized \code{numeric} vector: select the spectra that were recorded for the feature m/z within the
#'   specified minimum and maximum tolerance window (relative to the precursor m/z).
#'
#'   }
#'
#'   If no isolation was applied to record MS/MS data (\emph{e.g.} data-independent MS/MS), then all MS/MS spectra will
#'   be always be selected.
#'   
#'   \strong{NOTE}: Sometimes the isolation windows are not exported and cannot be deduced automatically (e.g. Agilent
#'   data). In that case, \code{fixedIsolationWidth} needs to be set to a numeric or \code{NA} value.
