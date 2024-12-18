# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
NULL

#' Report workflow data
#'
#' Functionality to report data produced by most workflow steps such as features, feature groups, formula and compound
#' annotations, and TPs.
#'
#' The reporting functionality is typically used at the very end of the workflow. It is used to overview the data
#' generated during the workflow, such as features, their annotations and TP screening results.
#' @name reporting
NULL

adjustReportSettings <- function(settings, adjSettings) return(modifyList(settings, adjSettings, keep.null = TRUE))

#' @details \code{report} reports all workflow data in an interactive \acronym{HTML} file. The reports include both
#'   tabular data (\emph{e.g.} retention times, annotation properties, screening results) and varios plots (\emph{e.g.}
#'   chromatograms, (annotated) mass spectra and many more). This function uses functionality from other \R packages,
#'   such as \CRANpkg{rmarkdown}, \CRANpkg{flexdashboard}, \CRANpkg{knitr} and \CRANpkg{bslib}.
#'
#' @param fGroups The \code{\link{featureGroups}} object that should be used for reporting data.
#' @param MSPeakLists,formulas,compounds,compsCluster,components,TPs Further objects (\code{\link{MSPeakLists}},
#'   \code{\link{formulas}}, \code{\link{compounds}}, \code{\link{compoundsCluster}}, \code{\link{components}},
#'   \code{\link{transformationProducts}}) that should be reported. Specify \code{NULL} to skip reporting a particular
#'   object. Note that \code{MSPeakLists} must be set if either \code{formulas} or \code{compounds} is set.
#' @param settingsFile The path to the report settings file used for report configuration (see \verb{Report settings}).
#' @param path The destination file path for files generated during reporting. Will be generated if needed. If
#'   \code{path=NULL} then the destination path is taken from the report settings (see below).
#' @param clearPath If \code{TRUE} then the report destination path will be (recursively) removed prior to reporting.
#' @param openReport If set to \code{TRUE} then the output report file will be opened with the system browser.
#' @param overrideSettings A \code{list} with settings that override those from the report settings file. Example:
#'   \code{overrideSettings=list(compounds=list(topMost=25))}.
#'
#' @template EICParams-arg
#' @template specSimParams-arg
#' @template parallel-arg
#'
#' @section Report settings: The report generation can be customized with a variety of settings that are read from a
#'   \file{YAML} file. This is especially useful if you want to change more advanced settings or want to add or remove
#'   the parts that are reported The report settings file is specified through the \code{settingsFile} argument. If not
#'   specified then default settings will be used. To ease creation of a new template settings file, the
#'   \code{genReportSettingsFile} function can be used.
#'
#'   The following settings are currently available: \itemize{
#'
#'   \item General \itemize{
#'
#'      \item \code{format}: the report format. Currently this can only be \code{"html"}.
#'
#'      \item \code{path}: the destination path (ignored if the \code{path} argument is specified).
#'
#'      \item \code{keepUnusedPlots}: the number of days that unused plot files are kept (see \verb{Plot file caching}).
#'
#'      \item \code{selfContained}: If \code{true} then the output \file{report.html} embeds all graphics and script
#'      dependencies. Otherwise these files are read from the \code{report_files/} directory. Self-contained reports are
#'      easily shared, since only the \file{report.html} needs to be copied. However, they may be slower to generate and
#'      render, especially when the report contains a lot of data.
#'
#'      \item \code{noDate} Set to \code{true} to omit the date from the report. Mainly used for internal purposes.
#'
#'   }
#'
#'   \item \code{summary}: defines the plots on the summary page: \code{chord}, \code{venn} and/or \code{upset}.
#'
#'   \item \code{features} \itemize{
#'
#'      \item \code{retMin}: if \code{true} then retention times are reported in minutes.
#'
#'      \item \code{chromatograms} \itemize{
#'
#'          \item \code{large}: inclusion of large chromatograms (used in feature group table and TP parent chromatogram
#'          view).
#'
#'          \item \code{small}: inclusion of small chromatograms (feature group table).
#'
#'          \item \code{features}: inclusion of chromatograms for individual features (features view). Set to \code{all}
#'          to also include plots for analyses in which a feature was not found (or removed afterwards).
#'          
#'          \item \code{intMax}: Method to determine the maximum intensity plot range: \code{eic} or \code{feature}.
#'          Sets the \code{intMax} argument to \code{plotChroms}.
#'
#'      }
#'
#'      \item \code{intensityPlots}: inclusion of intensity trend plots.
#'
#'   }
#'
#'   \item \code{MSPeakLists} \itemize{
#'
#'      \item \code{spectra}: inclusion of MS and MS/MS spectra (not annotated).
#'
#'   }
#'
#'   \item \code{formulas} \itemize{
#'
#'      \item \code{include}: whether formula results are reported (formula view). If \code{false} then the input
#'      \code{formulas} object is still used to amend \emph{e.g.} compound annotated spectra.
#'
#'      \item \code{normalizeScores}, \code{exclNormScores}: controls score normalization, sets the equally named
#'      arguments to \emph{e.g.} \code{\link{plotScores}}.
#'
#'      \item \code{topMost} only report this number of top ranked candidates. This number can be lowered to speed-up
#'      report generation.
#'
#'   }
#'
#'   \item \code{compounds} \itemize{
#'
#'      \item \code{normalizeScores}, \code{exclNormScores}, \code{topMost}: same as \code{formulas}, see above.
#'
#'   }
#'
#'   \item \code{TPs} \itemize{
#'
#'      \item \code{graphs}: inclusion of TP hierarchy graphs (generated with \code{\link{plotGraph}}).
#'
#'      \item \code{graphStructuresMax}: maximum number of structures to plot in hierarchy graphs (sets
#'      \code{structuresMax} argument of \code{\link[=plotGraph,transformationProductsStructure-method]{plotGraph}}).
#'
#'   }
#'
#'   \item \code{internalStandards} \itemize{
#'
#'      \item \code{graph}: inclusion of internal standard network plot
#'      (\code{\link[=plotGraph,featureGroups-method]{plotGraph}}).
#'
#'   }
#'
#'   }
#'
#' @section Plot file caching: When a new report is generated the plot files are stored inside the \code{report_files}
#'   sub-directory inside the destination path of the report. The plot files are kept so they can be reused to speed-up
#'   re-creation of reports (\emph{e.g.} with different report settings). After the report is generated, any unused plot
#'   files are removed unless they were recently created (controlled by the \code{keepUnusedPlots} setting, see previous
#'   section). The \code{clearPath} argument can be used to completely remove any old files.
#'
#' @note No data will be reported for feature groups in any of the reported objects (\code{formulas}, \code{compounds}
#'   etc) which are \emph{not} present in the input \code{\link{featureGroups}} object (\code{fGroups}).
#'
#'   The \code{topMost}, \code{topMostByRGroup} and \code{onlyPresent} \link[=EICParams]{EIC parameters} may be ignored,
#'   \emph{e.g.}, when generating overview plots.
#'
#' @references Creating MetFrag landing page URLs based on code from
#'   \href{https://github.com/Treutler/MetFamily}{MetFamily} R package. \cr\cr \addCitations{knitr}{2} \cr\cr
#'   \addCitations{knitr}{3}
#'
#' @aliases report
#' @rdname reporting
#' @export
setMethod("report", "featureGroups", function(fGroups, MSPeakLists, formulas, compounds, compsCluster, components,
                                              TPs, settingsFile, path, EICParams, specSimParams, clearPath, openReport,
                                              parallel, overrideSettings)
{
    ac <- checkmate::makeAssertCollection()
    if (!is.null(path))
        checkmate::assertPathForOutput(path, overwrite = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ clearPath + openReport + parallel, fixed = list(add = ac))
    aapply(checkmate::assertClass, . ~ MSPeakLists + formulas + compounds + compsCluster + components + TPs,
           c("MSPeakLists", "formulas", "compounds", "compoundsCluster", "components", "transformationProducts"),
           null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertFileExists(settingsFile, "r", add = ac)
    assertSpecSimParams(specSimParams, add = ac)
    checkmate::assertList(overrideSettings, any.missing = FALSE, names = "unique", add = ac)
    checkmate::reportAssertions(ac)
    
    if (is.null(MSPeakLists) && (!is.null(formulas) || !is.null(compounds)))
        stop("MSPeakLists is NULL, please specify when reporting formula and/or compounds")
    
    settings <- readYAML(settingsFile)
    settings <- adjustReportSettings(settings, overrideSettings)
    settings <- assertAndPrepareReportSettings(settings)
    
    if (is.null(path))
        path <- settings$general$path
    
    if (length(fGroups) == 0)
    {
        cat("No feature groups, nothing to report...\n")
        return(invisible(NULL))
    }
    
    prepareReportPath(path, clearPath)
    
    # UNDONE: check format setting here, when others are supported
    
    doReportHTML(fGroups, MSPeakLists, formulas, compounds, compsCluster, components, TPs, settings, path,
                 EICParams, specSimParams, openReport, parallel)
    
    invisible(NULL)
})

#' @details The \code{genReportSettingsFile} function generates a new template \file{YAML} file to configure report
#'   settings (see the next section).
#'
#' @param out The output file path.
#' @param baseFrom An existing report file to which the report settings should be based from. This is primarily used to
#'   update old settings files: the output settings file will be based on the old settings and amended with any missing.
#'
#' @rdname reporting
#' @export
genReportSettingsFile <- function(out = "report.yml", baseFrom = NULL)
{
    checkmate::assertPathForOutput(out, overwrite = TRUE)
    if (!is.null(baseFrom))
        checkmate::assertFileExists(baseFrom, "r")
    
    defFile <- system.file("report", "settings.yml", package = "patRoon")
    
    if (is.null(baseFrom))
        file.copy(defFile, out, overwrite = TRUE)
    else
    {
        settings <- readYAML(baseFrom)
        settings <- adjustReportSettings(readYAML(defFile), settings)
        settings <- assertAndPrepareReportSettings(settings, setAggr = FALSE)
        writeYAML(settings, out)
    }

    invisible(NULL)
}
