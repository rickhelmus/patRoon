### Feature groups

setGeneric("groups", function(object) standardGeneric("groups"))
setGeneric("getFeatures", function(fGroups) standardGeneric("getFeatures"))
setGeneric("groupFeatIndex", function(fGroups) standardGeneric("groupFeatIndex"))
setGeneric("groupInfo", function(fGroups) standardGeneric("groupInfo"))
setGeneric("removeAnalyses", function(fGroups, indices) standardGeneric("removeAnalyses"))
setGeneric("removeGroups", function(fGroups, indices) standardGeneric("removeGroups"))
setGeneric("removeEmptyGroups", function(fGroups) standardGeneric("removeEmptyGroups"))
setGeneric("removeEmptyAnalyses", function(fGroups) standardGeneric("removeEmptyAnalyses"))
setGeneric("averageGroups", function(fGroups) standardGeneric("averageGroups"))
setGeneric("updateFeatIndex", function(fGroups) standardGeneric("updateFeatIndex"))
setGeneric("export", function(fGroups, type, out) standardGeneric("export"))
setGeneric("groupTable", function(fGroups, average = FALSE) standardGeneric("groupTable"))
setGeneric("overlap", function(fGroups, which, exclusive = FALSE) standardGeneric("overlap"))
setGeneric("comparison", function(..., groupAlgo,
                                  groupArgs = list(rtalign = FALSE)) standardGeneric("comparison"), signature = "...")
setGeneric("regression", function(fGroups) standardGeneric("regression"))
setGeneric("groupFeatures", function(feat, algorithm, ...) standardGeneric("groupFeatures"))
setGeneric("generateConsensusXML", function(feat, out, rtalign, QT, maxAlignRT, maxAlignMZ, maxGroupRT, maxGroupMZ) standardGeneric("generateConsensusXML"))
setGeneric("groupFeaturesScreening", function(fGroups, scr) standardGeneric("groupFeaturesScreening"))
setGeneric("replicateGroupSubtract", function(fGroups, rGroups, threshold = 0) standardGeneric("replicateGroupSubtract"))

### utils (XCMS)

#' Conversion to xcmsSet objects
#'
#' Converts a \code{\link{features}} or \code{\link{featureGroups}} object to an
#' \code{\link{xcmsSet}} object.
#'
#' @param obj The object that should be converted.
#' @param exportedData,\dots Set to \code{TRUE} if analyses were exported as
#'   \code{mzXML} or \code{mzML} files (ignored for \code{featuresOpenMS} and
#'   \code{featuresXCMS} methods).
#'
setGeneric("getXcmsSet", function(obj, ...) standardGeneric("getXcmsSet"))

### MS Peak Lists

setGeneric("peakLists", function(obj) obj@peakLists)
setGeneric("generateMSPeakLists", function(fGroups, algorithm, ...) standardGeneric("generateMSPeakLists"))

### Components

setGeneric("componentTable", function(obj) standardGeneric("componentTable"))
setGeneric("componentInfo", function(obj) standardGeneric("componentInfo"))
setGeneric("findFGroup", function(obj, fGroup) standardGeneric("findFGroup"))
setGeneric("generateComponents", function(fGroups, algorithm, ...) standardGeneric("generateComponents"))

### Formulas

setGeneric("generateFormulas", function(fGroups, algorithm, ...) standardGeneric("generateFormulas"))

### Compounds

setGeneric("compoundTable", function(obj) standardGeneric("compoundTable"))
setGeneric("generateCompounds", function(fGroups, MSPeakLists, algorithm, ...) standardGeneric("generateCompounds"))
setGeneric("mergedCompoundNames", function(compounds) standardGeneric("mergedCompoundNames"))
setGeneric("identifiers", function(compounds) standardGeneric("identifiers"))
setGeneric("addFormulaScoring", function(compounds, formConsensus, updateScore = FALSE,
                                         formulaScoreWeight = 1) standardGeneric("addFormulaScoring"))
setGeneric("plotStructure", function(compounds, index, groupName, width = 500, height = 500, useGGPlot2 = FALSE) standardGeneric("plotStructure"))
setGeneric("plotScores", function(obj, index, groupName, normalizeScores = TRUE, useGGPlot2 = FALSE) standardGeneric("plotScores"))

### h-clustering

setGeneric("clusterProperties", function(cInfo) standardGeneric("clusterProperties"))
setGeneric("makeHCluster", function(fGroups, normFunc = max, metric = "euclidean", method = "ward.D2",
                                    average = TRUE) standardGeneric("makeHCluster"))
setGeneric("drawHeatMap", function(cInfo, col = colorRampPalette(c("black", "yellow"))(100),
                                   interactive = FALSE, ...) standardGeneric("drawHeatMap"))
setGeneric("getSilhouetteInfo", function(cInfo, ranges = seq(10, 100, 10)) standardGeneric("getSilhouetteInfo"))
setGeneric("hClusterFilter", function(fGroups, cInfo, k, c) standardGeneric("hClusterFilter"))

### target screening

#' Target and suspect screening
#'
#' Utilities to screen for targets with known or suspected identity.
#'
#' Besides 'full non-target analysis', where compounds may be identitified with
#' little to no prior knowledge, a common strategy is to screen for compounds
#' with known or suspected identity. This may be a generally favourable approach
#' if possible, as it can significantly reduce the load on data interpretation.
#'
#' @note Both \code{groupFeaturesScreening} and
#'   \code{importFeatureGroupsBrukerTASQ} use names from targets as feature
#'   group names, therefore, it is important that these are file-compatible
#'   names when \link[=reporting]{reporting data}.
#'
#' @name target-screening
NULL

#' @details \code{screenTargets} will screen a set of targets (provided as a
#'   \code{data.frame}) within an \code{\link{features}} or
#'   \code{\link{featureGroups}} object.
#'
#' @param obj The object that should be screened (\emph{i.e.}
#'   \code{\link{features}} or \code{\link{featureGroups}}).
#' @param targets A \code{data.frame} consisting of mandatory columns
#'   \code{"name"} (the targeted analyte name) and \code{"mz"} (the \emph{m/z}
#'   value of the target). In addition, a column \code{"ret"} can be used to
#'   specify the retention time (if unspecified no retention times are checked).
#' @param rtWindow,mzWindow The retention time window (in seconds) and
#'   \emph{m/z} window that will be used for matching a target.
#'
#' @return \code{screenTargets} will return a table (a \code{\link{data.table}})
#'   with detected targets and details such as retention and \emph{m/z}
#'   values. If a target is matched on multiple features/feature groups then
#'   each hit is reported as a separate row. In case concentration values were
#'   specified in the \link[=analysis-information]{analysis info} then the
#'   \code{featureGroups} method will also report correlation coefficients
#'   (\code{'RSQ'} column) calculated from the linear regression between
#'   intensity values and specified concentrations.
#'
#' @rdname target-screening
setGeneric("screenTargets", function(obj, targets, rtWindow = 12, mzWindow = 0.005) standardGeneric("screenTargets"))


### Misc.

#' Miscellaneous generics
#'
#' Various (S4) generic functions provide a common interface for common tasks
#' such as plotting and filtering data. The actual functionality and function
#' arguments are often specific for the implemented methods, for this reason,
#' please refer to the linked method documentation for each generic.
#'
#' @param obj The object the generic should be applied to.
#' @param \dots Any further method specific arguments. See method documentation
#'   for details.
#'
#' @name generics
NULL


#' @templateVar func analysisInfo
#' @templateVar desc returns the \link[=analysis-information]{analysis information} from an object.
#' @template generics
setGeneric("analysisInfo", function(obj) standardGeneric("analysisInfo"))

#' @templateVar func featureTable
#' @templateVar desc returns feature information.
#' @template generics
setGeneric("featureTable", function(obj) standardGeneric("featureTable"))

#' @templateVar func algorithm
#' @templateVar desc returns the algorithm that was used to generate an object.
#' @template generics
setGeneric("algorithm", function(obj) standardGeneric("algorithm"))

#' @templateVar func filter
#' @templateVar desc provides various functionality to do post-filtering of data.
#' @template generics
setGeneric("filter", function(obj, ...) standardGeneric("filter"))

#' @templateVar func plotInt
#' @templateVar desc plots the intensity of all contained features.
#' @template generics
setGeneric("plotInt", function(obj, ...) standardGeneric("plotInt"))

#' @templateVar func plotVenn
#' @templateVar desc plots a Venn diagram to assess unique and overlapping data.
#' @template generics
#'
#' @param which What should be plotted. See method documentation for specifics.
#'
setGeneric("plotVenn", function(obj, which = NULL, ...) standardGeneric("plotVenn"))

#' @templateVar func plotChord
#' @templateVar desc plots a Chord diagram to assess overlapping data.
#' @template generics
#'
#' @template plotChord-args
#'
setGeneric("plotChord", function(obj, addSelfLinks = FALSE, addRetMzPlots = TRUE, ...) standardGeneric("plotChord"))

#' @templateVar func plotSpec
#' @templateVar desc plots a (annotated) spectrum.
#' @template generics
setGeneric("plotSpec", function(obj, ...) standardGeneric("plotSpec"))

#' @templateVar func plotEIC
#' @templateVar desc plots extracted ion chromatogram(s).
#' @template generics
setGeneric("plotEIC", function(obj, ...) standardGeneric("plotEIC"))

#' @templateVar func consensus
#' @templateVar desc combines and merges data from various algorithms to
#'   generate a consensus.
#' @template generics
setGeneric("consensus", function(obj, ...) standardGeneric("consensus"))

#' @templateVar func formulaTable
#' @templateVar desc returns formula data.
#' @template generics
setGeneric("formulaTable", function(obj) standardGeneric("formulaTable"))

setGeneric("checkChromatograms", function(fGroups, mzWindow = 0.005, enabledFGroups = NULL) standardGeneric("checkChromatograms"))
setGeneric("compoundViewer", function(fGroups, MSPeakLists, compounds) standardGeneric("compoundViewer"))
setGeneric("reportCSV", function(fGroups, path = "report", reportFGroupsAsRows = TRUE, reportFGroupsAnalysisInfo = TRUE,
                                 reportFGroupsRetMz = TRUE, reportFeatures = FALSE, formConsensus = NULL,
                                 compounds = NULL, compoundNormalizeScores = TRUE, components = NULL,
                                 cInfo = NULL, clusterK = NULL, retMin = TRUE, clearPath = FALSE) standardGeneric("reportCSV"))
setGeneric("reportPDF", function(fGroups, path = "report", reportFGroups = TRUE,
                                 formConsensus = NULL, reportFormulaSpectra = TRUE, compounds = NULL, compoundNormalizeScores = TRUE,
                                 components = NULL, cInfo = NULL, clusterK = NULL, silInfo = NULL, clusterMaxLabels = 250,
                                 MSPeakLists = NULL, retMin = TRUE, EICGrid = c(2, 1), EICRtWindow = 20, EICMzWindow = 0.005,
                                 EICTopMost = NULL, EICOnlyPresent = TRUE, clearPath = FALSE) standardGeneric("reportPDF"))
setGeneric("reportMD", function(fGroups, path = "report", reportChord = TRUE, reportFGroups = TRUE,
                                formConsensus = NULL, reportFormulaSpectra = TRUE,
                                compounds = NULL, compoundNormalizeScores = TRUE,
                                components = NULL, cInfo = NULL, clusterK = NULL, silInfo = NULL,
                                interactiveHeat = FALSE, clusterMaxLabels = 250, MSPeakLists = NULL, retMin = TRUE,
                                EICRtWindow = 20, EICMzWindow = 0.005, EICTopMost = NULL, EICOnlyPresent = TRUE,
                                selfContained = TRUE, optimizePng = FALSE, maxProcAmount = getOption("patRoon.maxProcAmount"),
                                clearPath = FALSE) standardGeneric("reportMD"))
