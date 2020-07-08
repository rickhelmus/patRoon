if (!isGeneric("plot"))
    setGeneric("plot", function(x, y, ...) standardGeneric("plot"))


### Feature groups

setGeneric("groupTable", function(object, ...) standardGeneric("groupTable"))
setGeneric("groupFeatIndex", function(fGroups) standardGeneric("groupFeatIndex"))
setGeneric("groupInfo", function(fGroups) standardGeneric("groupInfo"))
setGeneric("removeEmptyAnalyses", function(fGroups) standardGeneric("removeEmptyAnalyses"))
setGeneric("averageGroups", function(fGroups, areas = FALSE) standardGeneric("averageGroups"))
setGeneric("export", function(fGroups, type, out) standardGeneric("export"))
setGeneric("unique", function(x, incomparables = FALSE, ...) standardGeneric("unique"))
setGeneric("overlap", function(fGroups, which, exclusive = FALSE) standardGeneric("overlap"))
setGeneric("comparison", function(..., groupAlgo,
                                  groupArgs = list(rtalign = FALSE)) standardGeneric("comparison"), signature = "...")
setGeneric("groupFeatures", function(feat, algorithm, ...) standardGeneric("groupFeatures"))
setGeneric("replicateGroupSubtract", function(fGroups, rGroups, threshold = 0) standardGeneric("replicateGroupSubtract"))
setGeneric("groupQualities", function(fGroups) standardGeneric("groupQualities"))
setGeneric("groupScores", function(fGroups) standardGeneric("groupScores"))

setGeneric("getPICSet", function(obj, ...) standardGeneric("getPICSet"))

### utils (XCMS)

#' Conversion to XCMS objects
#'
#' Converts a \code{\link{features}} or \code{\link{featureGroups}} object to an
#' \code{\link{xcmsSet}} or \code{\link{XCMSnExp}} object.
#'
#' @param obj The object that should be converted.
#' @param exportedData,\dots Set to \code{TRUE} if analyses were exported as
#'   \code{mzXML} or \code{mzML} files (ignored for \code{featuresOpenMS} and
#'   \code{featuresXCMS} methods).
#' @param verbose If \code{FALSE} then no text output is shown.
#' @rdname xcms-conv
setGeneric("getXCMSSet", function(obj, verbose = TRUE, ...) standardGeneric("getXCMSSet"))

#' @rdname xcms-conv
setGeneric("getXCMSnExp", function(obj, verbose = TRUE, ...) standardGeneric("getXCMSnExp"))

### MS Peak Lists

setGeneric("peakLists", function(obj) standardGeneric("peakLists"))
setGeneric("averagedPeakLists", function(obj) standardGeneric("averagedPeakLists"))
setGeneric("generateMSPeakLists", function(fGroups, algorithm, ...) standardGeneric("generateMSPeakLists"))
setGeneric("generateMSPeakListsMzR", function(fGroups, algorithm, ...) standardGeneric("generateMSPeakListsMzR"))

### Components

setGeneric("componentTable", function(obj) standardGeneric("componentTable"))
setGeneric("componentInfo", function(obj) standardGeneric("componentInfo"))
setGeneric("findFGroup", function(obj, fGroup) standardGeneric("findFGroup"))
setGeneric("plotGraph", function(obj, onlyLinked = TRUE) standardGeneric("plotGraph"))
setGeneric("generateComponents", function(fGroups, algorithm, ...) standardGeneric("generateComponents"))

### Formulas

setGeneric("formulaTable", function(obj, features = FALSE) standardGeneric("formulaTable"))
setGeneric("generateFormulas", function(fGroups, algorithm, ...) standardGeneric("generateFormulas"))

### Compounds

setGeneric("compoundTable", function(obj) standardGeneric("compoundTable"))
setGeneric("generateCompounds", function(fGroups, MSPeakLists, algorithm, ...) standardGeneric("generateCompounds"))
setGeneric("mergedCompoundNames", function(compounds) standardGeneric("mergedCompoundNames"))
setGeneric("identifiers", function(compounds) standardGeneric("identifiers"))
setGeneric("addFormulaScoring", function(compounds, formulas, updateScore = FALSE,
                                         formulaScoreWeight = 1) standardGeneric("addFormulaScoring"))
setGeneric("settings", function(compoundsMF) standardGeneric("settings"))

### clustering

setGeneric("makeHCluster", function(obj, method = "complete", ...) standardGeneric("makeHCluster"))

setGeneric("plotHeatMap", function(obj, ...) standardGeneric("plotHeatMap"))

### Sets

setGeneric("sets", function(obj) standardGeneric("sets"))
setGeneric("adducts", function(obj) standardGeneric("adducts"))
setGeneric("ionize", function(obj, sets = NULL) standardGeneric("ionize"))
setGeneric("groupAlgorithm", function(obj) standardGeneric("groupAlgorithm"))

### suspect screening

setGeneric("screenSuspects", function(fGroups, suspects, rtWindow = 12, mzWindow = 0.005, adduct = NULL,
                                      skipInvalid = TRUE, onlyHits = FALSE) standardGeneric("screenSuspects"))
setGeneric("screenInfo", function(obj) standardGeneric("screenInfo"))
setGeneric("annotateSuspects", function(fGroups, MSPeakLists = NULL, formulas = NULL,
                                        compounds = NULL, ...) standardGeneric("annotateSuspects"))


### Optimization

setGeneric("optimizedParameters", function(object, paramSet = NULL, DoEIteration = NULL) standardGeneric("optimizedParameters"))
setGeneric("optimizedObject", function(object, paramSet = NULL) standardGeneric("optimizedObject"))
setGeneric("scores", function(object, paramSet = NULL, DoEIteration = NULL) standardGeneric("scores"))
setGeneric("experimentInfo", function(object, paramSet, DoEIteration) standardGeneric("experimentInfo"))


### Misc.

#' Miscellaneous generics
#'
#' Various (S4) generic functions providing a common interface for common tasks
#' such as plotting and filtering data. The actual functionality and function
#' arguments are often specific for the implemented methods, for this reason,
#' please refer to the linked method documentation for each generic.
#'
#' @param obj The object the generic should be applied to.
#' @param \dots Any further method specific arguments. See method documentation
#'   for details.
#'
#' @section Other generics: Below are methods that are defined for existing
#'   generics (\emph{e.g.} defined in \code{base}). Please see method specific
#'   documentation for more details.
#'
#' @name generics
NULL


#' @templateVar func algorithm
#' @templateVar desc returns the algorithm that was used to generate the object.
#' @template generics
setGeneric("algorithm", function(obj) standardGeneric("algorithm"))

#' @templateVar func analysisInfo
#' @templateVar desc returns the \link[=analysis-information]{analysis information} from an object.
#' @template generics
setGeneric("analysisInfo", function(obj) standardGeneric("analysisInfo"))

#' @templateVar func analyses
#' @templateVar desc returns a \code{character} vector with the analyses for which data is present in this object.
#' @template generics
setGeneric("analyses", function(obj) standardGeneric("analyses"))

#' @templateVar func annotatedPeakList
#' @templateVar desc returns an annotated MS peak list.
#' @template generics
setGeneric("annotatedPeakList", function(obj, ...) standardGeneric("annotatedPeakList"))

#' @templateVar func calculatePeakQualities
#' @templateVar desc calculates chromatographic peak qualities and scores.
#' @template generics
setGeneric("calculatePeakQualities", function(obj, weights = NULL,
                                              flatnessFactor = 0.05, ...) standardGeneric("calculatePeakQualities"))

#' @templateVar func clusterProperties
#' @templateVar desc Obtain a list with properties of the generated cluster(s).
#' @template generics
setGeneric("clusterProperties", function(obj) standardGeneric("clusterProperties"))

#' @templateVar func clusters
#' @templateVar desc Obtain clustering object(s).
#' @template generics
setGeneric("clusters", function(obj) standardGeneric("clusters"))

#' @templateVar func cutClusters
#' @templateVar desc Returns assigned cluster indices of a cut cluster.
#' @template generics
setGeneric("cutClusters", function(obj) standardGeneric("cutClusters"))

#' @templateVar func consensus
#' @templateVar desc combines and merges data from various algorithms to generate a consensus.
#' @template generics
setGeneric("consensus", function(obj, ...) standardGeneric("consensus"))

#' @templateVar func export
#' @templateVar desc exports workflow data to a given format.
#' @template generics
setGeneric("export", function(obj, type, out, ...) standardGeneric("export"))

#' @templateVar func featureTable
#' @templateVar desc returns feature information.
#' @template generics
setGeneric("featureTable", function(obj, ...) standardGeneric("featureTable"))

#' @templateVar func featureTable
#' @templateVar desc updates feature information.
#' @param fTable New value for the \code{features} slot: a \code{list} with \code{data.table} items for each analysis.
#' @template generics
setGeneric("featureTable<-", function(obj, value) standardGeneric("featureTable<-"))

#' @templateVar func filter
#' @templateVar desc provides various functionality to do post-filtering of data.
#' @template generics
setGeneric("filter", function(obj, ...) standardGeneric("filter"))

#' @templateVar func getFeatures
#' @templateVar desc returns the object's \code{\link{features}} object.
#' @template generics
setGeneric("getFeatures", function(obj) standardGeneric("getFeatures"))

#' @templateVar func getMCS
#' @templateVar desc Calculcates the maximum common substructure.
#' @template generics
setGeneric("getMCS", function(obj, ...) standardGeneric("getMCS"))

#' @templateVar func groupNames
#' @templateVar desc returns a \code{character} vector with the names of the feature groups for which data is present in this object.
#' @template generics
setGeneric("groupNames", function(obj) standardGeneric("groupNames"))

#' @templateVar func plotChord
#' @templateVar desc plots a Chord diagram to assess overlapping data.
#' @template generics
#'
#' @template plotChord-args
#'
setGeneric("plotChord", function(obj, addSelfLinks = FALSE, addRetMzPlots = TRUE, ...) standardGeneric("plotChord"))

#' @templateVar func plotChroms
#' @templateVar desc plots extracted ion chromatogram(s).
#' @template generics
setGeneric("plotChroms", function(obj, ...) standardGeneric("plotChroms"))

#' @templateVar func plotInt
#' @templateVar desc plots the intensity of all contained features.
#' @template generics
setGeneric("plotInt", function(obj, ...) standardGeneric("plotInt"))

#' @templateVar func plotScores
#' @templateVar desc plots candidate scorings.
#' @template generics
setGeneric("plotScores", function(obj, ...) standardGeneric("plotScores"))

#' @templateVar func plotSilhouettes
#' @templateVar desc plots silhouette widths to evaluate the desired cluster size.
#' @template generics
#' @param kSeq An integer vector containing the sequence that should be used for
#'   average silhouette width calculation.
setGeneric("plotSilhouettes", function(obj, kSeq, ...) standardGeneric("plotSilhouettes"))

#' @templateVar func plotSpectrum
#' @templateVar desc plots a (annotated) spectrum.
#' @template generics
setGeneric("plotSpectrum", function(obj, ...) standardGeneric("plotSpectrum"))

#' @templateVar func plotStructure
#' @templateVar desc plots a chemical structure.
#' @template generics
setGeneric("plotStructure", function(obj, ...) standardGeneric("plotStructure"))

#' @templateVar func plotVenn
#' @templateVar desc plots a Venn diagram to assess unique and overlapping data.
#' @template generics
#'
setGeneric("plotVenn", function(obj, ...) standardGeneric("plotVenn"))

#' @templateVar func plotUpSet
#' @templateVar desc plots an UpSet diagram to assess unique and overlapping data.
#' @template generics
setGeneric("plotUpSet", function(obj, ...) standardGeneric("plotUpSet"))

#' @templateVar func delete
#' @templateVar desc Deletes results.
#' @template generics
setGeneric("delete", function(obj, ...) standardGeneric("delete"))

#' @templateVar func replicateGroups
#' @templateVar desc returns a \code{character} vector with the analyses for which data is present in this object.
#' @template generics
setGeneric("replicateGroups", function(obj) standardGeneric("replicateGroups"))

#' @templateVar func treeCut
#' @templateVar desc Manually cut a cluster.
#' @template generics
#' @param k,h Desired numbers of clusters. See \code{\link{cutree}}.
setGeneric("treeCut", function(obj, k = NULL, h = NULL, ...) standardGeneric("treeCut"))

#' @templateVar func treeCutDynamic
#' @templateVar desc Automatically cut a cluster.
#' @template generics
#' @template dynamictreecut
setGeneric("treeCutDynamic", function(obj, maxTreeHeight = 1, deepSplit = TRUE,
                                      minModuleSize = 1, ...) standardGeneric("treeCutDynamic"))


setGeneric("checkFeatures", function(fGroups, session, rtWindow = 30,
                                     clearSession = FALSE) standardGeneric("checkFeatures"))
setGeneric("getEICsForFGroups", function(fGroups, rtWindow, mzExpWindow, topMost, topMostByRGroup,
                                         onlyPresent) standardGeneric("getEICsForFGroups"))
setGeneric("compoundViewer", function(fGroups, MSPeakLists, compounds) standardGeneric("compoundViewer"))
setGeneric("reportCSV", function(fGroups, path = "report", reportFeatures = FALSE, formulas = NULL,
                                 formulasNormalizeScores = "max", formulasExclNormScores = NULL,
                                 compounds = NULL, compoundsNormalizeScores = "max",
                                 compoundsExclNormScores = c("score", "individualMoNAScore", "annoTypeCount"),
                                 compsCluster = NULL, components = NULL,
                                 retMin = TRUE, clearPath = FALSE) standardGeneric("reportCSV"))
setGeneric("reportPDF", function(fGroups, path = "report", reportFGroups = TRUE,
                                 formulas = NULL, formulasTopMost = 5,
                                 formulasNormalizeScores = "max", formulasExclNormScores = NULL,
                                 reportFormulaSpectra = TRUE,
                                 compounds = NULL, compoundsNormalizeScores = "max",
                                 compoundsExclNormScores = c("score", "individualMoNAScore", "annoTypeCount"),
                                 compoundsOnlyUsedScorings = TRUE, compoundsTopMost = 5,
                                 compsCluster = NULL, components = NULL, MSPeakLists = NULL, retMin = TRUE,
                                 EICGrid = c(2, 1), EICRtWindow = 20, EICMzExpWindow = 0.001, EICTopMost = NULL,
                                 EICTopMostByRGroup = FALSE, EICOnlyPresent = TRUE,
                                 clearPath = FALSE) standardGeneric("reportPDF"))
setGeneric("reportHTML", function(fGroups, path = "report", reportPlots = c("chord", "venn", "upset", "eics", "formulas"),
                                  formulas = NULL, formulasTopMost = 5,
                                  formulasNormalizeScores = "max", formulasExclNormScores = NULL,
                                  compounds = NULL, compoundsNormalizeScores = "max",
                                  compoundsExclNormScores = c("score", "individualMoNAScore", "annoTypeCount"),
                                  compoundsOnlyUsedScorings = TRUE, compoundsTopMost = 5, compsCluster = NULL,
                                  includeMFWebLinks = "compounds", components = NULL, interactiveHeat = FALSE,
                                  MSPeakLists = NULL, retMin = TRUE, EICRtWindow = 20, EICMzExpWindow = 0.001,
                                  EICTopMost = NULL, EICTopMostByRGroup = FALSE, EICOnlyPresent = TRUE,
                                  selfContained = TRUE, optimizePng = FALSE, clearPath = FALSE, openReport = TRUE,
                                  noDate = FALSE) standardGeneric("reportHTML"))

# Used for reporting
setGeneric("plotHash", function(x, ...) standardGeneric("plotHash"))
setGeneric("plotSpectrumHash", function(obj, ...) standardGeneric("plotSpectrumHash"))
setGeneric("plotScoresHash", function(obj, ...) standardGeneric("plotScoresHash"))
setGeneric("plotStructureHash", function(obj, ...) standardGeneric("plotStructureHash"))
setGeneric("plotChromsHash", function(obj, ...) standardGeneric("plotChromsHash"))
setGeneric("plotIntHash", function(obj, ...) standardGeneric("plotIntHash"))


#' @templateVar func [
#' @templateVar desc Subsets data within an object.
#' @template existing-generics
NULL

#' @templateVar func [[
#' @templateVar desc Extract data from an object.
#' @template existing-generics
NULL

#' @templateVar func $
#' @templateVar desc Extract data from an object.
#' @template existing-generics
NULL

#' @templateVar func as.data.table
#' @templateVar desc Converts an object to a table (\code{\link{data.table}}).
#' @template existing-generics
NULL

#' @templateVar func as.data.frame
#' @templateVar desc Converts an object to a table (\code{data.frame}).
#' @template existing-generics
NULL

#' @templateVar func length
#' @templateVar desc Returns the length of an object.
#' @template existing-generics
NULL

#' @templateVar func lengths
#' @templateVar desc Returns the lengths of elements within this object.
#' @template existing-generics
NULL

#' @templateVar func names
#' @templateVar desc Return names for this object.
#' @template existing-generics
NULL

#' @templateVar func plot
#' @templateVar desc Generates a plot for an object.
#' @template existing-generics
NULL

#' @templateVar func show
#' @templateVar desc Prints information about this object.
#' @template existing-generics
NULL
