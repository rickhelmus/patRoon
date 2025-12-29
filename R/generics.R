# SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @name feature-plotting
NULL

# NOTE/HACK: this file sets a lot of @name instead of @rdname roxygen tags, as it is parsed before other files and the
# @name tag is otherwise ignored.


### Non-exported generics

setGeneric("doGroupFeatures", function(feat, grouper, ..., verbose) standardGeneric("doGroupFeatures"))
setGeneric("groupFeaturesIMS", function(feat, ...) standardGeneric("groupFeaturesIMS"))
setGeneric("featureTable<-", function(obj, value) standardGeneric("featureTable<-"))
setGeneric("reorderAnalyses", function(obj, anas) standardGeneric("reorderAnalyses"))
setGeneric("removeEmptyAnalyses", function(fGroups) standardGeneric("removeEmptyAnalyses"))
setGeneric("averageGroups", function(fGroups, areas = FALSE, normalized = FALSE, by = "group",
                                     func = mean) standardGeneric("averageGroups"))
setGeneric("clearMobilities", function(obj) standardGeneric("clearMobilities"))
setGeneric("averageMSPeakLists", function(obj) standardGeneric("averageMSPeakLists"))
setGeneric("collapseComponents", function(obj) standardGeneric("collapseComponents"))
setGeneric("annScoreNames", function(obj, onlyNums) standardGeneric("annScoreNames"))
setGeneric("getFragInfo", function(obj, groupName, index, ...) standardGeneric("getFragInfo"))
setGeneric("prepareConsensusLabels", function(obj, ..., labels) standardGeneric("prepareConsensusLabels"))
setGeneric("mergedConsensusNames", function(obj, sets = TRUE) standardGeneric("mergedConsensusNames"))
setGeneric("identifiers", function(compounds) standardGeneric("identifiers"))
setGeneric("updateSetConsensus", function(obj) standardGeneric("updateSetConsensus"))
setGeneric("parentsFromScreening", function(obj) standardGeneric("parentsFromScreening"))
setGeneric("TPsFromScreening", function(obj) standardGeneric("TPsFromScreening"))
setGeneric("linkParentsToFGroups", function(TPs, fGroups) standardGeneric("linkParentsToFGroups"))
setGeneric("linkTPsToFGroups", function(TPs, fGroups) standardGeneric("linkTPsToFGroups"))
setGeneric("groupNamesResults", function(obj) standardGeneric("groupNamesResults"))

# Used for reporting
setGeneric("plotHash", function(x, ...) standardGeneric("plotHash"))
setGeneric("plotSpectrumHash", function(obj, ...) standardGeneric("plotSpectrumHash"))
setGeneric("plotScoresHash", function(obj, ...) standardGeneric("plotScoresHash"))
setGeneric("plotStructureHash", function(obj, ...) standardGeneric("plotStructureHash"))
setGeneric("plotChromsHash", function(obj, ...) standardGeneric("plotChromsHash"))
setGeneric("plotMobilogramsHash", function(obj, ...) standardGeneric("plotMobilogramsHash"))
setGeneric("plotIntHash", function(obj, ...) standardGeneric("plotIntHash"))
setGeneric("plotChordHash", function(obj, ...) standardGeneric("plotChordHash"))
setGeneric("plotVennHash", function(obj, ...) standardGeneric("plotVennHash"))
setGeneric("plotUpSetHash", function(obj, ...) standardGeneric("plotUpSetHash"))

setGeneric("getFeatureEIXInputTab", function(obj, type, ...) standardGeneric("getFeatureEIXInputTab"))
setGeneric("getFeatureEIXs", function(obj, type, ...) standardGeneric("getFeatureEIXs"))

setGeneric("initMSReadBackend", function(backend) standardGeneric("initMSReadBackend"))

### Features and feature groups

#' @name featureGroups-class
setGeneric("groupTable", function(object, ...) standardGeneric("groupTable"))

#' @name featureGroups-class
setGeneric("groupFeatIndex", function(fGroups) standardGeneric("groupFeatIndex"))

#' @name featureGroups-class
setGeneric("groupInfo", function(fGroups) standardGeneric("groupInfo"))

#' @name featureGroups-class
setGeneric("unique", function(x, incomparables = FALSE, ...) standardGeneric("unique"))

#' @name featureGroups-class
setGeneric("overlap", function(fGroups, which = NULL, aggregate = TRUE, exclusive = FALSE, ...) standardGeneric("overlap"))

#' @name featureGroups-class
setGeneric("selectIons", function(fGroups, components, prefAdduct, ...) standardGeneric("selectIons"))

#' @name featureGroups-compare
setGeneric("comparison", function(..., groupAlgo,
                                  groupArgs = list(rtalign = FALSE)) standardGeneric("comparison"), signature = "...")

#' @name feature-filtering
setGeneric("replicateSubtract", function(fGroups, replicates, threshold = 0) standardGeneric("replicateSubtract"))

#' @name featureGroups-class
setGeneric("groupQualities", function(fGroups) standardGeneric("groupQualities"))

#' @name featureGroups-class
setGeneric("groupScores", function(fGroups) standardGeneric("groupScores"))

#' @name featureGroups-class
setGeneric("internalStandards", function(fGroups) standardGeneric("internalStandards"))

#' @name featureGroups-class
setGeneric("internalStandardAssignments", function(fGroups, ...) standardGeneric("internalStandardAssignments"))

#' @name featureGroups-class
setGeneric("normInts", function(fGroups, featNorm = "none", groupNorm = FALSE, normFunc = max, standards = NULL,
                                ISTDRTWindow = 120, ISTDMZWindow = 300, minISTDs = 3, ...) standardGeneric("normInts"))

#' @name featureGroups-class
setGeneric("concentrations", function(fGroups, ...) standardGeneric("concentrations"))

#' @name featureGroups-class
setGeneric("calculateConcs", function(fGroups, ...) standardGeneric("calculateConcs"))

#' @name featureGroups-class
setGeneric("toxicities", function(fGroups, ...) standardGeneric("toxicities"))

#' @name featureGroups-class
setGeneric("calculateTox", function(fGroups, ...) standardGeneric("calculateTox"))

#' @name featureGroups-class
setGeneric("updateGroups", function(fGroups, ...) standardGeneric("updateGroups"))

#' @name feature-plotting
setGeneric("plotMobilograms", function(obj, ...) standardGeneric("plotMobilograms"))

#' @name kpic2-conv
setGeneric("getPICSet", function(obj, ...) standardGeneric("getPICSet"))

#' @name groupFeatures
setGeneric("groupFeatures", function(obj, algorithm, ...) standardGeneric("groupFeatures"))

#' @rdname groupFeaturesOpenMS
setGeneric("groupFeaturesOpenMS", function(feat, ...) standardGeneric("groupFeaturesOpenMS"))

#' @rdname groupFeaturesXCMS
setGeneric("groupFeaturesXCMS", function(feat, ...) standardGeneric("groupFeaturesXCMS"))

#' @rdname groupFeaturesXCMS3
setGeneric("groupFeaturesXCMS3", function(feat, ...) standardGeneric("groupFeaturesXCMS3"))

#' @rdname groupFeaturesKPIC2
setGeneric("groupFeaturesKPIC2", function(feat, ...) standardGeneric("groupFeaturesKPIC2"))

#' @rdname groupFeaturesGreedy
setGeneric("groupFeaturesGreedy", function(feat, ...) standardGeneric("groupFeaturesGreedy"))

#' @name check-GUI
setGeneric("checkFeatures", function(fGroups, session = "checked-features.yml", EICParams = getDefEICParams(),
                                     EIMParams = getDefEIMParams(),
                                     clearSession = FALSE) standardGeneric("checkFeatures"))


### utils (XCMS)

#' Conversion to XCMS objects
#'
#' Converts a \code{\link{features}} or \code{\link{featureGroups}} object to an \code{\link{xcmsSet}} or
#' \code{\link{XCMSnExp}} object.
#'
#' The conversion process will introduce some dummy values for metadata not present in \pkg{patRoon} objects. If the
#' \code{features} or \code{featureGroups} object was generated with \pkg{XCMS}, then no conversion is performed and the
#' original \pkg{XCMS} object will be returned, if possible. Conversion may still occur \emph{e.g.} due to the
#' application of some subsetting or filtering steps or the re-ordering of analyses.
#'
#' @param obj The object that should be converted.
#' @param \dots \setsWF Further arguments passed to non-sets method.
#'
#'   Otherwise ignored.
#' @param verbose If \code{FALSE} then no text output is shown.
#' @param set \setsWF The name of the set to be exported.
#'
#' @template loadrawdata-arg
#'
#' @templateVar consider for export
#' @templateVar append This should be kept \code{FALSE} as \pkg{XCMS} export currently does not support IMS features.
#' @template IMS-arg
#'
#' @section Sets workflows: In a \link[=sets-workflow]{sets workflow}, \code{\link{unset}} is used to convert the
#'   feature (group) data before the object is exported.
#'
#' @references reference \addCitations{xcms}
#'
#' @name xcms-conv
setGeneric("getXCMSSet", function(obj, verbose = TRUE, ...) standardGeneric("getXCMSSet"))

#' @rdname xcms-conv
setGeneric("getXCMSnExp", function(obj, verbose = TRUE, ...) standardGeneric("getXCMSnExp"))


### MS Peak Lists

#' @name MSPeakLists-class
setGeneric("peakLists", function(obj, ...) standardGeneric("peakLists"))

#' @name MSPeakLists-class
setGeneric("averagedPeakLists", function(obj, ...) standardGeneric("averagedPeakLists"))

#' @name MSPeakLists-class
setGeneric("spectrumSimilarity", function(obj, ...) standardGeneric("spectrumSimilarity"))

#' @name MSPeakLists-class
setGeneric("spectrumSimilarityMobility", function(obj, ...) standardGeneric("spectrumSimilarityMobility"))

#' @rdname generateMSPeakLists
setGeneric("generateMSPeakLists", function(fGroups, ...) standardGeneric("generateMSPeakLists"))

#' @rdname generateMSPeakListsMzR
setGeneric("generateMSPeakListsMzR", function(fGroups, ...) standardGeneric("generateMSPeakListsMzR"))

#' @rdname generateMSPeakListsDA
setGeneric("generateMSPeakListsDA", function(fGroups, ...) standardGeneric("generateMSPeakListsDA"))

#' @rdname generateMSPeakListsDAFMF
setGeneric("generateMSPeakListsDAFMF", function(fGroups, ...) standardGeneric("generateMSPeakListsDAFMF"))


### Components

#' @name components-class
setGeneric("componentTable", function(obj) standardGeneric("componentTable"))

#' @name components-class
setGeneric("componentInfo", function(obj) standardGeneric("componentInfo"))

#' @name components-class
setGeneric("findFGroup", function(obj, fGroup) standardGeneric("findFGroup"))

#' @name check-GUI
setGeneric("checkComponents", function(components, fGroups, session = "checked-components.yml",
                                       EICParams = getDefEICParams(),
                                       clearSession = FALSE) standardGeneric("checkComponents"))

#' @rdname generateComponents
setGeneric("generateComponents", function(fGroups, algorithm, ...) standardGeneric("generateComponents"))

#' @rdname generateComponentsRAMClustR
setGeneric("generateComponentsRAMClustR", function(fGroups, ...) standardGeneric("generateComponentsRAMClustR"))

#' @rdname generateComponentsCAMERA
setGeneric("generateComponentsCAMERA", function(fGroups, ...) standardGeneric("generateComponentsCAMERA"))

#' @rdname generateComponentsNontarget
setGeneric("generateComponentsNontarget", function(fGroups, ...) standardGeneric("generateComponentsNontarget"))

#' @rdname generateComponentsIntClust
setGeneric("generateComponentsIntClust", function(fGroups, ...) standardGeneric("generateComponentsIntClust"))

#' @rdname generateComponentsOpenMS
setGeneric("generateComponentsOpenMS", function(fGroups, ...) standardGeneric("generateComponentsOpenMS"))

#' @rdname generateComponentsCliqueMS
setGeneric("generateComponentsCliqueMS", function(fGroups, ...) standardGeneric("generateComponentsCliqueMS"))

#' @rdname generateComponentsSpecClust
setGeneric("generateComponentsSpecClust", function(fGroups, ...) standardGeneric("generateComponentsSpecClust"))

#' @param \dots Further arguments specified to the methods.
#' @rdname generateComponentsTPs
setGeneric("generateComponentsTPs", function(fGroups, ...) standardGeneric("generateComponentsTPs"))

#' @name components-class
setGeneric("expandForIMS", function(obj, ...) standardGeneric("expandForIMS"))


### Formulas

#' @rdname generateFormulas
setGeneric("generateFormulas", function(fGroups, MSPeakLists, algorithm,
                                        specSimParams = getDefSpecSimParams(removePrecursor = TRUE), ...) standardGeneric("generateFormulas"))

#' @rdname generateFormulasGenForm
setGeneric("generateFormulasGenForm", function(fGroups, ...) standardGeneric("generateFormulasGenForm"))

#' @rdname generateFormulasSIRIUS
setGeneric("generateFormulasSIRIUS", function(fGroups, ...) standardGeneric("generateFormulasSIRIUS"))

#' @rdname generateFormulasDA
setGeneric("generateFormulasDA", function(fGroups, ...) standardGeneric("generateFormulasDA"))

### Compounds

#' @name compounds-class
setGeneric("addFormulaScoring", function(compounds, formulas, updateScore = FALSE,
                                         formulaScoreWeight = 1) standardGeneric("addFormulaScoring"))

#' @name compoundsMF-class
setGeneric("settings", function(compoundsMF) standardGeneric("settings"))

#' @rdname generateCompounds
setGeneric("generateCompounds", function(fGroups, MSPeakLists, algorithm,
                                         specSimParams = getDefSpecSimParams(removePrecursor = TRUE), ...) standardGeneric("generateCompounds"))

#' @rdname generateCompoundsMetFrag
setGeneric("generateCompoundsMetFrag", function(fGroups, ...) standardGeneric("generateCompoundsMetFrag"))

#' @rdname generateCompoundsSIRIUS
setGeneric("generateCompoundsSIRIUS", function(fGroups, ...) standardGeneric("generateCompoundsSIRIUS"))

#' @rdname generateCompoundsLibrary
setGeneric("generateCompoundsLibrary", function(fGroups, ...) standardGeneric("generateCompoundsLibrary"))

#' @name patRoon-defunct
setGeneric("compoundViewer", function(fGroups, MSPeakLists, compounds) standardGeneric("compoundViewer"))

### clustering

#' @param \dots further arguments specified to methods.
#' @name compounds-cluster
setGeneric("makeHCluster", function(obj, method = "complete", ...) standardGeneric("makeHCluster"))

#' @name componentsIntClust-class
setGeneric("plotHeatMap", function(obj, ...) standardGeneric("plotHeatMap"))

### Sets

#' @name makeSet
setGeneric("makeSet", function(obj, ...) standardGeneric("makeSet"))

### TP generation

#' @name transformationProducts-class
setGeneric("parents", function(TPs) standardGeneric("parents"))

#' @name transformationProducts-class
setGeneric("products", function(TPs) standardGeneric("products"))

#' @param \dots Further arguments specified to the methods.
#' @rdname generateTPsLogic
setGeneric("generateTPsLogic", function(fGroups, minMass = 40, ...) standardGeneric("generateTPsLogic"))

### suspect screening

#' @name featureGroupsScreening-class
setGeneric("screenInfo", function(obj) standardGeneric("screenInfo"))

#' @param \dots Further arguments specified to the methods.
#' @name suspect-screening
setGeneric("screenSuspects", function(fGroups, suspects, rtWindow = defaultLim("retention", "medium"),
                                      mzWindow = defaultLim("mz", "medium"), IMSMatchParams = NULL,
                                      adduct = NULL, skipInvalid = TRUE, prefCalcChemProps = TRUE,
                                      neutralChemProps = FALSE,
                                      onlyHits = FALSE, ...) standardGeneric("screenSuspects"))

### Optimization

#' @name optimizationResult-class
setGeneric("optimizedParameters", function(object, paramSet = NULL, DoEIteration = NULL) standardGeneric("optimizedParameters"))

#' @name optimizationResult-class
setGeneric("optimizedObject", function(object, paramSet = NULL) standardGeneric("optimizedObject"))

#' @name optimizationResult-class
setGeneric("scores", function(object, paramSet = NULL, DoEIteration = NULL) standardGeneric("scores"))

#' @name optimizationResult-class
setGeneric("experimentInfo", function(object, paramSet, DoEIteration) standardGeneric("experimentInfo"))


### MS library

#' @name MSLibrary-class
setGeneric("records", function(obj) standardGeneric("records"))

#' @name MSLibrary-class
setGeneric("spectra", function(obj) standardGeneric("spectra"))

### Reporting

#' @name reporting-legacy
setGeneric("reportCSV", function(fGroups, path = "report", reportFeatures = FALSE, formulas = NULL,
                                 formulasNormalizeScores = "max", formulasExclNormScores = NULL,
                                 compounds = NULL, compoundsNormalizeScores = "max",
                                 compoundsExclNormScores = c("score", "individualMoNAScore", "annoTypeCount",
                                                             "annotHitCount", "libMatch"),
                                 compsCluster = NULL, components = NULL,
                                 retMin = TRUE, clearPath = FALSE) standardGeneric("reportCSV"))
#' @name reporting-legacy
setGeneric("reportPDF", function(fGroups, path = "report", reportFGroups = TRUE,
                                 formulas = NULL, formulasTopMost = 5,
                                 formulasNormalizeScores = "max", formulasExclNormScores = NULL,
                                 reportFormulaSpectra = TRUE,
                                 compounds = NULL, compoundsNormalizeScores = "max",
                                 compoundsExclNormScores = c("score", "individualMoNAScore", "annoTypeCount",
                                                             "annotHitCount", "libMatch"),
                                 compoundsOnlyUsedScorings = TRUE, compoundsTopMost = 5,
                                 compsCluster = NULL, components = NULL, MSPeakLists = NULL, retMin = TRUE,
                                 EICGrid = c(2, 1), EICParams = getDefEICParams(window = 20, topMost = 1,
                                                                                topMostByReplicate = TRUE),
                                 clearPath = FALSE) standardGeneric("reportPDF"))
#' @name reporting
setGeneric("report", function(fGroups, MSPeakLists = NULL, formulas = NULL, compounds = NULL, compsCluster = NULL,
                              components = NULL, TPs = NULL,
                              settingsFile = system.file("report", "settings.yml", package = "patRoon"),
                              path = NULL, EICParams = getDefEICParams(topMost = 1, topMostByReplicate = TRUE),
                              EIMParams = getDefEIMParams(topMost = 1, topMostByReplicate = TRUE),
                              specSimParams = getDefSpecSimParams(), clearPath = FALSE, openReport = TRUE,
                              parallel = TRUE, overrideSettings = list()) standardGeneric("report"))


### Misc.

#' Miscellaneous generics
#'
#' Various (S4) generic functions providing a common interface for common tasks
#' such as plotting and filtering data. The actual functionality and function
#' arguments are often specific for the implemented methods, for this reason,
#' please refer to the linked method documentation for each generic.
#'
#' @param obj The object the generic should be applied to.
#' @param TPs The \code{\link{transformationProducts}} derived object.
#' @param out Output file.
#' @param value The replacement value.
#' @param \dots Any further method specific arguments. See method documentation
#'   for details.
#'
#' @section Other generics: Below are methods that are defined for existing
#'   generics (\emph{e.g.} defined in \code{base}). Please see method specific
#'   documentation for more details.
#'
#' @name generics
NULL

#' @param \dots Method specific arguments.
#' @name id-conf
setGeneric("estimateIDConfidence", function(obj, ...) standardGeneric("estimateIDConfidence"))


#' @templateVar func adducts
#' @templateVar desc returns assigned adducts of the object.
#' @template generics
setGeneric("adducts", function(obj, ...) standardGeneric("adducts"))

#' @templateVar func adducts<-
#' @templateVar desc sets adducts of the object.
#' @template generics
setGeneric("adducts<-", function(obj, value, ...) standardGeneric("adducts<-"))

#' @templateVar func algorithm
#' @templateVar desc returns the algorithm that was used to generate the object.
#' @template generics
setGeneric("algorithm", function(obj) standardGeneric("algorithm"))

#' @templateVar func analysisInfo
#' @templateVar desc returns the \link[=analysis-information]{analysis information} of an object.
#' @template generics
#' @param df If \code{TRUE} then a \code{data.frame} is returned, otherwise a \code{data.table} is returned.
setGeneric("analysisInfo", function(obj, df = FALSE) standardGeneric("analysisInfo"))

#' @templateVar func analysisInfo<-
#' @templateVar desc modifies the \link[=analysis-information]{analysis information} of an object.
#' @template generics
setGeneric("analysisInfo<-", function(obj, value) standardGeneric("analysisInfo<-"))

#' @templateVar func analyses
#' @templateVar desc returns a \code{character} vector with the analyses for which data is present in this object.
#' @template generics
setGeneric("analyses", function(obj) standardGeneric("analyses"))

#' @templateVar func annotatedPeakList
#' @templateVar desc returns an annotated MS peak list.
#' @template generics
setGeneric("annotatedPeakList", function(obj, ...) standardGeneric("annotatedPeakList"))

#' @templateVar func annotations
#' @templateVar desc returns annotations.
#' @template generics
setGeneric("annotations", function(obj, ...) standardGeneric("annotations"))

#' @templateVar func assignMobilities
#' @templateVar desc assigns ion mobility and/or \acronym{CCS} values to workflow data.
#' @template generics
setGeneric("assignMobilities", function(obj, ...) standardGeneric("assignMobilities"))

#' @templateVar func calculatePeakQualities
#' @templateVar desc calculates chromatographic peak qualities and scores.
#' @template generics
#' @param weights,flatnessFactor See method documentation.
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

#' @templateVar func consensus
#' @templateVar desc combines and merges data from various algorithms to generate a consensus.
#' @template generics
setGeneric("consensus", function(obj, ...) standardGeneric("consensus"))

#' @templateVar func convertToMFDB
#' @templateVar desc Exports the object to a local database that can be used with \command{MetFrag}.
#' @template generics
setGeneric("convertToMFDB", function(TPs, out, ...) standardGeneric("convertToMFDB"))

#' @templateVar func convertToSuspects
#' @templateVar desc Converts an object to a suspect list.
#' @template generics
setGeneric("convertToSuspects", function(obj, ...) standardGeneric("convertToSuspects"))

#' @templateVar func cutClusters
#' @templateVar desc Returns assigned cluster indices of a cut cluster.
#' @template generics
setGeneric("cutClusters", function(obj) standardGeneric("cutClusters"))

#' @templateVar func defaultExclNormScores
#' @templateVar desc Returns default scorings that are excluded from normalization.
#' @template generics
setGeneric("defaultExclNormScores", function(obj) standardGeneric("defaultExclNormScores"))

#' @templateVar func export
#' @templateVar desc exports workflow data to a given format.
#' @template generics
#' @param type The export type.
setGeneric("export", function(obj, type, out, ...) standardGeneric("export"))

#' @templateVar func featureTable
#' @templateVar desc returns feature information.
#' @template generics
setGeneric("featureTable", function(obj, ...) standardGeneric("featureTable"))

#' @templateVar func filter
#' @templateVar desc provides various functionality to do post-filtering of data.
#' @template generics
setGeneric("filter", function(obj, ...) standardGeneric("filter"))

#' @templateVar func hasMobilities
#' @templateVar desc returns \code{TRUE} if the object was directly created from IMS data.
#' @template generics
setGeneric("fromIMS", function(obj) standardGeneric("fromIMS"))

#' @templateVar func getBPCs
#' @templateVar desc gets base peak chromatogram(s).
#' @template generics
setGeneric("getBPCs", function(obj, ...) standardGeneric("getBPCs"))

#' @templateVar func getFeatures
#' @templateVar desc returns the object's \code{\link{features}} object.
#' @template generics
setGeneric("getFeatures", function(obj) standardGeneric("getFeatures"))

#' @templateVar func getMCS
#' @templateVar desc Calculates the maximum common substructure.
#' @template generics
setGeneric("getMCS", function(obj, ...) standardGeneric("getMCS"))

#' @templateVar func getTICs
#' @templateVar desc gets total ion chromatogram(s).
#' @template generics
setGeneric("getTICs", function(obj, ...) standardGeneric("getTICs"))

#' @templateVar func groupNames
#' @templateVar desc returns a \code{character} vector with the names of the feature groups for which data is present in this object.
#' @template generics
setGeneric("groupNames", function(obj) standardGeneric("groupNames"))

#' @templateVar func hasMobilities
#' @templateVar desc returns \code{TRUE} if the object has ion mobility values
#' @template generics
setGeneric("hasMobilities", function(obj) standardGeneric("hasMobilities"))

#' @templateVar func plotBPCs
#' @templateVar desc plots base peak chromatogram(s).
#' @template generics
setGeneric("plotBPCs", function(obj, ...) standardGeneric("plotBPCs"))

#' @templateVar func plotChord
#' @templateVar desc plots a Chord diagram to assess overlapping data.
#' @template generics
#'
#' @template plotChord-args
setGeneric("plotChord", function(obj, addSelfLinks = FALSE, addRetMzPlots = TRUE, ...) standardGeneric("plotChord"))

#' @templateVar func plotChroms
#' @templateVar desc plots extracted ion chromatogram(s).
#' @template generics
setGeneric("plotChroms", function(obj, ...) standardGeneric("plotChroms"))

#' @templateVar func plotChroms3D
#' @templateVar desc plots a three dimensional chromatogram.
#' @template generics
setGeneric("plotChroms3D", function(obj, ...) standardGeneric("plotChroms3D"))

#' @templateVar func plotGraph
#' @templateVar desc Plots an interactive network graph.
#' @template generics
setGeneric("plotGraph", function(obj, ...) standardGeneric("plotGraph"))

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

#' @templateVar func plotTICs
#' @templateVar desc plots total ion chromatogram(s).
#' @template generics
setGeneric("plotTICs", function(obj, ...) standardGeneric("plotTICs"))

#' @templateVar func plotVenn
#' @templateVar desc plots a Venn diagram to assess unique and overlapping data.
#' @template generics
#'
setGeneric("plotVenn", function(obj, ...) standardGeneric("plotVenn"))

#' @templateVar func plotUpSet
#' @templateVar desc plots an UpSet diagram to assess unique and overlapping data.
#' @template generics
setGeneric("plotUpSet", function(obj, ...) standardGeneric("plotUpSet"))

#' @templateVar func predictRespFactors
#' @templateVar desc Prediction of response factors.
#' @template generics
setGeneric("predictRespFactors", function(obj, ...) standardGeneric("predictRespFactors"))

#' @templateVar func predictTox
#' @templateVar desc Prediction of toxicity values.
#' @template generics
setGeneric("predictTox", function(obj, ...) standardGeneric("predictTox"))

#' @templateVar func delete
#' @templateVar desc Deletes results.
#' @template generics
setGeneric("delete", function(obj, ...) standardGeneric("delete"))

#' @templateVar func plotVolcano
#' @templateVar desc plots a volcano plot.
#' @template generics
setGeneric("plotVolcano", function(obj, ...) standardGeneric("plotVolcano"))

#' @templateVar func replicates
#' @templateVar desc returns a \code{character} vector with the replicates for which data is present in this object.
#' @template generics
setGeneric("replicates", function(obj) standardGeneric("replicates"))

#' @templateVar func setObjects
#' @templateVar desc returns the \emph{set objects} of this object. See the documentation of \code{\link{workflowStepSet}}.
#' @template generics
setGeneric("setObjects", function(obj) standardGeneric("setObjects"))

#' @templateVar func sets
#' @templateVar desc returns the names of the sets inside this object. See the documentation for \link[=sets-workflow]{sets workflows}.
#' @template generics
setGeneric("sets", function(obj) standardGeneric("sets"))

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

#' @templateVar func unset
#' @templateVar desc Converts this object to a regular non-set object. See the documentation for \link[=sets-workflow]{sets workflows}.
#' @template generics
#' @param set The name of the set.
setGeneric("unset", function(obj, set) standardGeneric("unset"))


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

setGeneric("launchEICGUI", function(obj, ...) standardGeneric("launchEICGUI"))
