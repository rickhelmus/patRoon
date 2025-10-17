# SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#### Dependencies

#' @importFrom utils head tail modifyList setTxtProgressBar txtProgressBar write.csv write.table read.csv data getFromNamespace
#' @importFrom graphics axis barplot close.screen grconvertX grconvertY grid layout legend lines par plot.new points polygon rasterImage rect screen segments split.screen strwidth text title xinch yinch abline contour persp
#' @importFrom grDevices adjustcolor colorRampPalette dev.off png
#' @importFrom stats cutree dist hclust heatmap lm median rect.hclust sd setNames as.dendrogram order.dendrogram as.dist as.formula weighted.mean cor p.adjust t.test
#' @importFrom Rdpack reprompt
#' @importFrom magrittr %>%
#' @import methods
#' @import data.table
#' @import shiny
#' @import cluster
NULL # need this for doc generation

# Need to import these as functions generated with checkmate/withr don't always take namespace in to account
#' @importFrom checkmate makeAssertion vname
#' @importFrom withr defer
NULL

# UNDONE: rstudioapi optional?

# For Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib patRoon, .registration = TRUE
NULL


#' @include cache.R
NULL

# So it can be used as a S4 slot
setOldClass("hclust")
setOldClass(c("dissimilarity", "dist"))

# for method signatures
# doesn't work with devtools::load_all(): confuses xcmsRaw "[" operator (!?!?)
#setClassUnion("numChar", c("character", "numeric"))

setClassUnion("logChar", c("character", "logical"))

# this appears to fix the occasional error "cannot coerce type 'S4' to vector of type 'double'"
setGeneric("plot")

setGeneric("merge") # sometimes its not yet

#### Generics

#' @include generics.R
NULL

#### Document generation

#' Workflow solutions for mass-spectrometry based non-target analysis.
#'
#' \Sexpr[results=text,echo=FALSE]{packageDescription("patRoon", fields = "Description")}
#'
#' @section Package options:
#'
#'   The following package options (see \code{\link{options}}) can be set:
#'
#'   \itemize{
#'
#'   \item \code{patRoon.cache.mode}: A \code{character} setting the current caching mode: \code{"save"} and
#'   \code{"load"} will only save/load results to/from the cache, \code{"both"} (default) will do both and \code{"none"}
#'   to completely disable caching. This option can be changed anytime, which might be useful, for instance, to
#'   temporarily disable cached results before running a function.
#'
#'   \item \code{patRoon.cache.fileName}: a \code{character} specifying the name of the cache file (default is
#'   \file{cache.sqlite}).
#'
#'   \item \code{patRoon.MS.backends},\code{patRoon.MS.preferIMS},\code{patRoon.path.BrukerTIMS}:
#'   Options related to the \link[=msdata]{raw data interface}.
#'
#'   \item \code{patRoon.threads}: The number of threads to be used for parallelization. This is currently only used by
#'   the \link[=msdata]{raw data interface} and when the \code{piek} algorithm is used for \link[=getDefPeakParams]{peak
#'   detection}.
#'
#'   \item \code{patRoon.MP.maxProcs}: The maximum number of processes that should be initiated in parallel. A good
#'   starting point is the number of physical cores, which is the default as detected by
#'   \code{\link[parallel]{detectCores}}. This option is only used when \option{patRoon.MP.method="classic"}.
#'
#'   \item \code{patRoon.MP.method}: Either \code{"classic"} or \code{"future"}. The former is the default and uses
#'   \CRANpkg{processx} to execute multiple commands in parallel. When \code{"future"} the \code{\link{future.apply}}
#'   package is used for parallelization, which is especially useful for \emph{e.g.} cluster computing.
#'
#'   \item \code{patRoon.MP.futureSched}: Sets the \code{future.scheduling} function argument for
#'   \code{\link{future_lapply}}. Only used if \option{patRoon.MP.method="future"}.
#'
#'   \item \code{patRoon.MP.logPath}: The path used for logging of output from commands executed by multiprocess. Set to
#'   \code{FALSE} to disable logging.
#'
#'   \item \code{patRoon.path.pwiz}: The path in which the \command{ProteoWizard} binaries are installed. If unset an
#'   attempt is made to find this directory from the Windows registry and \option{PATH} environment variable.
#'
#'   \item \code{patRoon.path.GenForm}: The path to the \command{GenForm} executable. If not set (the default) the
#'   internal \code{GenForm} binary is used. Only set if you want to override the executable.
#'
#'   \item \code{patRoon.path.MetFragCL}: The complete file path to the MetFrag CL \file{jar} to be used by
#'   \code{\link{generateCompoundsMetFrag}}. Example: \code{"C:/MetFrag2.4.2-CL.jar"}.
#'
#'   \item \code{patRoon.path.MetFragCompTox}: The complete file path to the CompTox database \file{csv} file. See
#'   \code{\link{generateCompounds}} for more details.
#'
#'   \item \code{patRoon.path.MetFragPubChemLite}: The complete file path to the PubChemLite database \file{csv} file.
#'   See \code{\link{generateCompounds}} for more details.
#'
#'   \item \code{patRoon.path.SIRIUS}: The directory in which the SIRIUS binaries are installed. Used by all functions that interface with \command{SIRIUS}, such as \code{\link{generateFormulasSIRIUS}}
#'   and \code{\link{generateCompoundsSIRIUS}}. Example: \code{"C:/sirius-win64-3.5.1"}. Note that the location of the
#'   binaries differs for each operating system.
#'   are installed in different subdirectories for each location inside  this differs for each operating system
#'
#'   \item \code{patRoon.path.OpenMS}: The path in which the \command{OpenMS} binaries are installed.
#'
#'   \item \code{patRoon.path.obabel}: The path in which the \command{OpenBabel} binaries are installed.
#'
#'   \item \code{patRoon.path.BiotransFormer} The full file path to the \command{biotransformer} \file{.jar} command
#'   line utility. This needs to be set when \code{\link{generateTPsBioTransformer}} is used. For more details see
#'   \url{https://bitbucket.org/djoumbou/biotransformer/src/master}.
#'
#'   \item \code{patRoon.path.limits} A path to a customized \link[=limits]{limits YAML file}.
#'
#'   }
#'
#'   Most external dependencies are provided by \pkg{patRoonExt} or otherwise found in the system environment
#'   \option{PATH} variable. However, the \code{patRoon.path.*} options should be set if this fails or you want to
#'   override the location. The \code{\link{verifyDependencies}} function can be used to assess if dependencies are
#'   found.
#' 
"_PACKAGE"

#' Properties of sample analyses
#'
#' Properties for the sample analyses used in the workflow and utilities to automatically generate this information.
#'
#' In \pkg{patRoon} a \emph{sample analysis}, or simply \emph{analysis}, refers to a single MS analysis file (sometimes
#' also called \emph{sample} or \emph{file}). The \emph{analysis information} summarizes several properties for the
#' analyses, and is used in various steps throughout the workflow, such as \code{\link{findFeatures}}, averaging
#' intensities of feature groups and blank subtraction. The analysis information should be a \code{data.frame} or
#' \code{data.table} with a set of mandatory and optional columns (described below).
#'
#' @section Mandatory analysis information columns:
#'
#'   The following columns should be present in the analysis information:
#'
#' \itemize{
#'
#' \item \code{path_raw}, \code{path_centroid}, \code{path_profile}, \code{path_ims} Specifies the directory path for
#' the raw, centroided, profile and IMS data, respectively. See below for more details. At least one column should not
#' be empty for each row.
#'
#' \item \code{analysis} the file name \strong{without} extension and without directory path. Must be \strong{unique}
#' across all table rows.
#'
#' \item \code{replicate} name of the \emph{replicate}. Used to group analyses together that are
#' replicates of each other. Thus, the \code{replicate} column for all analyses considered to be belonging to the same
#' replicate should have an equal (but unique) value. Used for \emph{e.g.} averaging and
#' \code{\link[=filter,featureGroups-method]{filter}}.
#'
#' \item \code{blank} all analyses within this replicate are used by the \code{featureGroups} method of
#' \code{\link[=filter,featureGroups-method]{filter}} for blank subtraction. Multiple entries can be entered by
#' separation with a comma. May be empty (\code{""}) if no blank subtraction is desired.
#'
#' }
#'
#' @section Analysis paths, file types and file formats:
#'
#'   Depending on the workflow step, different \emph{file types} for the same analysis may be required.
#'
#' \itemize{
#'
#' \item \code{raw} Specifies the directory to raw HRMS files (\emph{e.g.} \file{.raw}, \file{.d}). This is used by
#' \emph{e.g.} \link[=MSConversion]{conversion of raw MS data} and the \link[=msdata]{OpenTIMS backend}.
#'
#' \item \code{centroid} Specifies the directory to centroided and exported HRMS files (\file{.mzML}, \file{.mzXML}).
#' These files are required by most feature finding algorithms.
#'
#' \item \code{profile} Specifies the directory to exported but not centroided (\emph{i.e.} profile) HRMS data files
#' (\file{.mzML}, \file{.mzXML}). This is currently only used by \code{\link{findFeaturesSAFD}}.
#'
#' \item \code{ims} Specifies the directory to exported IMS-HRMS data (\file{.mzML}). This is required in IMS workflows,
#' unless raw IMS-HRMS data is directly loaded with the \link[=msdata]{OpenTIMS backend}. See \emph{e.g.}
#' \code{\link[=assignMobilities_feat]{assignMobilities}} for more details.
#'
#' }
#'
#'   Some workflows may require multiple \emph{file formats} for a same \emph{file type}. In this case, the file formats
#'   should be stored within the same directory specified by the respective \code{path_*} column. For instance, if
#'   feature finding algorithms from \link[=findFeaturesOpenMS]{OpenMS} and \link[=findFeaturesEnviPick]{enviPick} are
#'   mixed then centroided \file{.mzML} and \file{.mzXML} files are needed, and files with both file formats must be
#'   stored in the directory specified by \code{path_centroid}.
#'
#'   If non-raw data files are not yet present and should be exported by \link[=MSConversion]{MS file conversion}, then
#'   \code{path_centroid}, \code{path_profile} and \code{path_ims} should specify the desired destination paths of the
#'   converted files.
#'
#' @section Optional columns and sample metadata:
#'
#'   The following columns may need to be present:
#'
#' \itemize{
#'
#' \item \code{conc} a numeric value specifying the 'concentration' for the analysis. This can be actually any kind of
#' numeric value such as exposure time, dilution factor or anything else which may be used to form a linear
#' relationship. This is used by the \code{\link[=as.data.table,featureGroups-method]{as.data.table}} method if
#' \code{regression=TRUE}. As of \pkg{patRoon} version 3.0, any other column than \code{"conc"} can be used by setting
#' its name with the \code{regression} argument.
#'
#' \item \code{norm_conc} a numeric value specifying the \emph{normalization concentration} for the analysis. See the
#' \verb{Feature intensity normalization} section in the \link[=featureGroups-class]{featureGroups documentation}) for
#' more details.
#'
#' }
#'
#'   Any other columns that are present will be added to the \code{features} and \code{featureGroups} objects as
#'   metadata. This metadata can be used \emph{e.g.} in various plotting and data subsetting functions.
#'
#' @name analysis-information
NULL

#' Clustering parameters
#'
#' Parameters for clustering data such as mass spectra and mobilograms.
#'
#' Different functionality within \pkg{patRoon} uses clustering to group similar data together, for instance, to average
#' mass spectra or mobilograms. A fast \code{C++} backend based on \CRANpkg{Rcpp} is used to perform the clustering.
#'
#' The clustering can be configured by the \emph{method} and \emph{window} parameter. The following clustering methods
#' are available: \itemize{
#'
#'   \item \code{"hclust"}: uses hierarchical clustering to find similar data points (using
#'   \href{https://github.com/cdalitz/hclust-cpp}{hclust-cpp}, which is based on the \CRANpkg{fastcluster} package).
#'
#'   \item \code{"distance"}: uses the distance between sorted data points to find close points.
#'
#'   \item \code{"bin"}: uses a simple binning approach to cluster similar mass peaks.
#'
#' }
#'
#' The \code{hclust} method may give more accurate results and was the default prior to \pkg{patRoon 3.0}, but is more
#' computationally demanding and generally unsuitable for IMS workflows due to excessive use of RAM. The \code{distance}
#' method is now default and suits most cases.
#'
#' The window parameter defines the clustering tolerance. For \code{method="hclust"} this corresponds to the cluster
#' height, for \code{method="distance"} this value is used to find nearby data points (+/- window) and for
#' \code{method="bin"} it corresponds to the bin width. Too small windows will prevent clustering close data points
#' (\emph{e.g.} resulting in split mass peaks in averaged spectra), whereas too big windows may cluster unrelated data
#' points together.
#'
#' @section Source: Averaging of mass spectra was originally based on algorithms from the
#'   \href{https://github.com/zeehio/msProcess}{msProcess} R package (now archived on CRAN).
#'
#' @references \addCitations{Rcpp} \cr\cr \addCitations{fastcluster} \cr\cr
#'
#' @name cluster-params
NULL

#' Optimization of feature finding and grouping parameters
#'
#' Automatic optimization of feature finding and grouping parameters through Design of Experiments (DoE).
#'
#' Many different parameters exist that may affect the output quality of feature finding and grouping. To avoid time
#' consuming manual experimentation, functionality is provided to largely automate the optimization process. The
#' methodology, which uses design of experiments (DoE), is based on the excellent
#' \href{https://github.com/rietho/IPO}{Isotopologue Parameter Optimization (IPO) R package}. The functionality of this
#' package is directly integrated in patRoon. Some functionality was added or changed, however, the principle algorithm
#' workings are nearly identical.
#'
#' Compared to IPO, the following functionality was added or changed: \itemize{
#'
#' \item The code was made more generic in order to include support for other feature finding/grouping algorithms
#' (\emph{e.g.} OpenMS, enviPick, XCMS3). \item The methodology of \command{FeatureFinderMetabo} (OpenMS) may be used to
#' find isotopes.
#'
#' \item The \code{maxModelDeviation} parameter was added to potentially avoid suboptimal results
#' (\href{https://github.com/rietho/IPO/issues/61}{issue discussed here}).
#'
#' \item The use of multiple 'parameter sets' (discussed below) which, for instance, allow optimizing qualitative
#' paremeters more easily (see \verb{examples}).
#'
#' \item More consistent optimization code for feature finding/grouping.
#'
#' \item More consistent output using S4 classes (\emph{i.e.} \code{\link{optimizationResult}} class).
#'
#' \item Parallelization is performed via the \CRANpkg{future} package instead of \pkg{BiocParallel}. If this is enabled
#' (\code{parallel=TRUE}) then any parallelization supported by the feature finding or grouping algorithm is disabled.
#'
#' }
#'
#'
#' @param anaInfo \link[=analysis-information]{Analysis info table} (passed to \code{\link{findFeatures}}).
#' @param algorithm The algorithm used for finding or grouping features (see \code{\link{findFeatures}} and
#'   \code{\link{groupFeatures}}).
#' @param \dots One or more lists with parameter sets (see below) (for \code{optimizeFeatureFinding} and
#'   \code{optimizeFeatureGrouping}). Alternatively, named arguments that set (and possibly override) the parameters
#'   that should be returned from \code{generateFeatureOptPSet} or \code{generateFGroupsOptPSet}.
#' @param templateParams Template parameter set (see below).
#' @param paramRanges A list with vectors containing absolute parameter ranges (minimum/maximum) that constrain numeric
#'   parameters choosen during experiments. See the \code{\link{getDefFeaturesOptParamRanges}} and
#'   \code{\link{getDefFGroupsOptParamRanges}} functions for defaults. Values should be \code{Inf} when no limit should
#'   be used.
#' @param maxIterations Maximum number of iterations that may be performed to find optimimum values. Used to restrict
#'   neededless long optimization procedures. In IPO this was fixed to \samp{50}.
#' @param maxModelDeviation See the \verb{Potential suboptimal results by optimization model} section below.
#'
#' @template parallel-arg
#'
#' @section Parameter sets: Which parameters should be optimized is determined by a \emph{parameter set}. A set is
#'   defined by a named \code{list} containing the minimum and maximum starting range for each parameter that should be
#'   tested. For instance, the set \code{list(chromFWHM = c(5, 10), mzPPM = c(5, 15))} specifies that the
#'   \code{chromFWHM} and \code{mzPPM} parameters (used by OpenMS feature finding) should be optimized within a range of
#'   \samp{5}-\samp{10} and \samp{5}-\samp{15}, respectively. Note that this range may be increased or decreased after a
#'   DoE iteration in order to find a better optimum. The absolute limits are controlled by the \code{paramRanges}
#'   function argument.
#'
#'   Multiple parameter sets may be specified (\emph{i.e.} through the \dots function argument). In this situation, the
#'   optimization algorithm is repeated for each set, and the final optimum is determined from the parameter set with
#'   the best response. The \code{templateParams} function argument may be useful in this case to define a template for
#'   each parameter set. Actual parameter sets are then constructed by joining each parameter set with the set specified
#'   for \code{templateParams}. When a parameter is defined in both a regular and template set, the parameter in the
#'   regular set takes precedence.
#'
#'   Parameters that should not be optimized but still need to be set for the feature finding/grouping functions should
#'   also be defined in a (template) parameter set. Which parameters should be optimized is determined whether its value
#'   is specified as a vector range or a single fixed value. For instance, when a set is defined as \code{list(chromFWHM
#'   = c(5, 10), mzPPM = 5)}, only the \code{chromFWHM} parameter is optimized, whereas \code{mzPPM} is kept constant at
#'   \samp{5}.
#'
#'   Using multiple parameter sets with differing fixed values allows optimization of qualitative values (see examples
#'   below).
#'
#'   The parameters specified in parameter sets are directly passed through the \code{\link{findFeatures}} or
#'   \code{\link{groupFeatures}} functions. Hence, grouping and retention time alignment parameters used by XCMS should
#'   (still) be set through the \code{groupArgs} and \code{retcorArgs} parameters.
#'
#'   \strong{NOTE:} For XCMS3, which normally uses parameter classes for settings its options, the parameters must be
#'   defined in a named list like any other algorithm. The set parameters are then used passed to the constructor of the
#'   right parameter class object (e.g. \code{\link{CentWaveParam}}, \code{\link{ObiwarpParam}}). For grouping/alignment
#'   sets, these parameters need to be specified in nested lists called \code{groupParams} and \code{retAlignParams},
#'   respectively (similar to \code{groupArgs}/\code{retcorArgs} for \code{algorithm="xcms"}). Finally, the underlying
#'   XCMS method to be used should be defined in the parameter set (\emph{i.e.} by setting the \code{method} field for
#'   feature parameter sets and the \code{groupMethod} and \code{retAlignMethod} for grouping/aligning parameter sets).
#'   See the examples below for more details.
#'
#'   \strong{NOTE:} Similar to IPO, the \code{peakwidth} and \code{prefilter} parameters for XCMS feature finding should
#'   be split in two different values: \itemize{
#'
#'   \item The minimum and maximum ranges for \code{peakwidth} are optimized by setting \code{min_peakwidth} and
#'   \code{max_peakwidth}, respectively.
#'
#'   \item The \code{k} and \code{I} parameters contained in \code{prefilter} are split in \code{prefilter} and
#'   \code{value_of_prefilter}, respectively.
#'
#'   }
#'
#'   \emph{Similary}, for KPIC2, the following parameters should be split: \itemize{
#'
#'   \item the \code{width} parameter (feature optimization) is optimized by specifying the \code{min_width} and
#'   \code{max_width} parameters.
#'
#'   \item the \code{tolerance} and \code{weight} parameters (feature grouping optimization) are optimized by setting
#'   \code{mz_tolerance}/\code{rt_tolerance} and \code{mz_weight}/\code{rt_weight} parameters, respectively.
#'
#'   }
#'
#'
#' @section Functions: The \code{optimizeFeatureFinding} and \code{optimizeFeatureGrouping} are the functions to be used
#'   to optimize parameters for feature finding and grouping, respectively. These functions are analogous to
#'   \code{\link[IPO]{optimizeXcmsSet}} and \code{\link[IPO]{optimizeRetGroup}} from \pkg{IPO}.
#'
#'   The \code{generateFeatureOptPSet} and \code{generateFGroupsOptPSet} functions may be used to generate a parameter
#'   set for feature finding and grouping, respectively. Some algorithm dependent default parameter optimization ranges
#'   will be returned. These functions are analogous to \code{\link[IPO]{getDefaultXcmsSetStartingParams}} and
#'   \code{\link[IPO]{getDefaultRetGroupStartingParams}} from \pkg{IPO}. However, unlike their IPO counterparts, these
#'   functions will not output default fixed values. The \code{generateFGroupsOptPSet} will only generate defaults for
#'   density grouping if \code{algorithm="xcms"}.
#'
#'   The \code{getDefFeaturesOptParamRanges} and \code{getDefFGroupsOptParamRanges} return the default absolute
#'   optimization parameter ranges for feature finding and grouping, respectively. These functions are useful if you
#'   want to set the \code{paramRanges} function argument.
#'
#' @section Potential suboptimal results by optimization model: After each experiment iteration an optimimum parameter
#'   set is found by generating a model containing the tested parameters and their responses. Sometimes the actual
#'   response from the parameters derived from the model is actually signficantly lower than expected. When the response
#'   is lower than the maximum reponse found during the experiment, the parameters belonging to this experimental
#'   maximum may be choosen instead. The \code{maxModelDeviation} argument sets the maximum deviation in response
#'   between the modelled and experimental maxima. The value is relative: \samp{0} means that experimental values will
#'   always be favored when leading to improved responses, whereas \code{1} will effectively disable this procedure (and
#'   return to 'regular' IPO behaviour).
#'
#' @return The \code{optimizeFeatureFinding} and \code{optimizeFeatureGrouping} return their results in a
#'   \code{\link{optimizationResult}} object.
#'
#' @section Source: The code and methodology is a direct adaptation from the \href{https://github.com/rietho/IPO}{IPO R
#'   package}.
#'
#' @references \insertRef{Libiseller2015}{patRoon}
#'
#' @example inst/examples/optimization.R
#'
#' @name feature-optimization
NULL

#' Identification confidence estimation
#'
#' Functions to estimate the identification confidence for suspects and annotation candidates.
#'
#' The \code{estimateIDConfidence} methods are used to estimate various properties to estimate the confidence of
#' identifications assigned to suspects and feature annotation candidates. These functions are typically executed after
#' running \code{\link{screenSuspects}}, \code{\link{generateFormulas}} and \code{\link{generateCompounds}}. Afterwards,
#' the following columns are added to the result tables (obtained with \emph{e.g.} \code{\link{screenInfo}},
#' \code{\link{annotations}} and \code{\link{as.data.table}}): \itemize{
#'
#'   \item \code{annSim} The \emph{annotation similarity}, defined as the similarity between the MS/MS peak list of a
#'   feature with (a) only the peaks that were annotated and (b) all the peaks. Thus, a value of one means that all
#'   MS/MS peaks were annotated. The similarity calculation is configured with the \code{specSimParams} argument to
#'   \code{estimateIDConfidence}.
#'
#'   \item \code{annSimForm} The annotation similarity specifically for formula annotations (equaling the \code{annSim}
#'   column from formula annotations). Only calculated for \link[=featureGroupsScreening]{suspects} and
#'   \code{\link{compounds}}.
#'
#'   \item \code{annSimBoth} The annotation similarity calculated with the combined set of annotated MS/MS peaks from
#'   formula and compound annotations. Only calculated for \link[=featureGroupsScreening]{suspects} and
#'   \code{\link{compounds}}.
#'
#'   \item \code{estIDLevel} Provides an \emph{estimation} of the identification level, roughly following that of
#'   \insertCite{Schymanski2014}{patRoon}. However, please note that this value is only an estimation, and manual
#'   interpretation is still necessary to assign final identification levels. The estimation is done through a set of
#'   rules, see the \verb{Identification level rules} section below.
#'
#' }
#'
#' In addition, the following columns are specifically added to suspect screening results: \itemize{
#'
#'   \item \code{annSimComp} The annotation similarity specifically for compound annotations (this equals the
#'   \code{annSim} column in compound annotations.
#'
#'   \item \code{formRank},\code{compRank} The rank of the suspect within the formula/compound annotation results.
#'
#'   \item \code{maxFrags} The maximum number of MS/MS fragments that can be matched for this suspect (based on the
#'   \code{fragments_*} columns from the suspect list).
#'
#'   \item \code{maxFragMatches},\code{maxFragMatchesRel} The absolute and relative amount of experimental MS/MS peaks
#'   that were matched from the fragments specified in the suspect list. The value for \code{maxFragMatchesRel} is
#'   relative to the value for \code{maxFrags}. The calculation of this column is influenced by the
#'   \code{checkFragments} argument to \code{estimateIDConfidence}.
#'
#' }
#'
#' The data for these columns is only calculated if \code{estimateIDConfidence} has the required data to do so. For
#' instance, \code{annSimForm} and \code{formRank} are only calculated if the \code{formulas} argument is set, and
#' levels for \code{estIDLevel} will be poor if no compound annotations are available.
#'
#' @templateVar normParam normalizeScores,compoundsNormalizeScores,formulasNormalizeScores
#' @templateVar noNone TRUE
#' @template norm-args
#'
#' @template specSimParams-arg
#' @template parallel-arg
#'
#' @param obj The object for which identification confidence should be estimated.
#' @param MSPeakLists,formulas,compounds Annotation data (\code{\link{MSPeakLists}}, \code{\link{formulas}} and
#'   \code{\link{compounds}}). All arguments can be \code{NULL}, but it is recommended to set them if possible to allow
#'   the most complete estimations.
#' @param absMzDev Maximum absolute \emph{m/z} deviation.
#' @param IDFile A file path to a YAML file with rules used for estimation of identification levels. See the
#'   \verb{Suspect annotation} section for more details. If not specified then a default rules file will be used.
#' @param logPath A directory path to store logging information. If \code{NULL} then logging is disabled. \strong{NOTE}:
#'   To avoid slowdowns by logging for potentially large number of candidates, logging is disabled for the
#'   \code{\link{formulas}} and \code{\link{compounds}} methods by default.
#'
#' @section Identification level rules: The estimation of identification levels is configured through a YAML file which
#'   specifies the rules for each level. The default file is shown below.
#'
#' @eval paste0("@@section Identification level rules: \\preformatted{", patRoon:::readAllFile(system.file("misc",
#'   "IDLevelRules.yml", package = "patRoon")), "}")
#'
#' @section Identification level rules: Most of the file should be self-explanatory. Some notes:
#'
#'   \itemize{
#'
#'   \item Each rule is either a field of \code{suspectFragments} (minimum number of MS/MS fragments matched from
#'   suspect list), \code{retention} (maximum retention deviation from suspect list), \code{rank} (the maximum
#'   annotation rank from formula or compound annotations), \code{all} (this level is always matched) or any of the
#'   scorings available from the formula or compound annotations.
#'
#'   \item In case any of the rules could be applied to either formula or compound annotations, the annotation type must
#'   be specified with the \code{type} field (\code{formula} or \code{compound}).
#'
#'   \item Identification levels should start with a number and may optionally be followed by a alphabetic character.
#'   The lowest levels are checked first.
#'
#'   \item If \code{relative=yes} then the relative scoring will be used for testing.
#'
#'   \item For \code{suspectFragments}: if the number of fragments from the suspect list (\code{maxFrags} column) is
#'   less then the minimum rule value, the minimum is adjusted to the number of available fragments.
#'
#'   \item The \code{or} and \code{and} keywords can be used to combine multiple conditions.
#'
#'   \item Any conditions that require suspect data (\emph{e.g.} \code{suspectFragments}) are only met with the suspects
#'   method for \code{estimateIDConfidence} method.
#'
#'   }
#'
#'   A template rules file can be generated with the \code{\link{genIDLevelRulesFile}} function, and this file can
#'   subsequently passed to \code{estimateIDConfidence}. The file format is highly flexible and (sub)levels can be added
#'   or removed if desired. Note that the default file is currently only suitable when annotation is performed with
#'   \command{GenForm} and \command{MetFrag}, for other algorithms it is crucial to modify the rules.
#'
#' @section Sets workflows: \code{estimateIDConfidence} performs its estimations per set. In addition, the
#'   following overall (not set specific) columns are calculated: \itemize{
#'
#'     \item \code{formRank} and \code{compRank} based on the ranking of the formula/compound in the set consensus data.
#'
#'     \item \code{estIDLevel}: based on the 'best' estimated identification level among the sets data (\emph{i.e.} the
#'     lowest). In case there is a tie between sub-levels (\emph{e.g.} \samp{3a} and \samp{3b}), then the sub-level is
#'     stripped (\emph{e.g.} \samp{3}).
#'
#'     \item Annotation similarities: taken as the maximum value from the data for each set.
#'
#'   }
#'
#' @return \code{estimateIDConfidence} amends the input object with aforementioned identification confidence properties.
#'
#' @author Rick Helmus <\email{r.helmus@@uva.nl}>, Emma Schymanski <\email{emma.schymanski@@uni.lu}> (contributions to
#'   identification level rules), Bas van de Velde (contributions to spectral similarity calculation).
#'
#' @references \insertAllCited{} \cr \cr \insertRef{Stein1994}{patRoon}
#'
#' @aliases estimateIDConfidence
#' @name id-conf
NULL

#' Sets workflows
#'
#' With sets workflows in \pkg{patRoon} a complete non-target (or suspect) screening workflow is performed with sample
#' analyses that were measured with different MS methods (typically positive and negative ionization).
#'
#' The analyses files that were measured with a different method are grouped in \emph{sets}. In the most typical case,
#' there is a \code{"positive"} and \code{"negative"} set, for the positively/negatively ionized data, respectively.
#' However, other distinctions than polarity are also possible (although currently the chromatographic method should be
#' the same between sets). A sets workflow is typically initiated with the \code{\link{makeSet}} method. The handbook
#' contains much more details about sets workflows.
#'
#' @seealso \code{\link{makeSet}} to initiate sets workflows, \code{\link{workflowStepSet}}, the \verb{Sets workflows}
#'   sections in other documentation pages and the \pkg{patRoon} handbook.
#' @name sets-workflow
NULL

#' Initiate sets workflows
#'
#' Initiate sets workflows from specified feature data.
#'
#' The \code{makeSet} method function is used to initiate a \link[=sets-workflow]{sets workflow}. The features from
#' input objects are combined and then \emph{neutralized} by replacing their \emph{m/z} values by neutral monoisotopic
#' masses. After neutralization features measured with \emph{e.g.} different ionization polarities can be grouped since
#' their neutral mass will be the same.
#'
#' The \link[=analysis-information]{analysis information} for this object is updated with all analyses, and a \code{set}
#' column is added to designate the set of each analysis. Note that currently, all analyses names \strong{must} be
#' unique across different sets.
#'
#' \code{makeSet} supports two types of input: \enumerate{
#'
#' \item \code{\link{features}} objects: \code{makeSet} combines the input objects into a \code{featuresSet} object,
#' which is then grouped in the 'usual way' with \code{\link{groupFeatures}}.
#'
#' \item \code{\link{featureGroups}} objects: In this case the features from the input objects are first neutralized and
#' feature groups between sets are then combined with \code{groupFeatures}.
#'
#' }
#'
#' The advantage of the \code{featureGroups} method is that it preserves any adduct annotations already present
#' (\emph{e.g.} as set by \code{selectIons} or \code{adducts<-}). Furthermore, this approach allows more advanced
#' workflows where the input \code{featureGroups} are first pre-treated with \emph{e.g.} filter before the sets object
#' is made. On the other hand, the \code{features} method is easier, as it doesn't require intermediate feature grouping
#' steps and is often sufficient since adduct annotations can be made afterwards with \code{selectIons}/\code{adducts<-}
#' and most \code{filter} operations do not need to be done per individual set.
#'
#' The adduct information used for feature neutralization is specified through the \code{adducts} argument.
#' Alternatively, when the \code{featureGroups} method of \code{makeSet} is used, then the adduct annotations already
#' present in the input objects can also by used by setting \code{adducts=NULL}. The adduct information is also used to
#' add adduct annotations to the output of \code{makeSet}.
#'
#' @param obj,\dots \code{\link{features}} or \code{\link{featureGroups}} objects that should be used for the
#'   \link[=sets-workflow]{sets workflow}.
#' @param adducts The adduct assignments to each set. Should either be a \code{list} with \code{\link{adduct}} objects
#'   or a \code{character} vector (\emph{e.g.} \code{"[M+H]+"}). The order should follow that of the objects given to
#'   the \code{obj} and \code{\dots} arguments.
#'
#'   For the \code{\link{featureGroups}} method: if \code{NULL} then adduct annotations are used.
#' @param labels The labels, or \emph{set names}, for each set to be created. The order should follow that of the
#'   objects given to the \code{obj} and \code{\dots} arguments. If \code{NULL}, then labels are automatically generated
#'   from the polarity of the specified \code{adducts} argument (\emph{e.g.} \code{"positive"}, \code{"negative"}).
#'
#' @note Initiating a sets workflow recursively, \emph{i.e.} with \code{featuresSet} or \code{featureGroupsSet} objects
#'   as input, is currently not supported.
#'
#' @return Either a \code{\link{featuresSet}} object (\code{features} method) or \code{\link{featureGroupsSet}} object
#'   (\code{featureGroups} method).
#'
#' @name makeSet
NULL

#' Bruker DataAnalysis utilities
#'
#' Miscellaneous utility functions which interface with Bruker DataAnalysis
#'
#' These functions communicate directly with Bruker DataAnalysis to provide
#' various functionality, such as calibrating and exporting data and adding
#' chromatographic traces. For this the \pkg{RDCOMClient} package is
#' required to be installed.
#'
#' @param anaInfo \link[=analysis-information]{Analysis info table}
#' @param mzWindow \emph{m/z} window (in Da) used for the chromatographic trace
#'   (if applicable).
#' @param ctype Type of the chromatographic trace. Valid options are:
#'   \code{"EIC"} (extracted ion chromatogram), \code{"TIC"} (total ion
#'   chromatogram, only for \code{addDAEIC}) and \code{"BPC"} (Base Peak
#'   Chromatogram).
#' @param bgsubtr If \code{TRUE} then background subtraction ('Spectral'
#'   algorithm) will be performed.
#' @param name For \code{addDAEIC}: the name for the chromatographic trace. For
#'   \code{addAllEICs}: \code{TRUE} to automatically set EIC names. Set to
#'   \code{NULL} for none.
#' @param hideDA Hides DataAnalysis while adding the chromatographic trace
#'   (faster).
#'
#' @template dasaveclose-args
#'
#' @name bruker-utils
#'
#' @seealso \code{\link{analysis-information}}
NULL

#' Interactive GUI utilities to check workflow data
#'
#' These functions provide interactive utilities to explore and review workflow data using a \pkg{\link{shiny}}
#' graphical user interface (GUI). In addition, unsatisfactory data (\emph{e.g.} noise identified as a feature and
#' unrelated feature groups in a component) can easily be selected for removal.
#'
#' The data selected for removal is stored in \emph{sessions}. These are \file{YAML} files to allow easy external
#' manipulation. The sessions can be used to restore the selections that were made for data removal when the GUI tool is
#' executed again. Furthermore, functionality is provided to import and export sessions. To actually remove the data the
#' \code{\link{filter}} method should be used with the session file as input.
#'
#' @param fGroups A \code{\link{featureGroups}} object.
#'
#'   This should be the 'new' object for \code{importCheckFeaturesSession} for which the session needs to be imported.
#' @param session The session file name.
#' @param clearSession If \code{TRUE} the session will be completely cleared before starting the GUI. This effectively
#'   removes all selections for data removal.
#'
#' @template EICParams-arg
#' 
#' @templateVar what these functions
#' @template uses-msdata
#' 
#' @note The \code{topMost} and \code{topMostByReplicate} \link[=EIXParams]{EIC parameters} are ignored.
#'
#' @references \insertRef{Chetnik2020}{patRoon}
#'
#' @name check-GUI
NULL

#' Retention order direction
#'
#' Calculation of the relative retention order between a parent and its transformation product (TP).
#'
#' The relative retention order between a parent and its TP (\code{retDir}) is used throughout TP screening workflows
#' for characterization and prioritization purposes. These are \code{numeric} values that hint what the the
#' chromatographic retention order of a TP might be compared to its parent: a value of \samp{-1} means it will elute
#' earlier, \samp{1} it will elute later and \samp{0} that there is no significant difference or the direction is
#' unknown.
#'
#' For TP data obtained with \code{\link{generateTPs}}, the missing \code{retDir} values are automatically calculated
#' based on the \verb{log P} difference between the parent and TP. Here, a typical reversed phase separation is assumed,
#' \emph{i.e.} compounds with (significantly) lower log P values likely elute earlier. The \code{minLogPDiff} parameter
#' of the \link[=getDefTPStructParams]{TPStructParams} argument sets the minimum log P difference to be considered
#' significant.
#'
#' For TP feature candidates that were linked by \code{\link{generateComponentsTPs}}, the \code{retDir} values are
#' calculated based on the retention time difference between the parent and TP feature groups. The \code{minRTDiff}
#' argument sets the minimum difference to be considered significant.
#'
#' @references \insertRef{Helmus2025}{patRoon}
#' 
#' @name retDir
NULL

#' Functionality to predict quantitative data
#'
#' Functions to predict response factors and feature concentrations from \acronym{SMILES} and/or
#' \command{SIRIUS+CSI:FingerID} fingerprints using the \pkg{MS2Quant} package.
#'
#' The \href{https://github.com/kruvelab/MS2Quant}{MS2Quant} \R package predicts concentrations from \acronym{SMILES}
#' and/or MS/MS fingerprints obtained with \command{SIRIUS+CSI:FingerID}. The \code{predictRespFactors} method functions
#' interface with this package to calculate response factors, which can then be used to calculate feature concentrations
#' with the \code{calculateConcs} method function.
#'
#' @param fGroups For \code{predictRespFactors} methods for feature annotations: The \code{\link{featureGroups}} object
#'   for which the annotations were performed.
#'
#'   For \code{calculateConcs}: The \code{\link{featureGroups}} object for which concentrations should be calculated.
#'
#'   For \code{getQuantCalibFromScreening}: A feature groups object screened for the calibrants with
#'   \code{\link{screenSuspects}}.
#' @param areas Set to \code{TRUE} to use peak areas instead of peak heights. Note: for \code{calculateConcs} this
#'   should follow what is in the \code{calibrants} table.
#' @param calibrants A \code{data.frame} with calibrants, see the \verb{Calibration} section below.
#'
#'   \setsWF Should be a \code{list} with the calibrants for each set.
#' @param eluent A \code{data.frame} that describes the LC gradient program. Should have a column \code{time} with the
#'   retention time in seconds and a column \code{B} with the corresponding percentage of the organic modifier
#'   (\samp{0-100}).
#' @param organicModifier The organic modifier of the mobile phase: either \code{"MeOH"} (methanol) or \code{"MeCN"}
#'   (acetonitrile).
#' @param pHAq The \acronym{pH} of the aqueous part of the mobile phase.
#' @param calibConcUnit The concentration unit used in the calibrants table. For possible values see the \code{concUnit}
#'   argument.
#' @param concs A \code{data.frame} with concentration data. See the \verb{Calibration} section below.
#' @param average Set to \code{TRUE} to average intensity values within replicates.
#'
#' @templateVar scoreName response factor
#' @templateVar scoreWeightName scoreWeight
#' @template update_comp_score-args
#'
#' @template parallel-arg
#'
#' @section Calibration: The \pkg{MS2Quant} package requires calibration to convert predicted ionization efficiencies to
#'   instrument/method specific response factors. The calibration data should be specified with the \code{calibrants}
#'   argument to \code{predictRespFactors}. This should be a \code{data.frame} with intensity observations at different
#'   concentrations for a set of calibrants. Each row specifies one intensity observation at one concentration. The
#'   table should have the following columns:
#'
#'   \itemize{
#'
#'   \item \code{name} The name of the calibrant. Can be freely chosen.
#'
#'   \item \code{SMILES} The \acronym{SMILES} of the calibrant.
#'
#'   \item \code{rt} The retention time of the calibrant (in seconds).
#'
#'   \item \code{intensity} The peak intensity (or area, see the \code{areas} argument) of the calibrant.
#'
#'   \item \code{conc} The concentration of the calibrant (see the \code{calibConcUnit} argument for specifying the unit).
#'
#'   }
#'
#'   It is recommended to include multiple calibrants (\emph{e.g.} \samp{>=10}) at multiple concentrations (\emph{e.g.}
#'   \samp{>=5}). The latter is achieved by adding multiple rows for the same calibrant (keeping the
#'   \code{name}/\code{SMILES}/\code{rt} columns constant). It is also possible to follow the column naming used by
#'   \pkg{MS2Quant} (however retention times should still be in seconds!). For more details and tips see
#'   \url{https://github.com/kruvelab/MS2Quant}.
#'
#'   The \code{getQuantCalibFromScreening} function can be used to automatically generate a calibrants table from a
#'   feature groups object with suspect screening results. Here, the idea is to perform a screening with
#'   \code{\link{screenSuspects}} with a suspect list that contain the calibrants, which is then used to construct the
#'   calibrant table. It is highly recommended to add retention times for the calibrants in the suspect list to ensure
#'   the calibrant is assigned to the correct feature. Furthermore, it is possible to simply add the calibrants to the
#'   'regular' suspect list in case a suspect screening was already part of the workflow. The
#'   \code{getQuantCalibFromScreening} function still requires you to specify concentration data, which is achieved via
#'   the \code{concs} argument. This should be a \code{data.frame} with a column \code{name} corresponding to the
#'   calibrant name (\emph{i.e.} same as used by \code{screenSuspects} above) and columns with concentration data. The
#'   latter columns specify the concentrations of a calibrant in different replicates (as defined in the
#'   \link[=analysis-information]{analysis information}). The concentration columns should be named after the
#'   corresponding replicate. Only those replicates that should be used for calibration need to be included.
#'   Furthermore, \code{NA} values can be used if a replicate should be ignored for a specific calibrant.
#'
#' @templateVar whatPred response factors
#' @templateVar predFunc predictRespFactors
#' @templateVar whatCalc concentrations
#' @templateVar calcFunc calculateConcs
#' @template pred-desc
#'
#' @section Assigning concentrations: In IMS workflows with post mobility assignments (see
#'   \code{\link[=assignMobilities_feat]{assignMobilities}}), the intensities of \emph{IMS parents} are used for the
#'   calculation of concentrations for mobility features (as calibration typically does not consider differences due to
#'   mobility filtering). However, the response factor assigned to mobility features are still used. This may yield
#'   small differences compared to workflows where concentrations are assigned prior to mobility assignments, as
#'   \code{assignMobilities} simply copies concentrations from IMS parents to mobility features.
#'
#' @return \code{predictRespFactors} returns an object amended with response factors (\code{RF_SMILES}/\code{LRF_SIRFP}
#'   columns).
#'
#'   \code{calculateConcs} returns a \code{\link{featureGroups}} based object amended with concentrations for each
#'   feature group (accessed with the \code{\link{concentrations}} method).
#'
#' @note \pkg{MS2Quant} currently \emph{only} supports \samp{M+H} and \samp{M+} adducts when performing predictions with
#'   \command{SIRIUS:FingerID} fingerprints. Predictions for candidates with other adducts, including \samp{M-H}, are
#'   skipped with a warning.
#'
#' @references \insertRef{Sepman2023}{patRoon}
#'
#' @seealso \link[=pred-tox]{Toxicity prediction}
#'
#' @name pred-quant
NULL

#' Functionality to predict toxicities
#'
#' Functions to predict toxicities from \acronym{SMILES} and/or \command{SIRIUS+CSI:FingerID} fingerprints using the
#' \pkg{MS2Tox} package.
#'
#' The \href{https://github.com/kruvelab/MS2Tox}{MS2Tox} \R package predicts toxicities from \acronym{SMILES} and/or
#' MS/MS fingerprints obtained with \command{SIRIUS+CSI:FingerID}. The \code{predictTox} method functions interface with
#' this package to predict toxicities, which can then be assigned to feature groups with the \code{calculateTox} method
#' function.
#'
#' @param fGroups For \code{predictTox} methods for feature annotations: The \code{\link{featureGroups}} object for
#'   which the annotations were performed.
#'
#'   For \code{calculateTox}: The \code{\link{featureGroups}} object for which toxicities should be assigned.
#' @param LC50Mode The mode used for predictions: should be \code{"static"} or \code{"flow"}.
#'
#' @templateVar scoreName response factor
#' @templateVar scoreWeightName scoreWeight
#' @template update_comp_score-args
#' 
#' @template parallel-arg
#'
#' @templateVar whatPred toxicities
#' @templateVar predFunc predictTox
#' @templateVar whatCalc toxicities
#' @templateVar calcFunc calculateTox
#' @template pred-desc
#'
#' @return \code{predictTox} returns an object amended with LC 50 values (\code{LC50_SMILES}/\code{LC50_SIRFP} columns).
#'
#'   \code{calculateTox} returns a \code{\link{featureGroups}} based object amended with toxicity values for each
#'   feature group (accessed with the \code{\link{toxicities}} method).
#'
#' @references \insertRef{Peets2022}{patRoon}
#'
#' @seealso \link[=pred-quant]{Concentration prediction}
#' 
#' @name pred-tox
NULL
