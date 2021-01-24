#### Dependencies

#' @importFrom utils head tail modifyList setTxtProgressBar txtProgressBar write.csv write.table read.csv
#' @importFrom graphics plot axis barplot close.screen grconvertX grconvertY grid layout legend lines par plot.new points polygon rasterImage rect screen segments split.screen strwidth text title xinch yinch abline contour persp
#' @importFrom grDevices adjustcolor colorRampPalette dev.off pdf png
#' @importFrom stats cutree dist hclust heatmap lm median rect.hclust sd setNames as.dendrogram order.dendrogram as.dist as.formula
#' @importFrom Rdpack reprompt
#' @importFrom magrittr %>%
#' @import data.table
#' @import shiny
#' @import cluster
#' @import ggplot2
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
setOldClass("dissimilarity")

# for method signatures
# doesn't work with devtools::load_all(): confuses xcmsRaw "[" operator (!?!?)
#setClassUnion("numChar", c("character", "numeric"))

#### Generics

#' @include generics.R
NULL

#### Document generation

#' Workflow solutions for mass-spectrometry based non-target analysis.
#'
#' \Sexpr[results=text,echo=FALSE]{packageDescription("patRoon", fields =
#' "Description")}
#'
#' @section Package options:
#'
#'   The following package options (see \code{\link{options}}) can be set:
#'
#'   \itemize{
#'
#'   \item \code{patRoon.cache.mode}: A \code{character} setting the current
#'   caching mode: \code{"save"} and \code{"load"} will only save/load results
#'   to/from the cache, \code{"both"} (default) will do both and \code{"none"}
#'   to completely disable caching. This option can be changed anytime, which
#'   might be useful, for instance, to temporarily disable cached results before
#'   running a function.
#'
#'   \item \code{patRoon.cache.fileName}: a \code{character} specifying the name
#'   of the cache file (default is \file{cache.sqlite}).
#'
#'   \item \code{patRoon.MP.maxProcs}: The maximum number of processes that
#'   should be initiated in parallel. A good starting point is the number of
#'   physical cores, which is the default as detected by
#'   \code{\link[parallel]{detectCores}}. This option is only used when
#'   \option{patRoon.MP.method="classic"}.
#'   
#'   \item \code{patRoon.MP.method}: Either \code{"classic"} or \code{"future"}.
#'   The former is the default and uses \code{\link{processx}} to execute
#'   multiple commands in parallel. When \code{"future"} the
#'   \code{\link{future.apply}} package is used for parallelization, which is
#'   especially useful for \emph{e.g.} cluster computing.
#'   
#'   \item \code{patRoon.MP.futureSched}: Sets the \code{future.scheduling} function
#'   argument for \code{\link{future_lapply}}. Only used if
#'   \option{patRoon.MP.method="future"}.
#'   
#'   \item \code{patRoon.MP.logPath}: The path used for logging of output from
#'   commands executed by multiprocess. Set to \code{FALSE} to disable logging.
#'   
#'   \item \code{patRoon.path.pwiz}: The path in which the
#'   \command{ProteoWizard} binaries are installed. If unset an attempt is made
#'   to find this directory from the Windows registry and \option{PATH}
#'   environment variable.
#'
#'   \item \code{patRoon.path.GenForm}: The path to the \command{GenForm}
#'   executable. If not set (the default) the internal \code{GenForm} binary is
#'   used. Only set if you want to override the executable.
#'
#'   \item \code{patRoon.path.MetFragCL}: The complete file path to the MetFrag
#'   CL \file{jar} file that \emph{must} be set when using
#'   \code{\link{generateCompoundsMetFrag}}. Example:
#'   \code{"C:/MetFrag2.4.2-CL.jar"}.
#'
#'   \item \code{patRoon.path.MetFragCompTox}: The complete file path to the
#'   CompTox database \file{csv} file. See \code{\link{generateCompounds}} for
#'   more details.
#'   
#'   \item \code{patRoon.path.MetFragPubChemLite}: The complete file path to the
#'   PubChem database \file{csv} file. See \code{\link{generateCompounds}} for
#'   more details.
#'
#'   \item \code{patRoon.path.SIRIUS}: The directory in which SIRIUS is
#'   installed. Unless the binaries can be located via the \option{PATH}
#'   environment variable, this \emph{must} be set when using
#'   \code{\link{generateFormulasSIRIUS}} or
#'   \code{\link{generateCompoundsSIRIUS}}. Example:
#'   \code{"C:/sirius-win64-3.5.1"}.
#'
#'   \item \code{patRoon.path.OpenMS}: The path in which the \command{OpenMS}
#'   binaries are installed. Usually the location is added to the \option{PATH}
#'   environment variable when OpenMS is installed, in which case this option
#'   can be left empty.
#'
#'   \item \code{patRoon.path.pngquant}: The path of the \command{pngquant}
#'   binary that is used when optimizing \file{.png} plots generated by
#'   \code{\link{reportHTML}} (with \code{optimizePng} set to \code{TRUE}). If the
#'   binary can be located through the \option{PATH} environment variable this
#'   option can remain empty. Note that some of the functionality of
#'   \code{reportHTML} only locates the binary through the \option{PATH}
#'   environment variable, hence, it is recommended to set up \option{PATH}
#'   instead.
#'
#'   \item \code{patRoon.path.obabel}: The path in which the \command{OpenBabel}
#'   binaries are installed. Usually the location is added to the \option{PATH}
#'   environment variable when OpenBabel is installed, in which case this option
#'   can be left empty.
#'   
#'   }
#'
"_PACKAGE"

#' Analysis information
#'
#' Required information for analyses that should be processed and utilities to
#' automatically generate this information.
#'
#' @details Several properties need to be known about analyses that should be
#'   processed during various workflow steps such as
#'   \link[=feature-finding]{finding features}, averaging intensities of feature
#'   groups and blank subtraction. This information should be made available
#'   with an 'analysis info' object, which is a \code{data.frame} containing the
#'   following columns:
#'
#' @template analysisInfo-list
#'
#' @details Most functionality requires the data files to be in either the
#'   \file{.mzXML} or \file{.mzML} format. Functionality that utilizes Bruker
#'   DataAnalysis (\emph{e.g.} \code{\link{findFeaturesBruker}}) may only work
#'   with data files in the proprietary \file{.d} format. Therefore, when tools
#'   with varying requirements are mixed (a common scenario), the data files
#'   also need to be present in the required formats. To deal with this
#'   situation the data files with varying formats should all be placed in the
#'   same path that was used to specify the location of the analysis.
#'
#'   Whether analysis data files are present in multiple formats or not, each
#'   analysis should only be entered \emph{once}. The \code{analysis} column is
#'   used as a basename to automatically find back the data file with the
#'   required format, hence, analysis names should be specified as the file name
#'   without its file extension.
#'
#'   The \code{group} column is \emph{mandatory} and needs to be filled in for
#'   each analysis. The \code{blank} column should also be present, however, these
#'   may contain empty character strings (\code{""}) for analyses where no blank
#'   subtraction should occur. The \code{conc} column is only required when
#'   obtaining regression information is desired with the
#'   \code{\link[=as.data.table,featureGroups-method]{as.data.table}} method.
#'
#' @name analysis-information
NULL


#' Finding features
#'
#' Functions and classes for collection of features.
#'
#' Several functions exist to collect features (\emph{i.e.} retention and MS
#' information that represent potential compounds) from a set of analyses. All
#' 'feature finders' return an object derived from the \code{\link{features}}
#' base class. The next step in a general workflow is to group and align these
#' features across analyses by \link[=feature-grouping]{feature groupers}. Note
#' that some feature finders have a plethora of options which sometimes may have
#' a large effect on the quality of results. Fine-tuning parameters is therefore
#' important, and the optimum is largely dependent upon applied analysis
#' methodology and instrumentation.
#'
#' @param analysisInfo \link[=analysis-information]{Analysis info table}.
#' @param ... further parameters passed to \code{\link[xcms]{xcmsSet}}
#'   (\code{findFeaturesXCMS}), \code{\link[enviPick]{enviPickwrap}}
#'   (\code{featurefinderEnviPick}) or to selected feature finding or importing
#'   algorithms (\code{findFeatures} and \code{importFeatures}).
#' @param verbose If set to \code{FALSE} then no text output is shown.
#'
#' @templateVar what \code{generateFeaturesOpenMS}
#' @template uses-multiProc
#'
#' @note The file format of analyses for \code{findFeaturesXCMS} and
#'   \code{findFeaturesXCMS3} must be \code{mzML} or \code{mzXML}.
#'
#' @name feature-finding
#' @return An object of a class which is derived from \code{\link{features}}.
#' @seealso \code{\link{features-class}} and \code{\link{analysis-information}}
NULL

#' Grouping of features
#'
#' Functions and classes for grouping of features across analyses.
#'
#' After \link[=feature-finding]{features have been found} the logical next step
#' is to align and group them across analyses. This process is necessary to
#' allow comparison of features between multiple analyses, which otherwise would
#' be difficult due to small deviations in retention and mass data. Thus,
#' algorithms of 'feature groupers' are used to collect features with similar
#' retention and mass data. In addition, advanced retention time alignment
#' algorithms exist to enhance grouping of features even with relative large
#' retention time deviations (\emph{e.g.} possibly observed from analyses
#' collected over a long period). Like \link[=feature-finding]{finding of
#' features}, various algorithms are supported which may have many parameters
#' that can be fine-tuned. This fine-tuning is likely to be necessary, since
#' optimal settings often depend on applied methodology and instrumentation.
#'
#' @param feat The \code{\link{features}} to be grouped.
#'   \code{importFeatureGroupsBrukerPA} and \code{importFeatureGroupsEnviMass}
#'   only support features generated by \code{\link{findFeaturesBruker}} and
#'   \code{\link{importFeaturesEnviMass}}, respectively.
#' @param analysisInfo A \code{data.frame} with \link[=analysis-information]{analysis info}.
#' @param rtalign Enable retention time alignment.
#' @param \dots Any parameters to be passed to the selected grouping/importing
#'   algorithm.
#' @param verbose if \code{FALSE} then no text output will be shown.
#' @param type Which file type should be imported or exported: \code{"brukerpa"}
#'   (Bruker ProfileAnalysis), \code{"brukertasq"} (Bruker TASQ),
#'   \code{envimass} (\pkg{enviMass}, only import) or \code{"mzmine"} (MZMine,
#'   only export).

#' @name feature-grouping
#' @return An object of a class which is derived from
#'   \code{\link{featureGroups}}.
#' @seealso \code{\link{featureGroups-class}}
NULL

#' Optimization of feature finding and grouping parameters
#'
#' Automatic optimization of feature finding and grouping parameters through
#' Design of Experiments (DoE).
#'
#' Many different parameters exist that may affect the output quality of feature
#' finding and grouping. To avoid time consuming manual experimentation,
#' functionality is provided to largely automate the optimization process. The
#' methodology, which uses design of experiments (DoE), is based on the
#' excellent \href{https://github.com/rietho/IPO}{Isotopologue Parameter
#' Optimization (IPO) R package}. The functionality of this package is directly
#' integrated in patRoon. Some functionality was added or changed, however, the
#' principle algorithm workings are nearly identical.
#'
#' Compared to IPO, the following functionality was added or changed:
#' \itemize{
#'   \item The code was made more generic in order to include support for other
#' feature finding/grouping algorithms (\emph{e.g.} OpenMS, enviPick, XCMS3).
#'   \item The methodology of \command{FeatureFinderMetabo} (OpenMS) may be used
#'   to find isotopes.
#'   \item The
#' \code{maxModelDeviation} parameter was added to potentially avoid
#' suboptimal results (\href{https://github.com/rietho/IPO/issues/61}{issue
#' discussed here}).
#'   \item The use of multiple 'parameter sets' (discussed
#' below) which, for instance, allow optimizing qualitative paremeters more
#' easily (see \verb{examples}).
#'   \item More consistent optimization code for feature
#' finding/grouping.
#'   \item More consistent output using S4 classes (\emph{i.e.}
#' \code{\link{optimizationResult}} class).
#'   \item Experiments are not (yet) executed in parallel (although feature
#'   finding or grouping may be if the algorithm supports it).
#' }
#'
#'
#' @param anaInfo \link[=analysis-information]{Analysis info table} (passed
#'   to \code{\link{findFeatures}}).
#' @param algorithm The algorithm used for finding or grouping features (see
#'   \code{\link{findFeatures}} and \code{\link{groupFeatures}}).
#' @param \dots One or more lists with parameter sets (see below) (for
#'   \code{optimizeFeatureFinding} and \code{optimizeFeatureGrouping}).
#'   Alternatively, named arguments that set (and possibly override) the
#'   parameters that should be returned from \code{generateFeatureOptPSet} or
#'   \code{generateFGroupsOptPSet}.
#' @param templateParams Template parameter set (see below).
#' @param paramRanges A list with vectors containing absolute parameter ranges
#'   (minimum/maximum) that constrain numeric parameters choosen during
#'   experiments. See the \code{\link{getDefFeaturesOptParamRanges}} and
#'   \code{\link{getDefFGroupsOptParamRanges}} functions for defaults. Values
#'   should be \code{Inf} when no limit should be used.
#' @param maxIterations Maximum number of iterations that may be performed to
#'   find optimimum values. Used to restrict neededless long
#'   optimization procedures. In IPO this was fixed to \samp{50}.
#' @param maxModelDeviation See the \verb{Potential suboptimal results by
#'   optimization model} section below.
#'
#' @section Parameter sets: Which parameters should be optimized is determined
#'   by a \emph{parameter set}. A set is defined by a named \code{list}
#'   containing the minimum and maximum starting range for each parameter that
#'   should be tested. For instance, the set \code{list(chromFWHM = c(5, 10),
#'   mzPPM = c(5, 15))} specifies that the \code{chromFWHM} and \code{mzPPM}
#'   parameters (used by OpenMS feature finding) should be optimized within a
#'   range of \samp{5}-\samp{10} and \samp{5}-\samp{15}, respectively. Note that
#'   this range may be increased or decreased after a DoE iteration in order to
#'   find a better optimum. The absolute limits are controlled by the
#'   \code{paramRanges} function argument.
#'
#'   Multiple parameter sets may be specified (\emph{i.e.} through the \dots
#'   function argument). In this situation, the optimization algorithm is
#'   repeated for each set, and the final optimum is determined from the parameter
#'   set with the best response. The \code{templateParams} function argument may
#'   be useful in this case to define a template for each parameter set. Actual
#'   parameter sets are then constructed by joining each parameter set with the
#'   set specified for \code{templateParams}. When a parameter is defined in
#'   both a regular and template set, the parameter in the regular set takes
#'   precedence.
#'
#'   Parameters that should not be optimized but still need to be set for the
#'   feature finding/grouping functions should also be defined in a (template)
#'   parameter set. Which parameters should be optimized is determined whether
#'   its value is specified as a vector range or a single fixed value.
#'   For instance, when a set is defined as \code{list(chromFWHM = c(5, 10),
#'   mzPPM = 5)}, only the \code{chromFWHM} parameter is optimized, whereas
#'   \code{mzPPM} is kept constant at \samp{5}.
#'
#'   Using multiple parameter sets with differing fixed values allows
#'   optimization of qualitative values (see examples below).
#'   
#'   The parameters specified in parameter sets are directly passed through
#'   the \code{\link{findFeatures}} or \code{\link{groupFeatures}} functions.
#'   Hence, grouping and retention time alignment parameters used by XCMS should
#'   (still) be set through the \code{groupArgs} and \code{retcorArgs}
#'   parameters.
#'   
#'   \strong{NOTE:} For XCMS3, which normally uses parameter classes for
#'   settings its options, the parameters must be defined in a named list like
#'   any other algorithm. The set parameters are then used passed to the
#'   constructor of the right parameter class object (e.g.
#'   \code{\link{CentWaveParam}}, \code{\link{ObiwarpParam}}). For
#'   grouping/alignment sets, these parameters need to be specified in nested
#'   lists called \code{groupParams} and \code{retAlignParams}, respectively
#'   (similar to \code{groupArgs}/\code{retcorArgs} for
#'   \code{algorithm="xcms"}). Finally, the underlying XCMS method to be used
#'   should be defined in the parameter set (\emph{i.e.} by setting the
#'   \code{method} field for feature parameter sets and the \code{groupMethod}
#'   and \code{retAlignMethod} for grouping/aligning parameter sets). See the
#'   examples below for more details.
#'
#'   \strong{NOTE:} Similar to IPO, the \code{peakwidth} and \code{prefilter}
#'   parameters for XCMS feature finding should be split in two different
#'   values:
#'   \itemize{
#'     \item The minimum and maximum ranges for \code{peakwidth} are optimized by
#'   setting \code{min_peakwidth} and \code{max_peakwidth}, respectively.
#'     \item The \code{k} and \code{I} parameters contained in \code{prefilter}
#'     are split in \code{prefilter} and \code{value_of_prefilter},
#'     respectively.
#'   }
#'
#' @section Functions: The \code{optimizeFeatureFinding} and
#'   \code{optimizeFeatureGrouping} are the functions to be used to optimize
#'   parameters for feature finding and grouping, respectively. These functions
#'   are analogous to \code{\link[IPO]{optimizeXcmsSet}} and
#'   \code{\link[IPO]{optimizeRetGroup}} from \pkg{IPO}.
#'
#'   The \code{generateFeatureOptPSet} and \code{generateFGroupsOptPSet} functions
#'   may be used to generate a parameter set for feature finding and grouping,
#'   respectively. Some algorithm dependent default parameter optimization ranges
#'   will be returned. These functions are analogous to
#'   \code{\link[IPO]{getDefaultXcmsSetStartingParams}} and
#'   \code{\link[IPO]{getDefaultRetGroupStartingParams}} from \pkg{IPO}. However,
#'   unlike their IPO counterparts, these functions will not output default fixed
#'   values. The \code{generateFGroupsOptPSet} will only generate defaults for
#'   density grouping if \code{algorithm="xcms"}.
#'
#'   The \code{getDefFeaturesOptParamRanges} and
#'   \code{getDefFGroupsOptParamRanges} return the default absolute optimization
#'   parameter ranges for feature finding and grouping, respectively. These
#'   functions are useful if you want to set the \code{paramRanges} function
#'   argument.
#'
#' @section Potential suboptimal results by optimization model: After each
#'   experiment iteration an optimimum parameter set is found by generating a
#'   model containing the tested parameters and their responses. Sometimes the
#'   actual response from the parameters derived from the model is actually
#'   signficantly lower than expected. When the response is lower than the
#'   maximum reponse found during the experiment, the parameters belonging to
#'   this experimental maximum may be choosen instead. The
#'   \code{maxModelDeviation} argument sets the maximum deviation in response
#'   between the modelled and experimental maxima. The value is relative:
#'   \samp{0} means that experimental values will always be favored when leading
#'   to improved responses, whereas \code{1} will effectively disable this
#'   procedure (and return to 'regular' IPO behaviour).
#'
#' @return The \code{optimizeFeatureFinding} and \code{optimizeFeatureGrouping}
#'   return their results in a \code{\link{optimizationResult}} object.
#'
#' @section Source: The code and methodology is a direct adaptation from the
#'   \href{https://github.com/rietho/IPO}{IPO R package}.
#'
#' @references \insertRef{Libiseller2015}{patRoon}
#'
#' @example inst/examples/optimization.R
#'
#' @name feature-optimization
NULL

#' Generation of MS Peak Lists
#'
#' Functionality to generate MS peak lists.
#'
#' Formula calculation and identification tools rely on mass spectra that belong
#' to features of interest. For processing, MS (and MS/MS) spectra are typically
#' reduced to a table with a column containing measured \emph{m/z} values and a
#' column containing their intensities. These 'MS peak lists' can then be used
#' for \link[=formula-generation]{formula generation} and
#' \link[=compound-generation]{compound generation}.
#'
#' MS and MS/MS peak lists are first generated for all features (or a subset, if
#' the \code{topMost} argument is set). During this step multiple spectra over
#' the feature elution profile are averaged. Subsequently, peak lists will be
#' generated for each feature group by averaging peak lists of the features
#' within the group. Functionality that uses peak lists will either use data
#' from individual features or from group averaged peak lists. For instance, the
#' former may be used by formulae calculation, while compound identification and
#' plotting functionality typically uses group averaged peak lists.
#'
#' Several functions exist to automatically extract MS peak lists for feature
#' groups.
#'
#' @param fGroups The \code{\link{featureGroups}} object from which MS peak
#'   lists should be extracted.
#' @param topMost Only extract MS peak lists from a maximum of \code{topMost}
#'   analyses with highest intensity. If \code{NULL} all analyses will be used.
#' @param maxMSRtWindow Maximum chromatographic peak window used for spectrum
#'   averaging (in seconds, +/- retention time). If \code{NULL} all spectra from
#'   a feature will be taken into account. Lower to decrease processing time.
#' @param minMSIntensity,minMSMSIntensity Minimum intensity for peak lists
#'   obtained with DataAnalysis. Highly recommended to set \samp{>0} as DA tends
#'   to report many very low intensity peaks.
#' @param avgFeatParams,avgFGroupParams A \code{list} with parameters used for
#'   averaging of peak lists for individual features and feature groups,
#'   respectively (see below).
#' @param \dots For \code{generateMSPeakLists}: Any parameters to be passed to
#'   the selected MS peak lists generation algorithm.
#'
#'   For \code{getDefAvgPListParams}: Optional named arguments that override
#'   defaults.
#'
#' @template dasaveclose-args
#'
#' @section Peak list averaging parameters: The parameters set used for
#'   averaging peak lists are set by the \code{avgFeatParams} and
#'   \code{avgFGroupParams} arguments. This should be a named \code{list} with
#'   the following values: \itemize{
#'
#'   \item \code{clusterMzWindow} \emph{m/z} window (in Da) used for clustering
#'   \emph{m/z} values when spectra are averaged. For \code{method="hclust"}
#'   this corresponds to the cluster height, while for \code{method="distance"}
#'   this value is used to find nearby masses (+/- window).  Too small windows
#'   will prevent clustering \emph{m/z} values (thus erroneously treating equal
#'   masses along spectra as different), whereas too big windows may cluster
#'   unrelated \emph{m/z} values from different or even the same spectrum
#'   together.
#'
#'   \item \code{topMost} Only retain this maximum number of MS peaks when
#'   generating averaged spectra. Lowering this number may exclude more
#'   irrelevant (noisy) MS peaks and decrease processing time, whereas higher
#'   values may avoid excluding lower intense MS peaks that may still be of
#'   interest.
#'
#'   \item \code{minIntensityPre} MS peaks with intensities below this value
#'   will be removed (applied prior to selection by \code{topMost}) before
#'   averaging.
#'
#'   \item \code{minIntensityPost} MS peaks with intensities below this value
#'   will be removed after averaging.
#'
#'   \item \code{avgFun} Function that is used to calculate average \emph{m/z}
#'   values.
#'
#'   \item \code{method} Method used for producing averaged MS spectra. Valid
#'   values are \code{"hclust"}, used for hierarchical clustering (using the
#'   \pkg{\link{fastcluster}} package), and \code{"distance"}, to use the
#'   between peak distance. The latter method may reduces processing time and
#'   memory requirements, at the potential cost of reduced accuracy.
#'
#'   \item \code{pruneMissingPrecursorMS} For MS data only: if \code{TRUE} then
#'   peak lists without a precursor peak are removed. Note that even when this
#'   is set to \code{FALSE}, functionality that relies on MS (not MS/MS) peak
#'   lists (\emph{e.g.} formulae calulcation) will still skip calculation if a
#'   precursor is not found.
#'
#'   \item \code{retainPrecursorMSMS} For MS/MS data only: if \code{TRUE} then
#'   always retain the precursor mass peak even if is not amongst the
#'   \code{topMost} peaks. Note that MS precursor mass peaks are always kept.
#'   Furthermore, note that precursor peaks in both MS and MS/MS data may still
#'   be removed by intensity thresholds (this is unlike the
#'   \code{\link[=filter,MSPeakLists-method]{filter}} method function).
#'
#'   }
#'
#'   Note that when Bruker algorithms are used these parameters only control
#'   generation of feature groups averaged peak lists: how peak lists for
#'   features are generated is controlled by DataAnalysis.
#'
#'   The \code{getDefAvgPListParams} function can be used to generate a default
#'   parameter list. The current defaults are:
#'
#' @eval paste0("@@section Peak list averaging parameters:",
#'   getDefAvgPListParamsRD())
#'
#' @return A \code{\link{MSPeakLists}} object that can be used for formulae
#'   calculation and compound identification.
#'
#' @section Source: Averaging of mass spectra algorithms used by are based on
#'   the \href{https://github.com/zeehio/msProcess}{msProcess} R package (now
#'   archived on CRAN).
#'
#' @references \addCitations{mzR}{1} \cr\cr
#'
#'   \addCitations{mzR}{2} \cr\cr
#'
#'   \addCitations{mzR}{3} \cr\cr
#'
#'   \addCitations{mzR}{4} \cr\cr
#'
#'   \addCitations{mzR}{5} \cr\cr
#'
#'   \addCitations{fastcluster}{1}
#'
#' @seealso \code{\link{MSPeakLists-class}}
#' @name MSPeakLists-generation
NULL

#' Automatic chemical formula generation
#'
#' Functionality to automatically calculate chemical formulae for all feature
#' groups.
#'
#' Several algorithms are provided to automatically generate formulae for given
#' feature groups. All tools use the accurate mass of a feature to
#' back-calculate candidate formulae. Depending on the algorithm and data
#' availability, other data such as isotopic pattern and MS/MS fragments may be
#' used to further improve formula assignment and ranking.
#'
#' When DataAnalysis is used for formula generation or
#' \code{calculateFeatures=TRUE} formulae are first calculated for each feature.
#' The results are then combined for final assignment of candidate formulae for
#' each feature group. If a formula was found in multiple features within the
#' group, the reported scorings and mass errors are averaged and other numeric
#' values are those from the feature in the analysis of the \code{"analysis"}
#' column. The calculation of formulae on 'feature level' might result in a
#' more thorough formula search and better removal of outliers (controlled by
#' \code{featThreshold} argument). In contrast, when calculations occur on
#' 'feature group level' (\emph{i.e.} \code{calculateFeatures=FALSE}), formulae
#' are directly assigned to each feature group (by using group averaged peak MS
#' lists), which significantly reduces processing time is, especially with many
#' analyses. Note that in both situations subsequent algorithms that use formula
#' data (\emph{e.g.} \code{\link{addFormulaScoring}} and \link{reporting}
#' functions) only use formula data that was eventually assigned to feature
#' groups. Furthermore, please note that calculation of formulae with
#' DataAnalysis always occurs on 'feature level'.
#'
#' @param fGroups \code{\link{featureGroups}} object for which formulae should
#'   be generated. This should be the same or a subset of the object that was
#'   used to create the specified \code{MSPeakLists} (only relevant for
#'   algorithms using \code{MSPeakLists}). In the case of a subset only the
#'   remaining feature groups in the subset are considered.
#' @param MSPeakLists An \code{\link{MSPeakLists}} object that was generated for
#'   the supplied \code{fGroups}.
#' @param MSMode Whether formulae should be generated only from MS data
#'   (\code{"ms"}), MS/MS data (\code{"msms"}) or both (\code{"both"}).
#'
#'   For \command{GenForm} selecting \code{"both"} will fall back to formula
#'   calculation with only MS data in case no MS/MS data is available.
#'
#'   For calulation with Bruker DataAnalysis selecting \code{"both"} will
#'   calculate formulae from MS data \emph{and} MS/MS data and combines the
#'   results (duplicated formulae are removed). This is useful when poor MS/MS
#'   data would exclude proper candidates.
#' @param calculateFeatures If \code{TRUE} fomulae are first calculated for all
#'   features prior to feature group assignment (see details).
#' @param featThreshold If \code{calculateFeatures=TRUE}: minimum presence
#'   (\samp{0-1}) of a formula in all features before it is considered as a
#'   candidate for a feature group. For instance, \code{featThreshold} dictates
#'   that a formula should be present in at least 75% of the features inside a
#'   feature group.
#' @param extraOpts An optional character vector with any other commandline
#'   options that will be passed to \command{GenForm} or \command{SIRIUS}. See
#'   the \verb{GenForm options} section/SIRIUS manual for all available
#'   commandline options.
#' @param topMost Only keep this number of candidates (per feature group) with
#'   highest score. For \command{SIRIUS}: Sets the \option{--candidates}
#'   commandline option.
#' 
#' @templateVar genForm TRUE
#' @template form-args
#'
#' @template adduct-arg
#'
#' @section Scorings: Each algorithm implements their own scoring system. Their
#'   names have been harmonized where possible. An overview is obtained with the
#'   \code{formulaScorings} function:
#'   \Sexpr[results=rd,echo=FALSE,stage=build]{patRoon:::tabularRD(patRoon::formulaScorings())}
#'
#' @templateVar what \code{generateFormulasGenForm} and \code{generateFormulasSIRIUS}
#' @template uses-multiProc
#'
#' @return A \code{\link{formulas}} object containing all generated formulae.
#'
#' @seealso \code{\link{formulas-class}}. The
#'   \href{https://www.researchgate.net/publication/307964728_MOLGEN-MSMS_Software_User_Manual}{GenForm
#'    manual} (also known as MOLGEN-MSMS).
#' @name formula-generation
NULL

#' Automatic compound identification
#'
#' Functionality to automatically identify chemical compounds from feature
#' groups.
#'
#' Several algorithms are provided to automatically identify compounds for given
#' feature groups. To this end, each measured masses for all feature groups are
#' searched within online database(s) (\emph{e.g.}
#' \href{https://pubchem.ncbi.nlm.nih.gov/}{PubChem}) to retrieve a list of
#' potential candidate chemical compounds. Depending on the algorithm and its
#' parameters, further scoring of candidates is then performed using, for
#' instance, matching of measured and theoretical isotopic patterns, presence
#' within other data sources such as patent databases and similarity of measured
#' and in-silico predicted MS/MS fragments. Note that this process is often
#' quite time consuming, especially for large feature group sets. Therefore,
#' this is often one of the last steps within the workflow and not performed
#' before feature groups have been prioritized.
#'
#' @param fGroups \code{\link{featureGroups}} object for which compounds should
#'   be identified. This should be the same or a subset of the object that was
#'   used to create the specified \code{MSPeakLists}. In the case of a subset
#'   only the remaining feature groups in the subset are considered.
#' @param MSPeakLists A \code{\link{MSPeakLists}} object that was generated for
#'   the supplied \code{fGroups}.
#' @param errorRetries Maximum number of retries after an error occurred. This
#'   may be useful to handle e.g. connection errors.
#' @param topMost Only keep this number of candidates (per feature group) with
#'   highest score. Set to \code{NULL} to always keep all candidates,
#'   however, please note that this may result in significant usage of CPU/RAM
#'   resources for large numbers of candidates.
#'
#' @param extraOpts For \command{MetFrag}: A named \code{list} containing
#'   further settings to be passed to \code{\link[metfRag]{run.metfrag}}. See
#'   the \href{http://ipb-halle.github.io/MetFrag/projects/metfragr/}{MetFragR}
#'   and \href{http://ipb-halle.github.io/MetFrag/projects/metfragcl/}{MetFrag
#'   CL} homepages for all available options.
#'
#'   For \command{SIRIUS}: a \code{character} vector with any extra commandline
#'   parameters for formula prediction. See the SIRIUS manual for more details.
#'
#'   Set to \code{NULL} to ignore.
#'
#' @template adduct-arg
#'
#' @section Scorings: Each algorithm implements their own scoring system. Their
#'   names have been simplified and harmonized where possible and are used for
#'   reporting and in the case \command{MetFrag} is used to specify how
#'   compounds should be scored (\code{scoreTypes} argument). The
#'   \code{compoundScorings} function can be used to get an overview of both the
#'   algorithm specific and generic scoring names. For instance, the table below
#'   shows all scorings for \command{MetFrag}: (some columns are omitted)
#'
#' @eval paste("@@section Scorings:",
#'   patRoon:::tabularRD(patRoon::compoundScorings("metfrag")[, 1:3]))
#'
#' @section Scorings: In addition, the \code{compoundScorings} function is also
#'   useful to programatically generate  a set of scorings to be used by
#'   \command{MetFrag}. For instance, the following can be given to the
#'   \code{scoreTypes} argument to use all default scorings for PubChem:
#'   \code{compoundScorings("metfrag", "pubchem", onlyDefault=TRUE)$name}.
#'
#'   For all \command{MetFrag} scoring types refer to the \verb{Candidate
#'   Scores} section on the
#'   \href{http://ipb-halle.github.io/MetFrag/projects/metfragr/}{MetFragR
#'   homepage}.
#'
#' @templateVar what \code{generateCompoundsMetFrag} and \code{generateCompoundsSIRIUS}
#' @template uses-multiProc
#'
#' @seealso \code{\link{compounds-class}}
#' @name compound-generation
NULL

#' Grouping feature groups in components
#'
#' Functionality to automatically group related feature groups (\emph{e.g.}
#' isotopes, adducts and homologues) to assist and simplify compound annotation.
#'
#' Several algorithms are provided to group feature groups that are related in
#' some (chemical) way to each other. These components generally include
#' adducts, isotopes, in-source fragments and homologues. The linking of this
#' data is generally useful to provide more information for compound annotation
#' and reduce the data size and thus complexity.
#'
#' @param fGroups \code{\link{featureGroups}} object for which components should
#'   be generated.
#' @param ionization Which ionization polarity was used to generate the data:
#'   should be \code{"positive"} or \code{"negative"}.
#' @param minSize The minimum size of a component. Smaller components than this
#'   size will be removed. For \code{RAMClustR}: sets the \code{minModuleSize}
#'   argument to \code{\link[RAMClustR]{ramclustR}}. See note below.
#' @param relMinReplicates Feature groups within a component are only kept when
#'   they contain data for at least this (relative) amount of replicate
#'   analyses. For instance, \samp{0.5} means that at least half of the
#'   replicates should contain data for a particular feature group in a
#'   component. In this calculation replicates that are fully absent within a
#'   component are not taken in to account. See note below.
#' @param extraOpts Named character vector with extra arguments directly passed
#'   to \code{\link{homol.search}} (\code{generateComponentsNontarget}) or
#'   \code{\link[CAMERA:annotate-methods]{CAMERA::annotate}}
#'   (\code{generateComponentsCAMERA}). Set to \code{NULL} to ignore.
#' @param rtDev,absMzDev,relMzDev Maximum deviation for retention time or
#'   absolute/relative \emph{m/z}.
#'
#'   For \code{generateComponentsRAMClustR}: Sets the \code{mzabs.error} and
#'   \code{ppm.error} arguments to \code{\link[RAMClustR]{do.findmain}}.
#'
#'   For \code{generateComponentsNontarget}: Sets the \code{rttol} and
#'   \code{mztol} arguments of \code{\link{homol.search}}.
#'
#' @return A \code{\link{components}} (derived) object containing all generated
#'   components.
#'
#' @note For \code{generateComponentsCAMERA} and
#'   \code{generateComponentsRAMClustR}: the \code{minSize} and
#'   \code{relMinReplicates} arguments provide additional filtering
#'   functionality not provided by \pkg{CAMERA} or \pkg{RAMClustR} (except
#'   \code{minSize}). Note that these filters are enabled by default, hence,
#'   final results may be different than what CAMERA/RAMClustR normally would
#'   return.
#'
#' @name component-generation
NULL

#' Target and suspect screening
#'
#' Utilities to screen for analytes with known or suspected identity.
#'
#' Besides 'full non-target analysis', where compounds may be identified with
#' little to no prior knowledge, a common strategy is to screen for compounds
#' with known or suspected identity. This may be a generally favorable approach
#' if possible, as it can significantly reduce the load on data interpretation.
#'
#' @note Both \code{screenSuspects} and \code{importFeatureGroupsBrukerTASQ} may
#'   use the suspect names to base file names used for reporting, logging etc.
#'   Therefore, it is important that these are file-compatible names.
#'
#' @name suspect-screening
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

#' Interactive GUI utilities
#'
#' Interactive utilities using \pkg{\link{shiny}} to provide a graphical user
#' interface (GUI).
#'
#' @param fGroups A \code{\link{featureGroups}} object.
#'
#' @name GUI-utils
NULL

