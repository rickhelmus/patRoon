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
#'   \item \code{patRoon.path.MetFragCL}: The complete file path to the MetFrag CL \file{jar} file that \emph{must} be
#'   set when using \code{\link{generateCompoundsMetFrag}}. Example: \code{"C:/MetFrag2.4.2-CL.jar"}.
#'
#'   \item \code{patRoon.path.MetFragCompTox}: The complete file path to the CompTox database \file{csv} file. See
#'   \code{\link{generateCompounds}} for more details.
#'
#'   \item \code{patRoon.path.MetFragPubChemLite}: The complete file path to the PubChem database \file{csv} file. See
#'   \code{\link{generateCompounds}} for more details.
#'
#'   \item \code{patRoon.path.SIRIUS}: The directory in which SIRIUS is installed. Unless the binaries can be located
#'   via the \option{PATH} environment variable, this \emph{must} be set when using \code{\link{generateFormulasSIRIUS}}
#'   or \code{\link{generateCompoundsSIRIUS}}. Example: \code{"C:/sirius-win64-3.5.1"}.
#'
#'   \item \code{patRoon.path.OpenMS}: The path in which the \command{OpenMS} binaries are installed. Usually the
#'   location is added to the \option{PATH} environment variable when OpenMS is installed, in which case this option can
#'   be left empty.
#'
#'   \item \code{patRoon.path.pngquant}: The path of the \command{pngquant} binary that is used when optimizing
#'   \file{.png} plots generated by \code{\link{reportHTMLLegacy}} (with \code{optimizePng} set to \code{TRUE}). If the
#'   binary can be located through the \option{PATH} environment variable this option can remain empty. Note that some
#'   of the functionality of \code{reportHTMLLegacy} only locates the binary through the \option{PATH} environment
#'   variable, hence, it is recommended to set up \option{PATH} instead.
#'
#'   \item \code{patRoon.path.obabel}: The path in which the \command{OpenBabel} binaries are installed. Usually the
#'   location is added to the \option{PATH} environment variable when OpenBabel is installed, in which case this option
#'   can be left empty.
#'
#'   \item \code{patRoon.path.BiotransFormer} The full file path to the \command{biotransformer} \file{.jar} command
#'   line utility. This needs to be set when \code{\link{generateTPsBioTransformer}} is used. For more details see
#'   \url{https://bitbucket.org/djoumbou/biotransformer/src/master}.
#'
#'   }
#'   
"_PACKAGE"

#' Properties of sample analyses
#'
#' Properties for the sample analyses used in the workflow and utilities to automatically generate this information.
#'
#' In \pkg{patRoon} a \emph{sample analysis}, or simply \emph{analysis}, refers to a single MS analysis file (sometimes
#' also called \emph{sample} or \emph{file}). The \emph{analysis information} summarizes several properties for the
#' analyses, and is used in various steps throughout the workflow, such as \code{\link{findFeatures}}, averaging
#' intensities of feature groups and blank subtraction. This information should be in a \code{data.frame}, with the
#' following columns:
#'
#' \itemize{
#'
#' \item \code{path} the full path to the directory of the analysis.
#'
#' \item \code{analysis} the file name \strong{without} extension. Must be \strong{unique}, even if the \code{path} is
#' different.
#'
#' \item \code{group} name of \emph{replicate group}. A replicate group is used to group analyses together that are
#' replicates of each other. Thus, the \code{group} column for all  analyses considered to be belonging to the same
#' replicate group should have an equal (but unique) value. Used for \emph{e.g.} averaging and
#' \code{\link[=filter,featureGroups-method]{filter}}.
#'
#' \item \code{blank} all analyses within this replicate group are used by the \code{featureGroups} method of
#' \code{\link[=filter,featureGroups-method]{filter}} for blank subtraction. Multiple entries can be entered by
#' separation with a comma.
#'
#' \item \code{conc} a numeric value specifying the 'concentration' for the analysis. This can be actually any kind of
#' numeric value such as exposure time, dilution factor or anything else which may be used to form a linear
#' relationship.
#'
#' \item \code{norm_conc} a numeric value specifying the \emph{normalization concentration} for the analysis. See the
#' \verb{Feature intensity normalization} section in the \link[=featureGroups-class]{featureGroups documentation}) for
#' more details.
#'
#' }
#'
#' Most workflows steps work with \file{mzXML} and \file{mzML} file formats. However, some algorithms only support
#' support one format (\emph{e.g.} \code{\link{findFeaturesOpenMS}}, \code{\link{findFeaturesEnviPick}}) or a
#' proprietary format (\code{\link{findFeaturesBruker}}). To mix such algorithms in the same workflow, the analyses
#' should be present in all required formats within the \emph{same} directory as specified by the \code{path} column.
#'
#' Each analysis should only be specified \emph{once} in the analysis information, even if multiple file formats are
#' available. The \code{path} and \code{analysis} columns are internally used by \pkg{patRoon} to automatically find the
#' path of analysis files with the required format.
#'
#' The \code{group} column is \emph{mandatory} and needs to be non-empty for each analysis. The \code{blank} column
#' should also be present, however, this may be empty (\code{""}) for analyses where no blank subtraction should occur.
#' The \code{conc} column is only required when obtaining regression information is desired with the
#' \code{\link[=as.data.table,featureGroups-method]{as.data.table}} method. Similarly, the \code{norm_conc} is only
#' necessary for the \code{\link{normInts}} method.
#'
#' @name analysis-information
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

#' Target and suspect screening
#'
#' Utilities to screen for analytes with known or suspected identity.
#'
#' Besides 'full non-target analysis', where compounds may be identified with little to no prior knowledge, a common
#' strategy is to screen for compounds with known or suspected identity. This may be a generally favorable approach if
#' possible, as it can significantly reduce the load on data interpretation.
#'
#' @section Sets workflows: In a \link[=sets-workflow]{sets workflow}, \code{screenSuspects} performs suspect screening
#'   for each set separately, and the screening results are combined afterwards. The \code{sets} column in the
#'   \code{screenInfo} data marks in which sets the suspect hit was found.
#'
#' @note Both \code{screenSuspects} may use the suspect names to base file names used for reporting, logging etc.
#'   Therefore, it is important that these are file-compatible names. For this purpose, \code{screenSuspects} will
#'   automatically try to convert long, non-unique and/or otherwise incompatible suspect names.
#'
#' @name suspect-screening
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
#' @note The \code{topMost} and \code{topMostByRGroup} EIC parameters (\code{\link{EICParams}}) are ignored.
#'
#' @references \insertRef{Chetnik2020}{patRoon}
#'
#' @name check-GUI
NULL

