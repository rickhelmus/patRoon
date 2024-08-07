% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/features.R, R/features-set.R,
%   R/features-bruker.R, R/features-envipick.R, R/features-kpic2.R,
%   R/features-openms.R, R/features-safd.R, R/features-sirius.R,
%   R/features-xcms.R, R/features-xcms3.R
\docType{class}
\name{features-class}
\alias{features-class}
\alias{features}
\alias{length,features-method}
\alias{show,features-method}
\alias{featureTable,features-method}
\alias{analysisInfo,features-method}
\alias{analyses,features-method}
\alias{replicateGroups,features-method}
\alias{as.data.table,features-method}
\alias{filter,features-method}
\alias{[,features,ANY,missing,missing-method}
\alias{[[,features,ANY,missing-method}
\alias{$,features-method}
\alias{delete,features-method}
\alias{calculatePeakQualities,features-method}
\alias{getTICs,features-method}
\alias{getBPCs,features-method}
\alias{plotTICs,features-method}
\alias{plotBPCs,features-method}
\alias{featuresSet-class}
\alias{featuresSet}
\alias{sets,featuresSet-method}
\alias{show,featuresSet-method}
\alias{as.data.table,featuresSet-method}
\alias{[,featuresSet,ANY,missing,missing-method}
\alias{filter,featuresSet-method}
\alias{featuresUnset-class}
\alias{featuresUnset}
\alias{unset,featuresSet-method}
\alias{featuresBruker-class}
\alias{featuresBruker}
\alias{featuresEnviPick-class}
\alias{featuresEnviPick}
\alias{featuresKPIC2-class}
\alias{featuresKPIC2}
\alias{delete,featuresKPIC2-method}
\alias{featuresOpenMS-class}
\alias{featuresOpenMS}
\alias{featuresSAFD-class}
\alias{featuresSAFD}
\alias{featuresSIRIUS-class}
\alias{featuresSIRIUS}
\alias{featuresXCMS-class}
\alias{featuresXCMS}
\alias{delete,featuresXCMS-method}
\alias{featuresXCMS3-class}
\alias{featuresXCMS3}
\alias{delete,featuresXCMS3-method}
\title{Base features class}
\usage{
\S4method{length}{features}(x)

\S4method{show}{features}(object)

\S4method{featureTable}{features}(obj)

\S4method{analysisInfo}{features}(obj)

\S4method{analyses}{features}(obj)

\S4method{replicateGroups}{features}(obj)

\S4method{as.data.table}{features}(x)

\S4method{filter}{features}(
  obj,
  absMinIntensity = NULL,
  relMinIntensity = NULL,
  retentionRange = NULL,
  mzRange = NULL,
  mzDefectRange = NULL,
  chromWidthRange = NULL,
  qualityRange = NULL,
  negate = FALSE
)

\S4method{[}{features,ANY,missing,missing}(x, i, j, ..., drop = TRUE)

\S4method{[[}{features,ANY,missing}(x, i)

\S4method{$}{features}(x, name)

\S4method{delete}{features}(obj, i = NULL, j = NULL, ...)

\S4method{calculatePeakQualities}{features}(obj, weights, flatnessFactor, parallel = TRUE)

\S4method{getTICs}{features}(obj, retentionRange = NULL, MSLevel = 1)

\S4method{getBPCs}{features}(obj, retentionRange = NULL, MSLevel = 1)

\S4method{plotTICs}{features}(
  obj,
  retentionRange = NULL,
  MSLevel = 1,
  retMin = FALSE,
  title = NULL,
  colourBy = c("none", "analyses", "rGroups"),
  showLegend = TRUE,
  xlim = NULL,
  ylim = NULL,
  ...
)

\S4method{plotBPCs}{features}(
  obj,
  retentionRange = NULL,
  MSLevel = 1,
  retMin = FALSE,
  title = NULL,
  colourBy = c("none", "analyses", "rGroups"),
  showLegend = TRUE,
  xlim = NULL,
  ylim = NULL,
  ...
)

\S4method{sets}{featuresSet}(obj)

\S4method{show}{featuresSet}(object)

\S4method{as.data.table}{featuresSet}(x)

\S4method{[}{featuresSet,ANY,missing,missing}(x, i, ..., sets = NULL, drop = TRUE)

\S4method{filter}{featuresSet}(obj, ..., negate = FALSE, sets = NULL)

\S4method{unset}{featuresSet}(obj, set)

\S4method{delete}{featuresKPIC2}(obj, i = NULL, j = NULL, ...)

\S4method{delete}{featuresXCMS}(obj, i = NULL, j = NULL, ...)

\S4method{delete}{featuresXCMS3}(obj, i = NULL, j = NULL, ...)
}
\arguments{
\item{obj, x, object}{\code{features} object to be accessed}

\item{absMinIntensity, relMinIntensity}{Minimum absolute/relative intensity for features to be kept. The relative
intensity is determined from the feature with highest intensity (within the same analysis). Set to \samp{0} or \code{NULL} to skip this step.}

\item{retentionRange, mzRange, mzDefectRange, chromWidthRange}{Range of retention time (in seconds), \emph{m/z}, mass
defect (defined as the decimal part of \emph{m/z} values) or chromatographic peak width (in seconds), respectively.
Features outside this range will be removed. Should be a numeric vector with length of two containing the min/max
values. The maximum can be \code{Inf} to specify no maximum range. Set to \code{NULL} to skip this step.}

\item{qualityRange}{Used to filter features by their peak qualities/scores
(see \code{\link{calculatePeakQualities}}). Should be a named \code{list} with min/max ranges for each
quality/score to be filtered (the \code{\link{featureQualityNames}} function can be used to obtain valid names).
Example: \code{qualityRange=list(ModalityScore=c(0.3, Inf),
SymmetryScore=c(0.5, Inf))}. Set to \code{NULL} to ignore.}

\item{negate}{If set to \code{TRUE} then filtering operations are performed in opposite manner.}

\item{i, j}{For \code{[}/\code{[[}: A numeric or character value which is used to select analyses by
their index or name, respectively (for the order/names see \code{analyses()}).\cr\cr For \code{[}: Can also be logical to perform logical selection
(similar to regular vectors). If missing all analyses are selected.\cr\cr For \code{[[}: should be a scalar value.\cr\cr For \code{delete}: The data to remove from. \code{i} are the
analyses as numeric index, logical or character, \code{j} the features as numeric index (row) of the feature. If either is
\code{NULL} then data for all is removed. \code{j} may also be a function: it will be called for each analysis,
with the feature table (a \code{data.table}), the analysis name and any other arguments passed as \code{\dots} to \code{delete}. The return value of this
function specifies the feature indices (rows) to be removed (specified as an \code{integer} or \code{logical} vector).}

\item{\dots}{For \code{delete}: passed to the function specified as \code{j}.

  For \code{plotTICs} and \code{plotBPCs}: further arguments passed to \code{\link[graphics]{plot}}.

  \setsPassedArgs1{features}}

\item{drop}{ignored.}

\item{name}{The analysis name (partially matched).}

\item{weights}{A named \code{numeric} vector that defines the weight for each score to calculate the
\verb{totalScore}. The names of the vector follow the score names. Unspecified weights are defaulted to \samp{1}.
Example: \code{weights=c(ApexBoundaryRatioScore=0.5, GaussianSimilarityScore=2)}.}

\item{flatnessFactor}{Passed to \pkg{MetaClean} as the \code{flatness.factor} argument to
\code{\link[MetaClean]{calculateJaggedness}} and \code{\link[MetaClean]{calculateModality}}.}

\item{parallel}{If set to \code{TRUE} then code is executed in parallel through the \CRANpkg{futures} package. Please
see the parallelization section in the handbook for more details.}

\item{MSLevel}{Integer vector with the ms levels (i.e., 1 for MS1 and 2 for MS2) 
to obtain TIC traces.}

\item{retMin}{Plot retention time in minutes (instead of seconds).}

\item{title}{Character string used for title of the plot. If \code{NULL} a title will be automatically generated.}

\item{colourBy}{Sets the automatic colour selection: "none" for a single colour or "analyses"/"rGroups" for a
distinct colour per analysis or analysis replicate group.}

\item{showLegend}{Plot a legend if TRUE.}

\item{xlim, ylim}{Sets the plot size limits used by
\code{\link[graphics]{plot}}. Set to \code{NULL} for automatic plot sizing.}

\item{sets}{\setsWF For \code{[} and \code{filter}: a \code{character} with name(s) of the sets to keep (or remove if
\code{negate=TRUE}).}

\item{set}{\setsWF The name of the set.}
}
\value{
\code{featureTable}: A \code{list} containing a
  \code{\link{data.table}} for each analysis with feature data

\code{analysisInfo}: A \code{data.frame} containing a column with
  analysis name (\code{analysis}), its path (\code{path}), and other columns
  such as replicate group name (\code{group}) and blank reference
  (\code{blank}).

\code{delete} returns the object for which the specified data was removed.

\code{calculatePeakQualities} returns a modified object amended with peak qualities and scores.
}
\description{
Holds information for all features present within a set of analysis.
}
\details{
This class provides a way to store intensity, retention times, \emph{m/z} and other data for all features in a set of
analyses. The class is \code{virtual} and derived objects are created by 'feature finders' such as
\code{findFeaturesOpenMS}, \code{findFeaturesXCMS} and \code{findFeaturesBruker}.
}
\section{Methods (by generic)}{
\itemize{
\item \code{length(features)}: Obtain total number of features.

\item \code{show(features)}: Shows summary information for this object.

\item \code{featureTable(features)}: Get table with feature information

\item \code{analysisInfo(features)}: Get analysis information

\item \code{analyses(features)}: returns a \code{character} vector with the names of the
analyses for which data is present in this object.

\item \code{replicateGroups(features)}: returns a \code{character} vector with the names of the
replicate groups for which data is present in this object.

\item \code{as.data.table(features)}: Returns all feature data in a table.

\item \code{filter(features)}: Performs common rule based filtering of features. Note
that this (and much more) functionality is also provided by the
\code{filter} method defined for \code{\link{featureGroups}}. However,
filtering a \code{features} object may be useful to avoid grouping large
amounts of features.

\item \code{x[i}: Subset on analyses.

\item \code{x[[i}: Extract a feature table for an analysis.

\item \code{$}: Extract a feature table for an analysis.

\item \code{delete(features)}: Completely deletes specified features.

\item \code{calculatePeakQualities(features)}: Calculates peak qualities for each feature. This uses
\href{https://github.com/KelseyChetnik/MetaClean/}{MetaClean} \R package to calculate the following metrics:
\verb{Apex-Boundary Ratio}, \verb{FWHM2Base}, \verb{Jaggedness}, \verb{Modality}, \verb{Symmetry}, \verb{Gaussian
Similarity}, \verb{Sharpness}, \verb{Triangle Peak Area Similarity Ratio} and \verb{Zig-Zag index}. Please see the
\pkg{MetaClean} publication (referenced below) for more details. For each metric, an additional score is calculated
by normalizing all feature values (unless the quality metric definition has a fixed range) and scale from \samp{0}
(worst) to \samp{1} (best). Then, a \verb{totalScore} for each feature is calculated by the (weighted) sum of all
score values.

\item \code{getTICs(features)}: Obtain the total ion chromatogram/s (TICs) of the analyses.

\item \code{getBPCs(features)}: Obtain the base peak chromatogram/s (BPCs) of the analyses.

\item \code{plotTICs(features)}: Plots the TICs of the analyses.

\item \code{plotBPCs(features)}: Plots the BPCs of the analyses.

}}
\section{Slots}{

\describe{
\item{\code{features}}{List of features per analysis file. Use the \code{featureTable} method for access.}

\item{\code{analysisInfo}}{Analysis group information. Use the \code{analysisInfo} method for access.}
}}

\note{
For \code{calculatePeakQualities}: sometimes \pkg{MetaClean} may return \code{NA} for the \verb{Gaussian
  Similarity} metric, in which case it will be set to \samp{0}.
}
\section{S4 class hierarchy}{
 \itemize{   \item{\code{\link{workflowStep}}}   \itemize{     \item{\strong{\code{\link{features}}}}     \itemize{       \item{\code{\link{featuresSet}}}       \item{\code{\link{featuresUnset}}}       \item{\code{\link{featuresFromFeatGroups}}}       \item{\code{\link{featuresConsensus}}}       \item{\code{\link{featuresBruker}}}       \item{\code{\link{featuresEnviPick}}}       \item{\code{\link{featuresKPIC2}}}       \item{\code{\link{featuresOpenMS}}}       \item{\code{\link{featuresSAFD}}}       \item{\code{\link{featuresSIRIUS}}}       \item{\code{\link{featuresBrukerTASQ}}}       \item{\code{\link{featuresXCMS}}}       \item{\code{\link{featuresXCMS3}}}     }   } }
}

\section{Sets workflows}{
 \setsWFClass{featuresSet}{features}

  \setsWFNewMethodsFeat{featuresUnset}{The adduct annotations for the selected set (\emph{e.g.} as passed to
  \code{makeSet}) are used to convert all feature masses to ionic \emph{m/z} values. }

  \setsWFChangedMethods{

  \item \code{filter} and the subset operator (\code{[}) have specific arguments to choose/filter by (feature
  presence in) sets. See the \code{sets} argument description.

  }
}

\references{
\insertRef{Chetnik2020}{patRoon}
}
\seealso{
\code{\link{findFeatures}}
}
\author{
Rick Helmus <\email{r.helmus@uva.nl}> and Ricardo Cunha <\email{cunha@iuta.de}> (\code{getTICs},
  \code{getBPCs}, \code{plotTICs} and \code{plotBPCs} functions)
}
