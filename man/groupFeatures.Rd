% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/feature_groups.R
\name{groupFeatures}
\alias{groupFeatures}
\alias{groupFeatures,features-method}
\alias{groupFeatures,data.frame-method}
\title{Grouping of features}
\usage{
groupFeatures(obj, algorithm, ...)

\S4method{groupFeatures}{features}(obj, algorithm, ..., verbose = TRUE)

\S4method{groupFeatures}{data.frame}(obj, algorithm, ..., verbose = TRUE)
}
\arguments{
\item{obj}{Either a \code{\link{features}} object to be grouped, or a \code{data.frame} with
\link[=analysis-information]{analysis info} to be passed to \code{groupFeaturesSIRIUS}}

\item{algorithm}{A \code{character} that specifies the algorithm to be used: either \code{"openms"}, \code{"xcms"},
\code{"xcms3"} or \code{"kpic2"} (\code{features method}), or \code{"sirius"} (\code{data.frame} method).}

\item{\dots}{Further parameters passed to the selected grouping algorithm.}

\item{verbose}{if \code{FALSE} then no text output will be shown.}
}
\value{
An object of a class which is derived from \code{\link{featureGroups}}.
}
\description{
Group equal features across analyses.
}
\details{
After \link[=findFeatures]{features have been found}, the next step is to align and group them across analyses. This
process is necessary to allow comparison of features between multiple analyses, which otherwise would be difficult
due to small deviations in retention and mass data. Thus, algorithms of 'feature groupers' are used to collect
features with similar retention and mass data. In addition, advanced retention time alignment algorithms exist to
enhance grouping of features even with relative large retention time deviations (\emph{e.g.} possibly observed from
analyses collected over a long period). Like \link{findFeatures}, various algorithms are supported which may have
many parameters that can be fine-tuned. This fine-tuning is likely to be necessary, since optimal settings often
depend on applied methodology and instrumentation.

\code{groupFeatures} is a generic function that will groupFeatures by one of the supported algorithms. The actual
  functionality is provided by algorithm specific functions such as \code{groupFeaturesOpenMS} and \code{groupFeaturesXCMS3}. While these
  functions may be called directly, \code{groupFeatures} provides a generic interface and is therefore usually preferred.

The \code{data.frame} method for \code{groupFeatures} is a special case that currently only supports the
  \code{"sirius"} algorithm.
}
\seealso{
The \code{\link{featureGroups}} output class and its methods and the algorithm specific functions:
  \code{\link{groupFeaturesOpenMS}}, \code{\link{groupFeaturesXCMS}}, \code{\link{groupFeaturesXCMS3}}, \code{\link{groupFeaturesKPIC2}}, \code{\link{groupFeaturesSIRIUS}}
}
