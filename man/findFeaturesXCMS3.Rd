% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/features-xcms3.R
\name{findFeaturesXCMS3}
\alias{findFeaturesXCMS3}
\title{Find features using XCMS (new interface)}
\usage{
findFeaturesXCMS3(
  analysisInfo,
  param = xcms::CentWaveParam(),
  ...,
  verbose = TRUE
)
}
\arguments{
\item{analysisInfo}{A \code{data.frame} with \link[=analysis-information]{Analysis information}.}

\item{param}{The method parameters used by XCMS peak finding, see
\code{\link[xcms:findChromPeaks]{xcms::findChromPeaks}}}

\item{\dots}{Further parameters passed to \code{\link[xcms:findChromPeaks]{xcms::findChromPeaks}}.}

\item{verbose}{If set to \code{FALSE} then no text output is shown.}
}
\value{
An object of a class which is derived from \code{\link{features}}.
}
\description{
Uses the new \code{xcms3} interface from the \pkg{xcms} package to find features.
}
\details{
This function uses XCMS3 to automatically find features. This function is called when calling \code{findFeatures} with
  \code{algorithm="xcms3"}.

The file format of analyses must be \code{mzML} or \code{mzXML}.

The input MS data files need to be centroided. The \code{\link{convertMSFiles}} function can be used to
  centroid data.
}
\references{
\addCitations{xcms}
}
\seealso{
\code{\link{findFeatures}} for more details and other algorithms.
}
