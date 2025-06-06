% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_groups-xcms3.R
\name{importFeatureGroupsXCMS3}
\alias{importFeatureGroupsXCMS3}
\title{Imports feature groups from XCMS (new interface)}
\usage{
importFeatureGroupsXCMS3(xdata, analysisInfo)
}
\arguments{
\item{xdata}{An \code{\link{XCMSnExp}} object.}

\item{analysisInfo}{A \code{data.frame} with \link[=analysis-information]{Analysis information}.}
}
\value{
An object of a class which is derived from \code{\link{featureGroups}}.
}
\description{
Imports grouped features from a \code{\link{XCMSnExp}} object from the \pkg{xcms} package.
}
\references{
\addCitations{xcms}
}
\seealso{
\code{\link{groupFeatures}}
}
