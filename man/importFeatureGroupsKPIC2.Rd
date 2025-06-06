% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_groups-kpic2.R
\name{importFeatureGroupsKPIC2}
\alias{importFeatureGroupsKPIC2}
\title{Imports feature groups from KPIC2}
\usage{
importFeatureGroupsKPIC2(picsSetGrouped, analysisInfo)
}
\arguments{
\item{picsSetGrouped}{A grouped \code{PIC set} object (\emph{e.g.} as returned by
\code{\link[KPIC:PICset.group]{KPIC::PICset.group}}).}

\item{analysisInfo}{A \code{data.frame} with \link[=analysis-information]{Analysis information}.}
}
\value{
An object of a class which is derived from \code{\link{featureGroups}}.
}
\description{
Imports grouped features from an \pkg{KPIC} object.
}
\references{
\insertRef{Ji2017}{patRoon}
}
\seealso{
\code{\link{groupFeatures}}
}
