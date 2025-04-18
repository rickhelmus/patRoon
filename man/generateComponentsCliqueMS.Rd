% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/components-cliquems.R
\name{generateComponentsCliqueMS}
\alias{generateComponentsCliqueMS}
\alias{generateComponentsCliqueMS,featureGroups-method}
\alias{generateComponentsCliqueMS,featureGroupsSet-method}
\title{Componentization of adducts, isotopes etc. with cliqueMS}
\usage{
generateComponentsCliqueMS(fGroups, ...)

\S4method{generateComponentsCliqueMS}{featureGroups}(
  fGroups,
  ionization = NULL,
  maxCharge = 1,
  maxGrade = 2,
  ppm = 10,
  adductInfo = NULL,
  absMzDev = 0.005,
  minSize = 2,
  relMinAdductAbundance = 0.75,
  adductConflictsUsePref = TRUE,
  NMConflicts = c("preferential", "mostAbundant", "mostIntense"),
  prefAdducts = c("[M+H]+", "[M-H]-"),
  extraOptsCli = NULL,
  extraOptsIso = NULL,
  extraOptsAnn = NULL,
  parallel = TRUE
)

\S4method{generateComponentsCliqueMS}{featureGroupsSet}(fGroups, ionization = NULL, ...)
}
\arguments{
\item{fGroups}{\code{\link{featureGroups}} object for which components should be generated.}

\item{\dots}{\setsWF Further arguments passed to the non-sets workflow method.}

\item{ionization}{Which ionization polarity was used to generate the data: should be \code{"positive"}
  or \code{"negative"}. If the \code{featureGroups} object has adduct annotations, and \code{ionization=NULL}, the
  ionization will be detected automatically.

  \setsWF This parameter is not supported for sets workflows, as the ionization will always be detected
  automatically.}

\item{maxCharge, maxGrade, ppm}{Arguments passed to \code{\link[cliqueMS:getIsotopes]{cliqueMS::getIsotopes}} and/or
\code{\link[cliqueMS:getAnnotation]{cliqueMS::getAnnotation}}.}

\item{adductInfo}{Sets the \code{adinfo} argument to \code{\link[cliqueMS:getAnnotation]{cliqueMS::getAnnotation}}.
If \code{NULL} then the default adduct information from \pkg{cliqueMS} is used (\emph{i.e.} the
\code{positive.adinfo}/\code{negative.adinfo} package datasets).}

\item{absMzDev}{Maximum absolute \emph{m/z} deviation.}

\item{minSize}{The minimum size of a component. Smaller components than this size will be removed. See note below.}

\item{relMinAdductAbundance}{The minimum relative abundance (\samp{0-1}) that an adduct should be assigned to
features within the same feature group. See the \verb{Feature components} section for more details.}

\item{adductConflictsUsePref}{If set to \code{TRUE}, and not all adduct assigments to the features within a feature
group are equal and at least one of those adducts is a preferential adduct (\code{prefAdducts} argument), then only
the features with (the lowest ranked) preferential adduct are considered. In all other cases or when
\code{adductConflictsUsePref=FALSE} only features with the most frequently assigned adduct is considered. See the
\verb{Feature components} section for more details.}

\item{NMConflicts}{The strategies to employ when not all neutral masses within a component are equal. Valid options
are: \code{"preferential"}, \code{"mostAbundant"} and \code{"mostIntense"}. Multiple strategies are possible, and
will be executed in the given order until one succeeds. See the \verb{Feature components} section for more details.}

\item{prefAdducts}{A \code{character} vector with one or more \emph{preferential adducts}. See the \verb{Feature
components} section for more details.}

\item{extraOptsCli, extraOptsIso, extraOptsAnn}{Named \code{list} with further arguments to be passed to
\code{\link[cliqueMS:getCliques]{cliqueMS::getCliques}}, \code{\link[cliqueMS:getIsotopes]{cliqueMS::getIsotopes}}
and \code{\link[cliqueMS:getAnnotation]{cliqueMS::getAnnotation}}, respectively. Set to \code{NULL} to ignore.}

\item{parallel}{If set to \code{TRUE} then code is executed in parallel through the \CRANpkg{futures} package. Please
see the parallelization section in the handbook for more details.}
}
\value{
A \code{\link{componentsFeatures}} derived object.
}
\description{
Uses \href{https://github.com/osenan/cliqueMS}{cliqueMS} to generate components using the
\code{\link[cliqueMS:getCliques]{cliqueMS::getCliques}} function.
}
\details{
This function uses cliqueMS to generate components. This function is called when calling \code{generateComponents} with
  \code{algorithm="cliquems"}.

The grouping of features in each component ('clique') is based on high similarity of chromatographic elution
  profiles. All features in each component are then annotated with the
  \code{\link[cliqueMS:getIsotopes]{cliqueMS::getIsotopes}} and
  \code{\link[cliqueMS:getAnnotation]{cliqueMS::getAnnotation}} functions.
}
\section{Feature components}{
 The returned components are based on so called \emph{feature components}. Unlike other
  algorithms, components are first made on a feature level (per analysis), instead of for complete feature groups. In
  the final step the feature components are converted to 'regular' components by employing a consensus approach with
  the following steps:

  \enumerate{

  \item If an adduct assigned to a feature only occurs as a minority compared to other adduct assigments within the
  same feature group, it is considered as an outlier and removed accordingly (controlled by the
  \code{relMinAdductAbundance} argument).

  \item For features within a feature group, only keep their adduct assignment if it occurs as the most frequent or
  is preferential (controlled by \code{adductConflictsUsePref} and \code{prefAdducts} arguments).

  \item Components are made by combining the feature groups for which at least one of their features are jointly
  present in the same feature component.

  \item Conflicts of neutral mass assignments within a component (\emph{i.e.} not all are the same) are dealt with.
  Firstly, all feature groups with an unknown neutral mass are split in another component. Then, if conflicts still
  occur, the feature groups with similar neutral mass (determined by \code{absMzDev} argument) are grouped. Depending
  on the \code{NMConflicts} argument, the group with one or more preferential adduct(s) or that is the largest or
  most intense is selected, whereas others are removed from the component. In case multiple groups contain
  preferential adducts, and \samp{>1} preferential adducts are available, the group with the adduct that matches
  first in \code{prefAdducts} 'wins'. In case of ties, one of the next strategies in \code{NMConflicts} is tried.

  \item If a feature group occurs in multiple components it will be removed completely.

  \item the \code{minSize} filter is applied.

  }
}

\section{Sets workflows}{
 In a \link[=sets-workflow]{sets workflow} the componentization is first performed for each
  set independently. The resulting components are then all combined in a \code{\link{componentsSet}} object. Note that
  the components themselves are never merged. The components are renamed to include the set name from which they were
  generated (\emph{e.g.} \code{"CMP1"} becomes \code{"CMP1-positive"}).
}

\references{
\insertRef{Senan2019}{patRoon}
}
\seealso{
\code{\link{generateComponents}} for more details and other algorithms.
}
