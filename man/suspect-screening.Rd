% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/main.R,
%   R/feature_groups-screening.R, R/feature_groups-screening-set.R,
%   R/feature_groups-tasq.R, R/features-tasq.R
\docType{class}
\name{screenSuspects}
\alias{screenSuspects}
\alias{suspect-screening}
\alias{screenSuspects,featureGroups-method}
\alias{screenSuspects,featureGroupsScreening-method}
\alias{numericIDLevel}
\alias{genIDLevelRulesFile}
\alias{screenSuspects,featureGroupsSet-method}
\alias{screenSuspects,featureGroupsScreeningSet-method}
\alias{featureGroupsBrukerTASQ-class}
\alias{featureGroupsBrukerTASQ}
\alias{featuresBrukerTASQ-class}
\alias{featuresBrukerTASQ}
\title{Target and suspect screening}
\usage{
screenSuspects(
  fGroups,
  suspects,
  rtWindow = 12,
  mzWindow = 0.005,
  adduct = NULL,
  skipInvalid = TRUE,
  prefCalcChemProps = TRUE,
  neutralChemProps = FALSE,
  onlyHits = FALSE,
  ...
)

\S4method{screenSuspects}{featureGroups}(
  fGroups,
  suspects,
  rtWindow,
  mzWindow,
  adduct,
  skipInvalid,
  prefCalcChemProps,
  neutralChemProps,
  onlyHits
)

\S4method{screenSuspects}{featureGroupsScreening}(
  fGroups,
  suspects,
  rtWindow,
  mzWindow,
  adduct,
  skipInvalid,
  onlyHits,
  amend = FALSE
)

numericIDLevel(level)

genIDLevelRulesFile(out, inLevels = NULL, exLevels = NULL)

\S4method{screenSuspects}{featureGroupsSet}(
  fGroups,
  suspects,
  rtWindow,
  mzWindow,
  adduct,
  skipInvalid,
  prefCalcChemProps,
  neutralChemProps,
  onlyHits
)

\S4method{screenSuspects}{featureGroupsScreeningSet}(
  fGroups,
  suspects,
  rtWindow,
  mzWindow,
  adduct,
  skipInvalid,
  prefCalcChemProps,
  neutralChemProps,
  onlyHits,
  amend = FALSE
)
}
\arguments{
\item{fGroups}{The \code{\link{featureGroups}} object that should be screened.}

\item{suspects}{A \code{data.frame} with suspect information. See the \verb{Suspect list format} section below.

  \setsWF Can also be a \code{list} with suspect lists to be used for each set (otherwise the same suspect lists is
  used for all sets). The \code{list} can be named with the sets names to mark which suspect list is to be used with
  which set (\emph{e.g.} \code{suspects=list(positive=suspsPos, negative=suspsNeg)}).}

\item{rtWindow, mzWindow}{The retention time window (in seconds) and \emph{m/z} window that will be used for matching
a suspect (+/- feature data).}

\item{adduct}{An \code{\link{adduct}} object (or something that can be converted to it with \code{\link{as.adduct}}).
Examples: \code{"[M-H]-"}, \code{"[M+Na]+"}. May be \code{NULL}, see \verb{Suspect list format} and \verb{Matching
of suspect masses} sections below.}

\item{skipInvalid}{If set to \code{TRUE} then suspects with invalid data (\emph{e.g.} missing names or other missing
data) will be ignored with a warning. Similarly, any suspects for which mass calculation failed (when no \code{mz}
column is present in the suspect list), for instance, due to invalid \code{SMILES}, will be ignored with a warning.}

\item{prefCalcChemProps}{If \code{TRUE} then calculated chemical properties such as the formula and
\acronym{InChIKey} are preferred over what is already present in the suspect list. For efficiency reasons it is
recommended to set this to \code{TRUE}. See the \verb{Validating and calculating chemical properties} section for
more details.}

\item{neutralChemProps}{If \code{TRUE} then the neutral form of the molecule is considered to calculate
\acronym{SMILES}, formulae etc. Enabling this may improve feature matching when considering common adducts
(\emph{e.g.} \code{[M+H]+}, \code{[M-H]-}). See the \verb{Validating and calculating chemical properties} section
for more details.}

\item{onlyHits}{If \code{TRUE} then all feature groups not matched by any of the suspects will be removed.}

\item{\dots}{Further arguments specified to the methods.}

\item{amend}{If \code{TRUE} then screening results will be \emph{amended} to the original object.}

\item{level}{The identification level to be converted.}

\item{out}{The file path to the target file.}

\item{inLevels, exLevels}{A \link[=regex]{regular expression} for the
identification levels to include or exclude, respectively. For instance,
\code{exLevels="4|5"} would exclude level 4 and 5 from the output file. Set
to \code{NULL} to ignore.}
}
\value{
\code{screenSuspects} returns a \code{\link{featureGroupsScreening}} object, which is a copy of the input
  \code{fGroups} object amended with additional screening information.
}
\description{
Utilities to screen for analytes with known or suspected identity.
}
\details{
Besides 'full non-target analysis', where compounds may be identified with little to no prior knowledge, a common
strategy is to screen for compounds with known or suspected identity. This may be a generally favorable approach if
possible, as it can significantly reduce the load on data interpretation.

\code{screenSuspects} is used to perform suspect screening. The input \code{\link{featureGroups}} object
  will be screened for suspects by \emph{m/z} values and optionally retention times. Afterwards, any feature groups
  not matched may be kept or removed, depending whether a full non-target analysis is desired.

\code{numericIDLevel} Extracts the numeric part of a given
  identification level (\emph{e.g.} \code{"3a"} becomes \samp{3}).

\code{genIDLevelRulesFile} Generates a template YAML file that is
  used to configure the rules for automatic estimation of identification
  levels. This file can then be used as input for
  \code{\link{annotateSuspects}}.
}
\note{
Both \code{screenSuspects} may use the suspect names to base file names used for reporting, logging etc.
  Therefore, it is important that these are file-compatible names. For this purpose, \code{screenSuspects} will
  automatically try to convert long, non-unique and/or otherwise incompatible suspect names.
}
\section{Sets workflows}{
 In a \link[=sets-workflow]{sets workflow}, \code{screenSuspects} performs suspect screening
  for each set separately, and the screening results are combined afterwards. The \code{sets} column in the
  \code{screenInfo} data marks in which sets the suspect hit was found.
}

\section{Suspect list format}{
 the \code{suspects} argument for \code{screenSuspects} should be a \code{data.frame}
  with the following mandatory and optional columns:

  \itemize{

  \item \code{name} The suspect name. Must be file-compatible. (\strong{mandatory})

  \item \code{rt} The retention time (in seconds) for the suspect. If specified the suspect will only be matched if
  its retention matches the experimental value (tolerance defined by the \code{rtWindow} argument).
  (\strong{optional})

  \item \code{neutralMass},\code{formula},\code{SMILES},\code{InChI} The neutral monoisotopic mass, chemical formula,
  SMILES or InChI for the suspect. (data from one of these columns are \strong{mandatory} in case no value from the
  \code{mz} column is available for a suspect)

  \item \code{mz} The ionized \emph{m/z} of the suspect. (\strong{mandatory} unless it can be calculated from one of
  the aforementioned columns)

  \item \code{adduct} A \code{character} that can be converted with \code{\link{as.adduct}}. Can be used to
  automatically calculate values for the \code{mz} column. (\strong{mandatory} unless data from the \code{mz} column
  is available, the \code{adduct} argument is set or \code{fGroups} has adduct annotations)

  \item \code{fragments_mz},\code{fragments_formula} One or more MS/MS fragments (specified as \emph{m/z} or
  formulae, respectively). Multiple values can be specified by separating them with a semicolon (\verb{;}). This data
  is used by \code{\link{annotateSuspects}} to report detected MS/MS fragments and calculate identification levels.
  (\strong{optional})

  }
}

\section{Matching of suspect masses}{
 How the mass of a suspect is matched with the mass of a feature depends on the
  available data: \itemize{

  \item If the suspect has data from the \code{mz} column of the suspect list, then this data is matched with the
  detected feature \emph{m/z}.

  \item Otherwise, if the suspect has data in the \code{adduct} column of the suspect list, this data is used to
  calculate its \code{mz} value, which is then used like above.

  \item In the last case, the neutral mass of the suspect is matched with the neutral mass of the feature. Hence,
  either the \code{adduct} argument needs to be specified, or the \code{featureGroups} input object must have adduct
  annotations.

  }
}

\section{Validating and calculating chemical properties}{
 Chemical properties such as \acronym{SMILES},
  \acronym{InChIKey} and formula in the suspect list are automatically validated and calculated if missing/invalid.

  The internal validation/calculation process performs the following steps: \itemize{

  \item Validation of \acronym{SMILES}, \acronym{InChI}, \acronym{InChIKey} and formula data (if present). Invalid
  entries will be set to \code{NA}.

  \item If \code{neutralChemProps=TRUE} then chemical data (\acronym{SMILES}, formulae etc.) is neutralized by
  (de-)protonation (using the \command{--neutralized} option of \command{OpenBabel}). An additional column
  \code{molNeutralized} is added to mark those molecules that were neutralized. Note that neutralization requires
  either \acronym{SMILES} or \acronym{InChI} data to be available.

  \item The \acronym{SMILES} and \acronym{InChI} data are used to calculate missing or invalid \acronym{SMILES},
  \acronym{InChI}, \acronym{InChIKey} and formula data. If \code{prefCalcChemProps=TRUE} then existing
  \acronym{InChIKey} and formula data is overwritten by calculated values whenever possible.

  \item The chemical formulae which were \emph{not} calculated are verified and normalized. This process may be time
  consuming, and is potentially largely avoided by setting \code{prefCalcChemProps=TRUE}.

  \item Neutral masses are calculated for missing values (\code{prefCalcChemProps=FALSE}) or whenever possible
  (\code{prefCalcChemProps=TRUE}).

  }

  Note that calculation of formulae for molecules that are isotopically labelled is currently only supported for
  deuterium (2H) elements.

  This functionality relies heavily on \href{https://github.com/openbabel/openbabel}{OpenBabel}, please make sure it
  is installed.
}

\references{
\insertRef{OBoyle2011}{patRoon}
}
\seealso{
\code{featureGroupsScreening}
}
