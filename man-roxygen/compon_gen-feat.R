#' @param relMinAdductAbundance The minimum relative abundance (\samp{0-1}) that an adduct should be assigned to
#'   features within the same feature group. See the \verb{Feature components} section for more details.
#' @param adductConflictsUsePref If set to \code{TRUE}, and not all adduct assigments to the features within a feature
#'   group are equal and at least one of those adducts is a preferential adduct (\code{prefAdducts} argument), then only
#'   the features with (the lowest ranked) preferential adduct are considered. In all other cases or when
#'   \code{adductConflictsUsePref=FALSE} only features with the most frequently assigned adduct is considered. See the
#'   \verb{Feature components} section for more details.
#' @param NMConflicts The strategies to employ when not all neutral masses within a component are equal. Valid options
#'   are: \code{"preferential"}, \code{"mostAbundant"} and \code{"mostIntense"}. Multiple strategies are possible, and
#'   will be executed in the given order until one succeeds. See the \verb{Feature components} section for more details.
#' @param prefAdducts A \code{character} vector with one or more \emph{preferential adducts}. See the \verb{Feature
#'   components} section for more details.
#'
#' @section Feature components: The returned components are based on so called \emph{feature components}. Unlike other
#'   algorithms, components are first made on a feature level (per analysis), instead of for complete feature groups. In
#'   the final step the feature components are converted to 'regular' components by employing a consensus approach with
#'   the following steps:
#'
#'   \enumerate{
#'
#'   \item If an adduct assigned to a feature only occurs as a minority compared to other adduct assigments within the
#'   same feature group, it is considered as an outlier and removed accordingly (controlled by the
#'   \code{relMinAdductAbundance} argument).
#'
#'   \item For features within a feature group, only keep their adduct assignment if it occurs as the most frequent or
#'   is preferential (controlled by \code{adductConflictsUsePref} and \code{prefAdducts} arguments).
#'
#'   \item Components are made by combining the feature groups for which at least one of their features are jointly
#'   present in the same feature component.
#'
#'   \item Conflicts of neutral mass assignments within a component (\emph{i.e.} not all are the same) are dealt with.
#'   Firstly, all feature groups with an unknown neutral mass are split in another component. Then, if conflicts still
#'   occur, the feature groups with similar neutral mass (determined by \code{absMzDev} argument) are grouped. Depending
#'   on the \code{NMConflicts} argument, the group with one or more preferential adduct(s) or that is the largest or
#'   most intense is selected, whereas others are removed from the component. In case multiple groups contain
#'   preferential adducts, and \samp{>1} preferential adducts are available, the group with the adduct that matches
#'   first in \code{prefAdducts} 'wins'. In case of ties, one of the next strategies in \code{NMConflicts} is tried.
#'
#'   \item If a feature group occurs in multiple components it will be removed completely.
#'
#'   \item the \code{minSize} filter is applied.
#'
#'   }
