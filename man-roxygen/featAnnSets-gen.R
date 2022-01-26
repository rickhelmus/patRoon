#' @section Sets workflows: With a \link[=sets-workflow]{sets workflow}, annotation is first performed for each set.
#'   This is important, since the annotation algorithms typically cannot work with data from mixed ionization modes. The
#'   annotation results are then combined to generate a \emph{sets consensus}: \itemize{
#'
#'   \item The annotation tables for each feature group from the set specific data are combined. Rows with overlapping
#'   candidates (determined by the <%=UID%>) are merged.
#'
#'   \item Set specific data (\emph{e.g.} the ionic formula) is retained by renaming their columns with set specific
#'   names.
#'
#'   \item The MS/MS fragment annotations (\code{fragInfo} column) from each set are combined.
#'
#'   \item The scorings for each set are averaged to calculate overall scores.
#'
#'   \item The candidates are re-ranked based on their average ranking among the set data (if a candidate is absent in a
#'   set it is assigned the poorest rank in that set).
#'
#'   \item The coverage of each candidate among sets is calculated. Depending on the \code{setThreshold} and
#'   \code{setThresholdAnn} arguments, candidates with low abundance are removed.
#'
#'   }

