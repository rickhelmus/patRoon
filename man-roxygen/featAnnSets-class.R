#' @param set \setsWF The name of the set.
#' @param sets \setsWF A \code{character} with name(s) of the sets to keep (or remove if \code{negate=TRUE}). Note: if
#'   \code{updateConsensus=FALSE} then the \code{setCoverage} column of the annotation results is not updated.
#' @param updateConsensus \setsWF If \code{TRUE} then the annonation consensus among set results is updated. See the
#'   \verb{Sets workflows} section for more details.
#' @param perSet,mirror \setsWF If \code{perSet=TRUE} then the set specific mass peaks are annotated separately.
#'   Furthermore, if \code{mirror=TRUE} (and there are two sets in the object) then a mirror plot is generated.
#' @param setThreshold,setThresholdAnn \setsWF Thresholds used to create the annotation set consensus. See
#'   \code{\link{<%=generator%>}}.
#'
#' @slot setThreshold,setThresholdAnn \setsWF A copy of the equally named arguments that were passed when this object
#'   was created by \code{\link{<%=generator%>}}.
#' @slot origFGNames \setsWF The original (order of) names of the \code{\link{featureGroups}} object that was used to
#'   create this object.
#'
#' @section Sets workflows: \setsWFClass{<%=class%>}{<%=parent%>}
#'
#'   \setsWFNewMethodsSO{<%=classUnset%>}{Only the annotation results that are present in the specified set are kept
#'   (based on the set consensus, see below for implications).}
#'
#'   \setsWFChangedMethods{
#'
#'   \item \code{filter} and the subset operator (\code{[}) Can be used to select data that is only present for selected
#'   sets. Depending on the \code{updateConsenus}, both either operate on set consensus or original data (see below for
#'   implications).
#'
#'   \item \code{annotatedPeakList} Returns a combined annotation table with all sets.
#'
#'   \item \code{plotSpectrum} Is able to highlight set specific mass peaks (\code{perSet} and \code{mirror} arguments).
#'
#'   \item \code{consensus} Creates the algorithm consensus based on the original annotation data (see below for
#'   implications). Furthermore, like the sets workflow method for \code{\link{<%=generator%>}}, this method supports
#'   the \code{setThreshold} and \code{setThresholdAnn} arguments.
#'
#'   <%= if (exists("extraMethods")) extraMethods %>
#'
#'   }
#'
#'   Two types of annotation data are stored in a \code{<%=class%>} object: \enumerate{
#'
#'   \item Annotations that are produced from a consensus between set results (see \code{<%=generator%>}).
#'
#'   \item The 'original' annotation data per set, prior to when the set consensus was made. This includes candidates
#'   that were filtered out because of the thresholds set by \code{setThreshold} and \code{setThresholdAnn}. However,
#'   when \code{filter} or subsetting (\code{[}) operations are performed, the original data is also updated.
#'
#'   }
#'
#'   In most cases the first data is used. However, in a few cases the original annotation data is used (as indicated
#'   above), for instance, to re-create the set consensus. It is important to realize that the original annotation data
#'   may have \emph{additional} candidates, and a newly created set consensus may therefore have 'new' candidates. For
#'   instance, when the object consists of the sets \code{"positive"} and \code{"negative"} and \code{setThreshold=1}
#'   was used to create it, then \code{<%=exObj%>[, sets = "positive", updateConsensus = TRUE]} may now have additional
#'   candidates, \emph{i.e.} those that were not present in the \code{"negative"} set and were previously removed due to
#'   the consensus threshold filter.
