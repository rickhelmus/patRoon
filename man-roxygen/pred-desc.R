#' @param obj The workflow object for which predictions should be performed, \emph{e.g.} feature groups with screening
#'   results (\code{\link{featureGroupsScreening}}) or compound annotations (\code{\link{compounds}}).
#' @param featureAnn A \code{\link{featureAnnotations}} object (\emph{e.g.} \code{\link{formulasSIRIUS}} or
#'   \code{\link{compounds}}) which contains <%=whatPred%>. Optional if \code{<%=calcFunc%>} is called on suspect
#'   screening results (\emph{i.e.} \code{\link{featureGroupsScreening}} method).
#' @param concUnit The concentration unit for calculated <%=whatCalc%>. Can be molar based (\code{"nM"}, \code{"uM"},
#'   \code{"mM"}, \code{"M"}) or mass based (\code{"ngL"}, \code{"ugL"}, \code{"mgL"}, \code{"gL"}). Furthermore, can be
#'   prefixed with \code{"log "} for logarithmic concentrations (\emph{e.g.} \code{"log mM"}).
#' @param type Which types of predictions should be performed: should be \code{"FP"} (\command{SIRIUS-CSI:FingerID}
#'   fingerprints), \code{"SMILES"} or \code{"both"}. Only relevant for \code{\link{compoundsSIRIUS}} method.
#' @param \dots \setsWF Further arguments passed to the non-sets workflow method.
#'
#' @section Predicting <%=whatPred%>: The <%=whatPred%> are predicted with the \code{<%=predFunc%>} generic functions,
#'   which accepts the following input:
#'
#' \itemize{
#'
#' \item \link[=suspect-screening]{Suspect screening results}. The \acronym{SMILES} data is used to predict
#' <%=whatPred%> for suspect hits.
#'
#' \item Formula annotation data obtained with \code{"sirius"} algorithm (\code{\link{generateFormulasSIRIUS}}). The
#' predictions are performed for each formula candidate using \command{SIRIUS+CSI:FingerID} fingerprints. For this
#' reason, the \code{getFingerprint} argument must be set to \code{TRUE} when generating the formula data.
#'
#' \item Compound annotation data obtained with the \code{"sirius"} algorithm (\code{\link{generateCompoundsSIRIUS}}).
#' The predictions are performed for each annotation candidate using its \acronym{SMILES} and/or
#' \command{SIRIUS+CSI:FingerID} fingerprints. The predictions are performed on a per formula basis, hence,
#' <%=whatPred%> for isomers will be equal.
#'
#' \item Compound annotation data obtained with algorithms other than \code{"sirius"}. The <%=whatPred%> are predicted
#' from \acronym{SMILES} data.
#'
#' }
#'
#'   When \acronym{SMILES} data is used then predictions of <%=whatPred%> are generally more accurate. However,
#'   calculations with \command{SIRIUS+CSI:FingerID} fingerprints are faster and only require the formula and MS/MS
#'   spectrum, \emph{i.e.} not the full structure. Hence, calculations with \acronym{SMILES} are mostly useful in
#'   suspect screening workflows, or with high confidence compound annotation data, whereas MS/MS fingerprints are
#'   suitable with unknowns.
#'
#'   For annotation data the calculations are performed for \emph{all} candidates. This can especially lead to long
#'   running calculations when \acronym{SMILES} data is used. Hence, it is \strong{strongly} recommended to first
#'   prioritize the annotation results, \emph{e.g.} with the \code{topMost} argument to the
#'   \link[=filter,featureAnnotations-method]{filter method}.
#'
#'   When <%=whatPred%> are predicted from \command{SIRIUS+CSI:FingerID} fingerprints then only formula and MS/MS
#'   spectra are used, even if compound annotations are used for input. The major difference is that with formula
#'   annotation input \emph{all} formula candidates for which a fingerprint could be generated are considered, whereas
#'   with compound annotations only candidate formulae are considered for which also a structure could be assigned.
#'   Hence, the formula annotation input could be more comprehensive, whereas predictions from structure annotations
#'   could lead to more representative results as only formulae are considered for which at least one structure could be
#'   assigned.
#'
#' @section Assigning <%=whatCalc%>: The \code{<%=calcFunc%>} generic function is used to assign <%=whatCalc%> for each
#'   feature using the <%=whatPred%> discussed in the previous section. The function takes <%=whatPred%> from suspect
#'   screening results and/or feature annotation data. If multiple <%=whatPred%> were predicted for the same feature
#'   group, for instance when multiple annotation candidates or suspect hits for this feature group are present, then a
#'   <%=whatCalc%> is assigned for all <%=whatPred%>. These values can later be easily aggregated with \emph{e.g.} the
#'   \link[=as.data.table,featureGroups-method]{as.data.table} function.
#'
#' @note The \CRANpkg{rcdk} package and \href{https://github.com/openbabel/openbabel}{OpenBabel} tool are used
#'   internally to calculate molecular weights. Please make sure that \command{OpenBabel} is installed.
#'
#' @references \insertRef{OBoyle2011}{patRoon} \cr \cr \addCitations{rcdk}{1}
