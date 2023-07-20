#' @section Predicting <%=what%>: The <%=what%> are predicted with the \code{<%=predFunc%>} generic functions, which
#'   accepts the following input:
#'
#' \itemize{
#'
#' \item \link[=suspect-screening]{Suspect screening results}. The \acronym{SMILES} data is used to predict <%=what%>
#' for suspect hits.
#'
#' \item Formula annotation data obtained with \code{"sirius"} algorithm (\code{\link{generateFormulasSIRIUS}}).
#' The predictions are performed for each formula candidate using \command{SIRIUS+CSI:FingerID} fingerprints. For this
#' reason, the \code{getFingerprint} argument must be set to \code{TRUE} when generating the formula data.
#'
#' \item Compound annotation data obtained with the \code{"sirius"} algorithm (\code{\link{generateCompoundsSIRIUS}}).
#' The predictions are performed for each annotation candidate using its \acronym{SMILES} and/or
#' \command{SIRIUS+CSI:FingerID} fingerprints. The predictions are performed on a per formula basis, hence, <%=what%>
#' for isomers will be equal.
#'
#' \item Compound annotation data obtained with algorithms other than \code{"sirius"}. The <%=what%> are predicted from
#' \acronym{SMILES} data.
#'
#' }
#'
#'   When \acronym{SMILES} data is used then predictions of <%=what%> are generally more accurate. However, calculations
#'   with \command{SIRIUS+CSI:FingerID} fingerprints are faster and only require the formula and MS/MS spectrum,
#'   \emph{i.e.} not the full structure. Hence, calculations with \acronym{SMILES} are mostly useful in suspect
#'   screening workflows, or with high confidence compound annotation data, whereas MS/MS fingerprints are suitable with
#'   unknowns.
#'
#'   For annotation data the calculations are performed for \emph{all} candidates. This can especially lead to long
#'   running calculations when \acronym{SMILES} data is used. Hence, it is \strong{strongly} recommended to first
#'   prioritize the annotation results, \emph{e.g.} with the \code{topMost} argument to the
#'   \link[=filter,featureAnnotations-method]{filter method}.
#'
#'   When <%=what%> are predicted from \command{SIRIUS+CSI:FingerID} fingerprints then only formula and MS/MS spectra
#'   are used, even if compound annotations are used for input. The major difference is that with formula annotation
#'   input \emph{all} formula candidates for which a fingerprint could be generated are considered, whereas with
#'   compound annotations only candidate formulae are considered for which also a structure could be assigned. Hence,
#'   the formula annotation input could be more comprehensive, whereas predictions from structure annotations could lead
#'   to more representative results as only formulae are considered for which at least one structure could be assigned.
