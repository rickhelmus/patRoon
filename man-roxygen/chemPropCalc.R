#' @param prefCalcChemProps If \code{TRUE} then calculated chemical properties such as the formula and
#'   \acronym{InChIKey} are preferred over what is already present in the <%=whatCP%>. For efficiency reasons it is
#'   recommended to set this to \code{TRUE}. See the \verb{Validating and calculating chemical properties} section for
#'   more details.
#' @param neutralChemProps If \code{TRUE} then the neutral form of the molecule is considered to calculate
#'   \acronym{SMILES}, fomrulae etc. Enabling this may improve feature matching when considering common adducts
#'   (\emph{e.g.} \code{[M+H]+}, \code{[M-H]-}). See the \verb{Validating and calculating chemical properties} section
#'   for more details.
#'
#' @section Validating and calculating chemical properties: Chemical properties such as \acronym{SMILES},
#'   \acronym{InChIKey} and formula in the <%=whatCP%> are automatically validated and calculated if missing/invalid.
#'
#'   The internal validation/calculation process performs the following steps: \itemize{
#'
#'   \item Validation of \acronym{SMILES}, \acronym{InChI}, \acronym{InChIKey} and formula data (if present). Invalid
#'   entries will be set to \code{NA}.
#'
#'   \item If \code{neutralChemProps=TRUE} then the \acronym{SMILES}\acronym{InChI} data is neutralized by
#'   (de-)protonation (using the \command{--neutralized} option of \command{OpenBabel}). If
#'   \code{prefCalcChemProps=FALSE} then only charged molecules are considered.
#'
#'   \item The \acronym{SMILES} and \acronym{InChI} data are used to calculate missing or invalid \acronym{SMILES},
#'   \acronym{InChI}, \acronym{InChIKey} and formula data. If \code{prefCalcChemProps=TRUE} then existing
#'   \acronym{InChIKey} and formula data is overwritten by calculated values whenever possible.
#'
#'   \item The chemical formulae which were \emph{not} calculated are verified and normalized. This process may be time
#'   consuming, and is potentially largely avoided by setting \code{prefCalcChemProps=TRUE}.
#'
#'   \item Neutral masses are calculated for missing values (\code{prefCalcChemProps=FALSE}) or whenever possible
#'   (\code{prefCalcChemProps=TRUE}).
#'
#'   }
#'
#'   Note that calculation of formulae for molecules that are isotopically labelled is currently only supported for
#'   deuterium (2H) elements.
#'
#'   This functionality relies heavily on \href{http://openbabel.org/wiki/Main_Page}{OpenBabel}, please make sure it is
#'   installed.
#'
#' @references \insertRef{OBoyle2011}{patRoon}
