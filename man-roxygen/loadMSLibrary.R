#' @param file A \code{character} string that specifies the file path to the <%=format%> library.
#' @param potAdducts,potAdductsLib If and how missing adducts (\code{Precursor_type} data) are guessed,
#'   \code{potAdducts} should be either: \itemize{
#'
#'   \item \code{FALSE}: do not perform adduct guessing.
#'
#'   \item \code{TRUE}: guesses adducts based on a common set of known adducts (currently based on
#'   \code{\link{GenFormAdducts}} and \code{\link{MetFragAdducts}}). If \code{potAdductsLib} is \code{TRUE} then also
#'   any adducts specified in the library are used.
#'
#'   \item A \code{list} with \code{\link{adduct}} objects or \code{character} vector that can be converted with
#'   \code{\link{as.adduct}}. Only the specified adducts will be used for guessing missing values.
#'
#'   }
#' @param absMzDev The maximum absolute \emph{m/z} deviation when guessing missing adducts.
#' @param calcSPLASH If set to \code{TRUE} then missing \acronym{SPLASH} values will be calculated (see below).
#'
#' @section Automatic curation of library data: Several strategies are applied to automatically verify and improve
#'   library data. This is important, since library records may have inconsistent or erroneous data, which makes them
#'   unsuitable in automated workflows such as compounds annotation with \code{\link{generateCompoundsLibrary}}.
#'
#'   The loaded library data is post-treated as follows: \itemize{
#'
#'   \item The \code{DB#} field is renamed to \code{DB_ID} to improve compatibility with \R column names.
#'   
#'   \item Synonyms (\code{Synon} fields) are merged together, mainly to save memory usage.
#'
#'   \item Inconsistently formatted \code{NA} data (\emph{e.g.} \code{"n/a"}, \code{"N/A"} or empty strings) are set to
#'   regular \R \code{NA} values.
#'
#'   \item The case of record field names are made consistent.
#'
#'   \item The \code{Formula} and \code{ExactMass} fields are renamed to \code{formula} and \code{neutralMass},
#'   respectively. This is for consistency with other data generated with \pkg{patRoon}.
#'
#'   \item \code{character} field data is trimmed from leading/trailing whitespace.
#'
#'   \item Mass data is verified to be properly numeric, and set to \code{NA} otherwise.
#'
#'   \item The format of formulae data is made consistent: ionic species (with or without square brackets) or converted
#'   to a regular formula format.
#'
#'   \item Chemical identifiers such as \acronym{SMILES} and formulae are verified and missing values are calculated if
#'   possible. See XXX for more details.
#'
#'   \item Shortened data in the \code{Ion_mode} field (\acronym{P}/\acronym{N}) is converted to the long format
#'   (\code{POSITIVE}/\code{NEGATIVE}).
#'
#'   \item Many different adduct flavors typically found as \code{Precursor_type} data are converted and normalized to
#'   the generic textual format used by \pkg{patRoon} (see \code{\link{as.adduct}}).
#'
#'   \item If \code{potAdducts!=FALSE} then missing or invalid adduct data in \code{Precursor_type} is guessed based on
#'   the difference between the neutral and ionic mass. If multiple adducts explain the mass difference the result is
#'   \code{NA}.
#'
#'   \item Missing ion \emph{m/z} data (\code{PrecursorMZ} field) is calculated from adduct data, if possible.
#'
#'   \item Missing \href{https://splash.fiehnlab.ucdavis.edu/}{SPLASH} data is calculated with the \pkg{splashR} package
#'   if \code{calcSPLASH=TRUE}.
#'
#'   }
#'
#' @return The loaded data is returned in an \code{\link{MSLibrary}} object.
#'
#' @seealso The \code{\link{MSLibrary}} documentation for various methods to post-process the data and
#'   \code{\link{generateCompoundsLibrary}} for annotation of features with the library data.
#'
#' @source Guessing adducts from neutral/ionic mass differences was inspired from
#'   \href{http://ipb-halle.github.io/MetFrag/}{MetFrag}.
#'
#' @references \insertRef{Wohlgemuth2016}{patRoon} \cr\cr \insertRef{Ruttkies2016}{patRoon} \cr\cr
#'   \addCitations{Rcpp}{1} \cr\cr \addCitations{Rcpp}{2} \cr\cr \addCitations{Rcpp}{3}
#'   
