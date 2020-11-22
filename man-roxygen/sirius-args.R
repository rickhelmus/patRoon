#' @param profile Name of the configuration profile, for example:
#'   \option{"qtof"}, \option{"orbitrap"}, \option{"fticr"}. Sets the
#'   \option{--profile} commandline option.
#' @param <%=if (ident) "formulaDatabase" else "database" %> If not \code{NULL},
#'   use a database for retrieval of formula candidates. Possible values are:
#'   \option{"pubchem"}, \option{"bio"}, \option{"kegg"}, \option{"hmdb"}. Sets
#'   the \option{--database} commandline option.
#' @param noise Median intensity of the noise (\code{NULL} ignores this
#'   parameter). Sets the \option{--noise} commandline option.
#' @param cores The number of cores \command{SIRIUS} will use. If \code{NULL}
#'   then the default of all cores will be used.
#' @param extraOptsGeneral,extraOptsFormula a \code{character} vector with any
#'   extra commandline parameters for \command{SIRIUS}. For \command{SIRIUS}
#'   versions \code{<4.4} there is no distinction between general and formula
#'   options. Otherwise commandline options specified in \code{extraOptsGeneral}
#'   are added prior to the \code{formula} command, while options specified in
#'   \code{extraOptsFormula} are added in afterwards. See the \command{SIRIUS}
#'   manual for more details. Set to \code{NULL} to ignore.
#' @param splitBatches If \code{TRUE} then the calculations done by
#'   \command{SIRIUS} will be evenly split over multiple \command{SIRIUS} calls
#'   (which may be run in parallel depending on the \link[=patRoon-package]{set
#'   package options}). If \code{splitBatches=FALSE} then all feature
#'   calculations are performed from a single \command{SIRIUS} execution, which
#'   is often the fastest if calculations are performed on a single computer.
#'
#' @note For annotations performed with \command{SIRIUS} it is often the fastest
#'   to keep the default \code{SIRBatchSize=0}. In this case, the
#'   \code{maxProcAmount} argument will be ignored and all \command{SIRIUS}
#'   output will be printed to the terminal (unless \code{verbose=FALSE}).
#'
#' @references \insertRef{Dhrkop2019}{patRoon} \cr\cr
#'   \insertRef{Duhrkop2015}{patRoon} \cr\cr \insertRef{Duhrkop2015-2}{patRoon}
#'   \cr\cr \insertRef{Bcker2008}{patRoon}
