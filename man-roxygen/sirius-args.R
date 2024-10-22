#' @param profile Name of the configuration profile, for example: \option{"qtof"}, \option{"orbitrap"},
#'   \option{"fticr"}. Sets the \option{--profile} commandline option.
#' @param <%=if (ident) "formulaDatabase" else "database" %> If not \code{NULL}, use a database for retrieval of formula
#'   candidates. Possible values are: \option{"pubchem"}, \option{"bio"}, \option{"kegg"}, \option{"hmdb"}. Sets the
#'   \option{--database} commandline option.
#' @param noise Median intensity of the noise (\code{NULL} ignores this parameter). Sets the \option{--noise}
#'   commandline option.
#' @param cores The number of cores \command{SIRIUS} will use. If \code{NULL} then the default of all cores will be
#'   used.
#' @param login,alwaysLogin Specifies if and how account logging of SIRIUS should be handled:
#' 
#'   \code{login=FALSE}: no automatic login is performed and the active login status is not checked.
#'
#'   \code{login="check"}: aborts if no active login is present.
#'
#'   \code{login="interactive"}: interactively ask for login (using \CRANpkg{getPass}).
#'
#'   \code{login=c(username="...", password="...")}: perform the login with the given details. For security reasons,
#'   please do not enter the details directly, but use e.g. environment variables or store/retrieve them with the
#'   \CRANpkg{keyring} package.
#'   
#'   if \code{alwaysLogin=TRUE} then a login is always performed, otherwise only if SIRIUS reports no active login.
#' 
#'   See the \href{https://boecker-lab.github.io/docs.sirius.github.io/account-and-license/}{SIRIUS website} and
#'   \pkg{patRoon} handbook for more information.
#' @param extraOptsGeneral,extraOptsFormula a \code{character} vector with any extra commandline parameters for
#'   \command{SIRIUS}. For \command{SIRIUS} versions \code{<4.4} there is no distinction between general and formula
#'   options. Otherwise commandline options specified in \code{extraOptsGeneral} are added prior to the \code{formula}
#'   command, while options specified in \code{extraOptsFormula} are added in afterwards. See the \command{SIRIUS}
#'   manual for more details. Set to \code{NULL} to ignore.
#' @param splitBatches If \code{TRUE} then the calculations done by \command{SIRIUS} will be evenly split over multiple
#'   \command{SIRIUS} calls (which may be run in parallel depending on the \link[=patRoon-package]{set package
#'   options}). If \code{splitBatches=FALSE} then all feature calculations are performed from a single \command{SIRIUS}
#'   execution, which is often the fastest if calculations are performed on a single computer.
#' @param projectPath,dryRun These are mainly for internal purposes. \code{projectPath} sets the output directory for
#'   the \command{SIRIUS} output (a temporary directory if \code{NULL}). If \code{dryRun} is \code{TRUE} then no
#'   computations are done and only the results from \code{projectPath} are processed.
#'   
#'   \setsWF \code{projectPath} should be a \code{character} specifying the paths for each set.
#'
#' @note For annotations performed with \command{SIRIUS} it is often the fastest to keep the default
#'   \code{splitBatches=FALSE}. In this case, all \command{SIRIUS} output will be printed to the terminal (unless
#'   \code{verbose=FALSE} or \option{patRoon.MP.method="future"}). Furthermore, please note that only annotations to be
#'   performed for the same adduct are grouped in a single batch execution.
#'
#' @references \insertRef{Dhrkop2019}{patRoon} \cr\cr \insertRef{Duhrkop2015}{patRoon} \cr\cr
#'   \insertRef{Duhrkop2015-2}{patRoon} \cr\cr \insertRef{Bcker2008}{patRoon}
