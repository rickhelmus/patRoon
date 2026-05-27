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
#'   See the \href{https://v6.docs.sirius-ms.io/account-and-license/}{SIRIUS website} and \pkg{patRoon} handbook for
#'   more information.
#' @param projectPath These are mainly for internal purposes. \code{projectPath} sets the output directory for
#'   the \command{SIRIUS} output (a temporary directory if \code{NULL}).
#'   
#'   \setsWF \code{projectPath} should be a \code{character} specifying the paths for each set.
#'
#' @references \insertRef{Hoffmann2021}{patRoon} \cr\cr \insertRef{Dhrkop2020}{patRoon} \cr\cr
#'   \insertRef{Dhrkop2019}{patRoon} \cr\cr \insertRef{Duhrkop2015}{patRoon} \cr\cr
#'   \insertRef{Duhrkop2015-2}{patRoon} \cr\cr \insertRef{Bcker2008}{patRoon}
