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
#' @param config An \R6 \code{RSirius::JobSubmission} configuration object, typically obtained with
#'   \code{getSIRIUSConfig()}. If \code{NULL}, the default \command{SIRIUS} configuration is used.
#' @param runMode,projectPath Whether to execute a \command{SIRIUS} processing job (\code{runMode="execute"}) or load
#'   results from an existing \command{SIRIUS} project (\code{runMode"read"}). If \code{runMode="execute"} then
#'   \code{projectPath} can be \code{NULL} and a temporary project will be used, otherwise \code{projectPath} must point
#'   to an existing project.
#'
#'   \strong{NOTE:} if \code{runMode="execute"} then any existing project at \code{projectPath} will be removed.
#'
#'   \strong{NOTE:} This is primarily intended for internal purposes, but may be of interest to e.g. re-import SIRIUS
#'   results.
#'
#'   \setsWF \code{projectPath} should be a \code{character} specifying the paths for each set.
#' @param SIRIUSAPI An \code{rsirius_api} object for connecting to the \command{SIRIUS} API. If \code{NULL}, a new
#'   connection will be started automatically.
#' @param SIRIUSPath The full path to the \command{SIRIUS} command-line interface executable. Only used when
#'   \code{SIRIUSAPI} is not provided.
#'
#' @references \insertRef{Hoffmann2021}{patRoon} \cr\cr \insertRef{Dhrkop2020}{patRoon} \cr\cr
#'   \insertRef{Dhrkop2019}{patRoon} \cr\cr \insertRef{Duhrkop2015}{patRoon} \cr\cr
#'   \insertRef{Duhrkop2015-2}{patRoon} \cr\cr \insertRef{Bcker2008}{patRoon}
