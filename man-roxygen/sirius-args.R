#' @param profile Name of the configuration profile, for example:
#'   \option{"qtof"}, \option{"orbitrap"}, \option{"fticr"}. Sets the
#'   \option{--profile} commandline option.
#' @param <%=if (ident) "formulaDatabase" else "database" %> If not \code{NULL}, use a database for retrieval of formula
#'   candidates. Possible values are: \option{"pubchem"}, \option{"bio"},
#'   \option{"kegg"}, \option{"hmdb"}. Sets the \option{--database} commandline
#'   option.
#' @param noise Median intensity of the noise (\code{NULL} ignores this
#'   parameter). Sets the \option{--noise} commandline option.
