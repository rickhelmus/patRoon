#' @details \itemize{
#'    \item \code{path} the full path to the analysis.
#'    \item \code{analysis} the filename \strong{without} extension. Must be \strong{unique},
#'    even if the \code{path} is different.
#'    \item \code{group} name of \emph{replicate group}. A replicate group is used
#'    to group analyses together that are replicates of each other. Thus, the
#'    \code{group} column for all  analyses considered to be belonging to the same
#'    replicate group should have an equal (but unique) value. Used for \emph{e.g.}
#'    avaraging and \code{\link[=filter,featureGroups-method]{filter}}.
#'    \item \code{blank}}: all analyses within this replicate group are used by
#'    \code{\link[=filter,featureGroups-method]{filter}} for blank subtraction.
#'    Multiple entries can be entered by separation with a comma.
#'    \item \code{conc} a numeric value specifying the 'concentration' of the analysis.
#'    This can be actually any kind of quantitative value such as exposure time,
#'    dilution factor or anything else which may be used to form a linear relationship.
#' }
