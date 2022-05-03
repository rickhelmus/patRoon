#' @describeIn <%=class%> Exports this object as a \file{.csv} file that can be used as a \command{MetFrag} local
#'   database. Any duplicate TPs (formed by different pathways or parents) will be merged based on their
#'   \acronym{InChIKey}.
#' @param out The file name of the the output \command{MetFrag} database.
#' @param includeParents Set to \code{TRUE} to include the parents in the database.
