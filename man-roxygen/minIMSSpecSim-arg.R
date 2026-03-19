#' @param minIMSSpecSim \IMSWF If the spectrum similarity of an IMS feature group compared to its IMS precursor (see
#'   \code{\link[=assignMobilities_feat]{assignMobilities}}) is at least this value, then the IMS feature group will not
#'   be subjected to the annotation algorithm and all feature annotation properties will be copied from its precursor.
#'   This assumes that feature annotation is primarily influenced by the MS/MS spectrum, and can be used to speed up the
#'   feature annotation process. All scorings, annotation similarities etc. are copied from the IMS precursor. The
#'   fragment annotations are also copied (\code{fragInfo} result column), however, these are adjusted based on the peak
#'   list data of the IMS feature group.
#'
#'   <%= if (is.character(append)) append else "" %>
