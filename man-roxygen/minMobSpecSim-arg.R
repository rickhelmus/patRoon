#' @param minMobSpecSim \IMSWF If the spectrum similarity of a mobility feature group compared to its IMS parent (see
#'   \code{\link[=assignMobilities_feat]{assignMobilities}}) is at least this value, then the mobility feature group
#'   will not be subjected to the annotation algorithm and all feature annotation properties will be copied from its
#'   parent. This assumes that feature annotation is primarily influenced by the MS/MS spectrum, and can be used to
#'   speed up the feature annotation process. All scorings, annotation similarities etc. are copied from the IMS parent.
#'   The fragment annotations are also copied (\code{fragInfo} result column), however, these are adjusted based on the
#'   peak list data of the mobility feature group.
#'
#'   <%= if (is.character(append)) append else "" %>
