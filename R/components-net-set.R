# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
#' @include components-net.R
#' @include components-set.R
NULL

#' @param set \setsWF The name of the set.
#'
#' @section Sets workflows: \setsWFClass{componentsNetSet}{componentsNet}
#'
#'   \setsWFNewMethodsSO{componentsNetUnset}{Only the components in the specified set are kept. Furthermore, the
#'   component names are restored to non-set specific names (see \code{\link{generateComponents}} for more details).}
#'
#'   \setsWFChangedMethods{
#'
#'   \item \code{plotGraph} Currently can only create graph networks from one set (specified by the \code{set}
#'   argument).
#'
#'   }
#'
#'   Note that the \code{componentsNetSet} class does not have \code{featureComponents} and \code{featureGraphs} slots.
#'   Instead, the \code{\link{setObjects}} method can be used to access this data for a specific set.
#'
#' @rdname componentsNet-class
#' @export
componentsNetSet <- setClass("componentsNetSet", contains = "componentsSet")

#' @rdname componentsNet-class
#' @param \dots \setsWF Further arguments passed to the non-sets workflow method.
#' @export
setMethod("plotGraph", "componentsNetSet", function(obj, analysis, set, ...) plotGraph(unset(obj, set), analysis = analysis, ...))

#' @rdname components-class
#' @export
setMethod("expandForIMS", "componentsNetSet", function(obj, ...) cannotExpandComponIMS(obj))

#' @rdname componentsNet-class
#' @export
componentsNetUnset <- setClass("componentsNetUnset", contains = "componentsNet")

#' @rdname componentsNet-class
#' @export
setMethod("unset", "componentsNetSet", function(obj, set)
{
    cu <- callNextMethod(obj, set)
    so <- setObjects(obj)[[set]]
    return(componentsNetUnset(components = componentTable(cu), componentInfo = componentInfo(cu),
                              featureComponents = so@featureComponents, featureGraphs = so@featureGraphs,
                              annotationObjects = so@annotationObjects,
                              algorithm = paste0(algorithm(so), "_unset")))
})
