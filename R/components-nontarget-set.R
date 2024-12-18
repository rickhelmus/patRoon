# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
#' @include components-nontarget.R
#' @include components-set.R
NULL

#' @param set \setsWF The name of the set.
#'
#' @section Sets workflows: \setsWFClass{componentsNTSet}{componentsNT}
#'
#'   \setsWFNewMethodsSO{componentsNTUnset}{Only the components in the specified set are kept. Furthermore, the
#'   component names are restored to non-set specific names (see \code{\link{generateComponents}} for more details).}
#'
#'   \setsWFChangedMethods{
#'
#'   \item \code{plotGraph} Currently can only create graph networks from one set (specified by the \code{set}
#'   argument).
#'
#'   }
#'
#'   Note that the \code{componentsNTSet} class does not have a \code{homol} slot. Instead, the \code{\link{setObjects}}
#'   method can be used to access this data for a specific set.
#'
#' @rdname componentsNT-class
#' @export
componentsNTSet <- setClass("componentsNTSet", contains = "componentsSet")

setMethod("collapseComponents", "componentsNTSet", function(obj)
{
    obj@setObjects <- lapply(obj@setObjects, collapseComponents)
    obj <- syncComponentsSetObjects(obj)
    return(obj)
})

#' @rdname componentsNT-class
#' @param \dots \setsWF Further arguments passed to the non-sets workflow method.
#' @export
setMethod("plotGraph", "componentsNTSet", function(obj, onlyLinked = TRUE, set, ...) plotGraph(unset(obj, set), onlyLinked = onlyLinked, ...))

#' @rdname componentsNT-class
#' @export
componentsNTUnset <- setClass("componentsNTUnset", contains = "componentsNT")


#' @rdname componentsNT-class
#' @export
setMethod("unset", "componentsNTSet", function(obj, set)
{
    cu <- callNextMethod(obj, set)
    so <- setObjects(obj)[[set]]
    return(componentsNTUnset(components = componentTable(cu), componentInfo = componentInfo(cu),
                             homol = so@homol, algorithm = paste0(algorithm(so), "_unset")))
})
