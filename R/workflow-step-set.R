# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
NULL

#' (Virtual) base class for sets related workflow objects
#'
#' This class is the base for many \link[=sets-workflow]{sets workflows} related classes. This class is virtual, and
#' therefore never created directly.
#'
#' The most important purpose of this class is to hold data that is specific for a set. These \emph{set objects} are
#' typically objects with classes from a regular non-sets workflow (\emph{e.g.} \code{\link{components}},
#' \code{\link{compounds}}), and are used by the sets workflow object to \emph{e.g.} form a consensus. Since the set
#' objects may contain additional data, such as algorithm specific slots, it may in some cases be of interest to access
#' them directly with the \code{setObjects} method (described below).
#'
#' @param obj,object An object that is derived from \code{workflowStepSet}.
#'
#' @slot setObjects A \code{list} with the \emph{set objects} (see the \verb{Details} section). The \code{list} is named
#'   with the set names.
#'
#' @templateVar class workflowStepSet
#' @template class-hierarchy
#'
#' @export
workflowStepSet <- setClass("workflowStepSet", slots = c(setObjects = "list"),
                            contains = "VIRTUAL")

#' @describeIn workflowStepSet Accessor for the \code{setObjects} slot.
#' @export
setMethod("setObjects", "workflowStepSet", function(obj) obj@setObjects)

#' @describeIn workflowStepSet Returns the names for each set in this object.
#' @export
setMethod("sets", "workflowStepSet", function(obj) names(setObjects(obj)))

#' @describeIn workflowStepSet Shows summary information for this object.
#' @export
setMethod("show", "workflowStepSet", function(object)
{
    printf("Sets: %s\n", paste0(sets(object), collapse = ", "))
})
