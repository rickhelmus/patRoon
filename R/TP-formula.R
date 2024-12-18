# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
#' @include TP.R
NULL

#' Base transformation products (TP) class with formula information
#'
#' Holds information for all TPs for a set of parents, including chemical formulae.
#'
#' This (virtual) class is derived from the \code{\link{transformationProducts}} base class, please see its
#' documentation for more details. Objects from this class are returned by \link[=generateTPs]{TP generators}. More
#' specifically, algorithms that works with chemical formulae (\emph{e.g.} \code{library_formula}), uses this class to
#' store their results. The methods defined for this class extend the functionality for the base
#' \code{\link{transformationProducts}} class.
#'
#' @param obj \code{transformationProductsFormula} derived object to be accessed
#'
#' @seealso The base class \code{\link{transformationProducts}} for more relevant methods and \code{\link{generateTPs}}
#'
#' @templateVar class transformationProductsFormula
#' @template class-hierarchy
#'
#' @export
transformationProductsFormula <- setClass("transformationProductsFormula",
                                          contains = c("VIRTUAL", "transformationProducts"))

setMethod("linkParentsToFGroups", "transformationProductsFormula", function(TPs, fGroups)
{
    return(screenInfo(fGroups)[name %in% names(TPs), c("name", "group"), with = FALSE])
})

#' @templateVar class transformationProductsFormula
#' @template plotGraph-TPs
#'
#' @template plotGraph
#'
#' @export
setMethod("plotGraph", "transformationProductsFormula", function(obj, which, components = NULL, prune = TRUE,
                                                                 onlyCompletePaths = FALSE, width = NULL, height = NULL)
{
    checkmate::assert(
        checkmate::checkSubset(which, names(obj), empty.ok = FALSE),
        checkmate::checkIntegerish(which, lower = 1, upper = nrow(parents(obj)), any.missing = FALSE, min.len = 1),
        .var.name = "which"
    )
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(components, "componentsTPs", null.ok = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ prune + onlyCompletePaths, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    doPlotTPGraph(as.data.table(obj[which]), parents(obj),
                  cmpTab = if (!is.null(components)) as.data.table(components) else NULL, structuresMax = 0,
                  prune = prune, onlyCompletePaths = onlyCompletePaths, width = width, height = height)
})
