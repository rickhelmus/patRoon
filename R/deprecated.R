# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

# nocov start

#' Deprecated and renamed functions.
#'
#' Please do not use these functions anymore since they may be removed in the
#' future.
#'
#' @param \dots Passed to successor function.
#'
#' @name patRoon-deprecated
#' @keywords internal
NULL

#' @details \code{reportHTML} performs HTML reporting, please use
#'   \code{\link{report}} instead.
#' @rdname patRoon-deprecated
#' @export
#' @keywords internal
reportHTML <- function(...)
{
    .Deprecated("report")
    report(...)
}

#' @details Please use \code{\link{estimateIDConfidence}} instead.
#' @rdname patRoon-deprecated
#' @export
#' @keywords internal
annotateSuspects <- function(...)
{
    .Deprecated("estimateIDConfidence")
    estimateIDConfidence(...)
}

# nocov end
