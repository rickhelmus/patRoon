#' @include main.R
NULL

#' (Virtual) Base class for all workflow objects.
#'
#' All workflow objects (\emph{e.g.} \code{\link{featureGroups}},
#' \code{\link{compounds}}, etc) are derived from this class. Objects from this
#' class are never created directly.
#'
#' @param obj,x An object (derived from) this class.
#'
#' @templateVar class workflowStep
#' @template class-hierarchy
#'
#' @export
workflowStep <- setClass("workflowStep", slots = c(algorithm = "character"),
                         contains = "VIRTUAL")

setValidity("workflowStep", function(object)
{
    if (length(object@algorithm) > 0 && nzchar(object@algorithm))
        return(TRUE)
    return("Algorithm not set!")
})

setMethod("initialize", "workflowStep", function(.Object, ...)
{
    .Object <- callNextMethod()
    # have to do this manually: default initialize method doesn't call validObject() if no arguments were given ...
    if (length(list(...)) == 0)
        validObject(.Object)
    .Object
})

#' @describeIn workflowStep Returns the algorithm that was used to generate an
#'   object.
#' @export
setMethod("algorithm", "workflowStep", function(obj) obj@algorithm)

#' @describeIn workflowStep Summarizes the data in this object and returns this
#'   as a \code{\link{data.table}}.
#' @param \dots Method specific arguments. Please see the documentation of the
#'   derived classes.
#' @export
setMethod("as.data.table", "workflowStep", function(x, ...) stop("VIRTUAL"))

#' @describeIn workflowStep This method simply calls \code{as.data.table} and
#'   converts the result to a classic a \code{data.frame}.
#' @param row.names,optional Ignored.
#' @export
setMethod("as.data.frame", "workflowStep", function(x, row.names = NULL,
                                                      optional = FALSE, ...) as.data.frame(as.data.table(x, ...)))

#' @describeIn workflowStep Shows summary information for this object.
#' @export
setMethod("show", "workflowStep", function(object)
{
    printf("A %s object (derived from %s)\n", class(object),
           paste0(getAllSuperClasses(getClassDef(class(object))), collapse = " -> "))
    showObjectSize(object)
    printf("Algorithm: %s\n", algorithm(object))
    invisible(NULL)
})
