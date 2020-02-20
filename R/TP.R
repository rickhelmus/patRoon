#' @include main.R
NULL

# UNDONE: precursor/suspect --> parent?

#' Base transformation product (TP) predictions class
#'
#' Holds information for all predicted TPs for a set of suspects.
#'
#' This class holds all generated data for predicted transformation products for
#' a set of suspects. The class is \code{virtual} and derived objects are
#' created by \link[=TP-prediction]{TP predictors}.
#'
#' @param obj,x,object \code{TPPredictions} object to be accessed
#'
#' @seealso \code{\link{TP-prediction}}
#'
#' @slot suspects Table with all suspects with predictions. Use the
#'   \code{suspects} method for access.
#' @slot predictions List of predicted TPs for each suspect. Use the
#'   \code{predictions} method for access.
#'
#' @templateVar seli suspects
#' @templateVar selOrderi names()
#' @templateVar dollarOpName suspect
#' @template sub_op-args
#'
#' @templateVar class TPPredictions
#' @template class-hierarchy
#'
#' @export
TPPredictions <- setClass("TPPredictions",
                          slots = c(suspects = "data.table", predictions = "list"),
                          contains = c("VIRTUAL", "workflowStep"))

#' @describeIn TPPredictions Accessor method for the \code{suspects} slot of a
#'   \code{TPPredictions} class. This is a \code{data.table} with all suspects
#'   used for predictions.
#' @aliases suspects
#' @export
setMethod("suspects", "TPPredictions", function(pred) pred@suspects)

#' @describeIn TPPredictions Accessor method for the \code{predictions} slot of
#'   a \code{TPPredictions} class. Each TP result is stored as a
#'   \code{\link{data.table}}.
#' @aliases predictions
#' @export
setMethod("predictions", "TPPredictions", function(pred) pred@predictions)

#' @describeIn components Obtain total number of predictions.
#' @export
setMethod("length", "TPPredictions", function(x) sum(sapply(predictions(x), nrow)))

#' @describeIn TPPredictions Obtain the names of all suspects with predicted
#'   TPs.
#' @export
setMethod("names", "TPPredictions", function(x) suspects(x)$name)

#' @describeIn TPPredictions Show summary information for this object.
#' @export
setMethod("show", "TPPredictions", function(object)
{
    callNextMethod()
    
    printf("Suspects with predictions: %s (%d total)\n", getStrListWithMax(names(object), 6, ", "), nrow(suspects(object)))
    printf("Total TPs: %d\n", length(object))
})

#' @describeIn TPPredictions Subset on suspects.
#' @export
setMethod("[", c("TPPredictions", "ANY", "missing", "missing"), function(x, i, ...)
{
    if (!missing(i))
    {
        i <- assertSubsetArgAndToChr(i, names(x))
        x@predictions <- x@predictions[i]
        x@suspects <- x@suspects[name %in% i]
    }
    
    return(x)
})

#' @describeIn TPPredictions Extracts a table with TPs for a suspect.
#' @export
setMethod("[[", c("TPPredictions", "ANY", "missing"), function(x, i, j)
{
    assertExtractArg(i)
    return(x@predictions[[i]])
})

#' @describeIn TPPredictions Extracts a table with TPs for a suspect.
#' @export
setMethod("$", "TPPredictions", function(x, name)
{
    eval(substitute(x@predictions$NAME_ARG, list(NAME_ARG = name)))
})

#' @describeIn TPPredictions Returns all prediction data in a table.
#' @export
setMethod("as.data.table", "TPPredictions", function(x) rbindlist(predictions(x), idcol = "suspect"))

#' @export
setMethod("convertToSuspects", "TPPredictions", function(pred, includePrec, tidy)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFlag(includePrec, add = ac)
    checkmate::assertFlag(tidy, add = ac)
    checkmate::reportAssertions(ac)
    
    predAll <- rbindlist(predictions(pred))
    predAll <- predAll[, c("name", "mass")]
    
    if (includePrec)
        predAll <- rbind(suspects(pred), predAll, fill = TRUE)
    
    if (tidy)
        predAll <- predAll[, c("name", "mz"), with = FALSE]
    
    return(predAll)
})
