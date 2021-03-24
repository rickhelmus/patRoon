#' @include main.R
#' @include workflow-step.R
NULL

#' Base transformation product (TP) predictions class
#'
#' Holds information for all predicted TPs for a set of parents.
#'
#' This class holds all generated data for predicted transformation products for
#' a set of parents. The class is \code{virtual} and derived objects are
#' created by \link[=TP-prediction]{TP predictors}.
#'
#' @param obj,x,object \code{TPPredictions} object to be accessed
#'
#' @seealso \code{\link{TP-prediction}}
#'
#' @slot parents Table with all parents with predictions. Use the
#'   \code{parents} method for access.
#' @slot predictions List of predicted TPs for each suspect. Use the
#'   \code{predictions} method for access.
#'
#' @templateVar seli parents
#' @templateVar selOrderi names()
#' @templateVar dollarOpName suspect
#' @template sub_op-args
#'
#' @templateVar class TPPredictions
#' @template class-hierarchy
#'
#' @export
TPPredictions <- setClass("TPPredictions",
                          slots = c(parents = "data.table", predictions = "list"),
                          contains = c("VIRTUAL", "workflowStep"))

#' @describeIn TPPredictions Accessor method for the \code{parents} slot of a
#'   \code{TPPredictions} class. This is a \code{data.table} with all parents
#'   used for predictions.
#' @aliases parents
#' @export
setMethod("parents", "TPPredictions", function(pred) pred@parents)

#' @describeIn TPPredictions Accessor method for the \code{predictions} slot of
#'   a \code{TPPredictions} class. Each TP result is stored as a
#'   \code{\link{data.table}}.
#' @aliases predictions
#' @export
setMethod("predictions", "TPPredictions", function(pred) pred@predictions)

#' @describeIn components Obtain total number of predictions.
#' @export
setMethod("length", "TPPredictions", function(x) if (length(predictions(x)) == 0) 0 else sum(sapply(predictions(x), nrow)))

#' @describeIn TPPredictions Obtain the names of all parents with predicted
#'   TPs.
#' @export
setMethod("names", "TPPredictions", function(x) parents(x)$name)

#' @describeIn TPPredictions Show summary information for this object.
#' @export
setMethod("show", "TPPredictions", function(object)
{
    callNextMethod()
    
    printf("Parents: %s (%d total)\n", getStrListWithMax(names(object), 6, ", "), nrow(parents(object)))
    printf("Total TPs: %d\n", length(object))
})

#' @describeIn TPPredictions Subset on parents.
#' @export
setMethod("[", c("TPPredictions", "ANY", "missing", "missing"), function(x, i, ...)
{
    if (!missing(i))
    {
        i <- assertSubsetArgAndToChr(i, names(x))
        x@predictions <- x@predictions[i]
        x@parents <- x@parents[name %in% i]
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
setMethod("convertToSuspects", "TPPredictions", function(pred, includeParents)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFlag(includeParents, add = ac)
    checkmate::reportAssertions(ac)
    
    predAll <- rbindlist(predictions(pred))
    keepCols <- c("name", "SMILES", "InChI", "InChIKey", "formula", "neutralMass")
    predAll <- predAll[, intersect(keepCols, names(predAll)), with = FALSE]
    predAll <- prepareSuspectList(predAll, NULL, FALSE, FALSE)
    
    if (includeParents)
        predAll <- rbind(parents(pred), predAll, fill = TRUE)
    
    return(predAll)
})

setMethod("needsScreening", "TPPredictions", function(pred) TRUE)

setMethod("linkTPsToFGroups", "TPPredictions", function(pred, fGroups)
{
    TPNames <- as.data.table(pred)$name
    ret <- screenInfo(fGroups)[name %in% TPNames, c("group", "name"), with = FALSE]
    setnames(ret, "name", "TP_name")
    return(ret)
})
