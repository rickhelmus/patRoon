#' @include main.R
NULL

workflowStepSet <- setClass("workflowStepSet", slots = c(setObjects = "list"),
                            contains = "VIRTUAL")

setMethod("setObjects", "workflowStepSet", function(obj) obj@setObjects)
setMethod("sets", "workflowStepSet", function(obj) names(setObjects(obj)))

setMethod("show", "workflowStepSet", function(object)
{
    printf("Sets: %s\n", paste0(sets(object), collapse = ", "))
})
