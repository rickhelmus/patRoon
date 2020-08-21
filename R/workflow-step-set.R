#' @include main.R
NULL

workflowStepSet <- setClass("workflowStepSet", slots = c(adducts = "list", setObjects = "list"),
                            contains = "VIRTUAL")

setMethod("setObjects", "workflowStepSet", function(obj) obj@setObjects)
setMethod("sets", "workflowStepSet", function(obj) names(setObjects(obj)))
setMethod("adducts", "workflowStepSet", function(obj) obj@adducts)

setMethod("show", "workflowStepSet", function(object)
{
    printf("Sets: %s\n", paste0(sets(object), collapse = ", "))
    printf("Adducts: %s\n", paste0(sapply(object@adducts, as.character), collapse = ", "))
})
