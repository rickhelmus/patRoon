#' @include main.R
#' @include components-nontarget.R
#' @include components-set.R
NULL

componentsNTSet <- setClass("componentsNTSet", contains = "componentsSet")

#' @export
setMethod("plotGraph", "componentsNTSet", function(obj, onlyLinked, sets) plotGraph(unset(obj, sets), onlyLinked = onlyLinked))

componentsNTUnset <- setClass("componentsNTUnset", contains = "componentsNT")
setMethod("unset", "componentsNTSet", function(obj, sets)
{
    cu <- callNextMethod(obj, sets)
    so <- setObjects(obj)[[sets]] # sets can only be length 1
    return(componentsNTUnset(components = componentTable(cu), componentInfo = componentInfo(cu),
                             homol = so@homol, algorithm = paste0(algorithm(so), "_unset")))
})
