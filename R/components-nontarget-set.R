#' @include main.R
#' @include components-nontarget.R
#' @include components-set.R
NULL

#' @export
componentsNTSet <- setClass("componentsNTSet", contains = "componentsSet")

setMethod("collapseComponents", "componentsNTSet", function(obj)
{
    obj@setObjects <- lapply(obj@setObjects, collapseComponents)
    obj <- syncComponentsSetObjects(obj)
    return(obj)
})

#' @export
setMethod("plotGraph", "componentsNTSet", function(obj, onlyLinked, set) plotGraph(unset(obj, set), onlyLinked = onlyLinked))

#' @export
componentsNTUnset <- setClass("componentsNTUnset", contains = "componentsNT")
setMethod("unset", "componentsNTSet", function(obj, set)
{
    cu <- callNextMethod(obj, set)
    so <- setObjects(obj)[[set]]
    return(componentsNTUnset(components = componentTable(cu), componentInfo = componentInfo(cu),
                             homol = so@homol, algorithm = paste0(algorithm(so), "_unset")))
})
