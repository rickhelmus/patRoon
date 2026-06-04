#' @include workflow.R
NULL

# Wrapper method functions.
# NOTE: these are defined here to avoid circular dependencies between the workflow and the other classes.

doWfStep <- function(obj, func, slotNameIn, slotNameOut, param, paramClass, ...)
{
    if (is.null(param))
        param <- new(paramClass, template = templateDir(obj))
    slot(obj, slotNameOut) <- do.call(func, list(slot(obj, slotNameIn), param = param, ...))
    return(obj)
}

doWfFeat <- function(..., algo)
{
    doWfStep(func = paste0("findFeaturesP", algo), slotNameIn = "analysisInfo", slotNameOut = "features",
             paramClass = paste0("Features", algo, "Param"), ...)
}

#' @rdname findFeaturesOpenMS
setMethod("findFeaturesP", c("workflow", "FeaturesOpenMSParam"),
          \(obj, param = NULL, ...) doWfFeat(obj, algo = "OpenMS", param = param, ...))

#' @rdname findFeaturesOpenMS
setMethod("findFeaturesPOpenMS", "workflow",
          \(obj, param = NULL, ...) doWfFeat(obj, algo = "OpenMS", param = param, ...))

#' @rdname findFeaturesXCMS3
setMethod("findFeaturesP", c("workflow", "FeaturesXCMS3Param"),
          \(obj, param = NULL, ...) doWfFeat(obj, algo = "XCMS3", param = param, ...))

#' @rdname findFeaturesXCMS3
setMethod("findFeaturesPXCMS3", "workflow",
          \(obj, param = NULL, ...) doWfFeat(obj, algo = "XCMS3", param = param, ...))

#' @rdname findFeaturesEnviPick
setMethod("findFeaturesP", c("workflow", "FeaturesEnviPickParam"),
          \(obj, param = NULL, ...) doWfFeat(obj, algo = "EnviPick", param = param, ...))

#' @rdname findFeaturesEnviPick
setMethod("findFeaturesPEnviPick", "workflow",
          \(obj, param = NULL, ...) doWfFeat(obj, algo = "EnviPick", param = param, ...))

#' @rdname findFeaturesKPIC2
setMethod("findFeaturesP", c("workflow", "FeaturesKPIC2Param"),
          \(obj, param = NULL, ...) doWfFeat(obj, algo = "KPIC2", param = param, ...))

#' @rdname findFeaturesKPIC2
setMethod("findFeaturesPKPIC2", "workflow",
          \(obj, param = NULL, ...) doWfFeat(obj, algo = "KPIC2", param = param, ...))

#' @rdname findFeaturesPiek
setMethod("findFeaturesP", c("workflow", "FeaturesPiekParam"),
          \(obj, param = NULL, ...) doWfFeat(obj, algo = "Piek", param = param, ...))

#' @rdname findFeaturesPiek
setMethod("findFeaturesPPiek", "workflow",
          \(obj, param = NULL, ...) doWfFeat(obj, algo = "Piek", param = param, ...))
