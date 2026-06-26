#' @include workflow.R
NULL

# Wrapper method functions.
# NOTE: these are defined here to avoid circular dependencies between the workflow and the other classes.

doWfStep <- function(obj, func, slotNameIn, slotNameOut, param, paramClass, ...)
{
    if (is.null(param))
        param <- new(paramClass, template = templateDir(obj))
    args <- c(sapply(slotNameIn, slot, object = obj, simplify = FALSE), list(param = param, ...))
    names(args)[1] <- "obj" # HACK: first should always be args, any other slot inputs should remain
    slot(obj, slotNameOut) <- do.call(func, args)
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

#' @rdname findFeaturesSAFD
setMethod("findFeaturesP", c("workflow", "FeaturesSAFDParam"),
          \(obj, param = NULL, ...) doWfFeat(obj, algo = "SAFD", param = param, ...))

#' @rdname findFeaturesSAFD
setMethod("findFeaturesPSAFD", "workflow",
          \(obj, param = NULL, ...) doWfFeat(obj, algo = "SAFD", param = param, ...))

#' @rdname findFeaturesPiek
setMethod("findFeaturesP", c("workflow", "FeaturesPiekParam"),
          \(obj, param = NULL, ...) doWfFeat(obj, algo = "Piek", param = param, ...))

#' @rdname findFeaturesPiek
setMethod("findFeaturesPPiek", "workflow",
          \(obj, param = NULL, ...) doWfFeat(obj, algo = "Piek", param = param, ...))


doWfGroupFeat <- function(..., algo)
{
    doWfStep(func = paste0("groupFeaturesP", algo), slotNameIn = "features", slotNameOut = "fGroups",
             paramClass = paste0("FeatureGroups", algo, "Param"), ...)
}

#' @rdname groupFeaturesOpenMS
setMethod("groupFeaturesP", c("workflow", "FeatureGroupsOpenMSParam"),
          \(obj, param = NULL, ...) doWfGroupFeat(obj, algo = "OpenMS", param = param, ...))

#' @rdname groupFeaturesOpenMS
setMethod("groupFeaturesPOpenMS", "workflow",
          \(obj, param = NULL, ...) doWfGroupFeat(obj, algo = "OpenMS", param = param, ...))

#' @rdname groupFeaturesXCMS3
setMethod("groupFeaturesP", c("workflow", "FeatureGroupsXCMS3Param"),
          \(obj, param = NULL, ...) doWfGroupFeat(obj, algo = "XCMS3", param = param, ...))

#' @rdname groupFeaturesXCMS3
setMethod("groupFeaturesPXCMS3", "workflow",
          \(obj, param = NULL, ...) doWfGroupFeat(obj, algo = "XCMS3", param = param, ...))

#' @rdname groupFeaturesKPIC2
setMethod("groupFeaturesP", c("workflow", "FeatureGroupsKPIC2Param"),
          \(obj, param = NULL, ...) doWfGroupFeat(obj, algo = "KPIC2", param = param, ...))

#' @rdname groupFeaturesKPIC2
setMethod("groupFeaturesPKPIC2", "workflow",
          \(obj, param = NULL, ...) doWfGroupFeat(obj, algo = "KPIC2", param = param, ...))

#' @rdname groupFeaturesGreedy
setMethod("groupFeaturesP", c("workflow", "FeatureGroupsGreedyParam"),
          \(obj, param = NULL, ...) doWfGroupFeat(obj, algo = "Greedy", param = param, ...))

#' @rdname groupFeaturesGreedy
setMethod("groupFeaturesPGreedy", "workflow",
          \(obj, param = NULL, ...) doWfGroupFeat(obj, algo = "Greedy", param = param, ...))


doWfMSPL <- function(...)
{
    doWfStep(func = "generateMSPeakListsP", slotNameIn = "fGroups", slotNameOut = "MSPeakLists",
             paramClass = "MSPeakListsParam", ...)
}

#' @rdname generateMSPeakLists
setMethod("generateMSPeakListsP", "workflow", \(obj, param = NULL, ...) doWfMSPL(obj, param = param, ...))


doWfFormulas <- function(..., algo)
{
    doWfStep(func = paste0("generateFormulasP", algo), slotNameIn = c("fGroups", "MSPeakLists"),
             slotNameOut = "formulas", paramClass = paste0("Formulas", algo, "Param"), ...)
}

#' @rdname generateFormulasGenForm
setMethod("generateFormulasP", c("workflow", "FormulasGenFormParam"),
          \(obj, param = NULL, ...) doWfFormulas(obj, algo = "GenForm", param = param, ...))

#' @rdname generateFormulasGenForm
setMethod("generateFormulasPGenForm", "workflow",
          \(obj, param = NULL, ...) doWfFormulas(obj, algo = "GenForm", param = param, ...))


doWfCompounds <- function(..., algo)
{
    doWfStep(func = paste0("generateCompoundsP", algo), slotNameIn = c("fGroups", "MSPeakLists"),
             slotNameOut = "compounds", paramClass = paste0("Compounds", algo, "Param"), ...)
}

#' @rdname generateCompoundsMetFrag
setMethod("generateCompoundsPMetFrag", "workflow",
          \(obj, param = NULL, ...) doWfCompounds(obj, algo = "MetFrag", param = param, ...))

#' @rdname generateCompoundsMetFrag
setMethod("generateCompoundsP", c("workflow", "CompoundsMetFragParam"),
          \(obj, param = NULL, ...) doWfCompounds(obj, algo = "MetFrag", param = param, ...))


doWfCompon <- function(..., algo)
{
    doWfStep(func = paste0("generateComponentsP", algo), slotNameIn = "fGroups", slotNameOut = "components",
             paramClass = paste0("Components", algo, "Param"), ...)
}

#' @rdname generateComponentsRAMClustR
setMethod("generateComponentsP", c("workflow", "ComponentsRAMClustRParam"),
          \(obj, param = NULL, ...) doWfCompon(obj, algo = "RAMClustR", param = param, ...))

#' @rdname generateComponentsRAMClustR
setMethod("generateComponentsPRAMClustR", "workflow",
          \(obj, param = NULL, ...) doWfCompon(obj, algo = "RAMClustR", param = param, ...))

#' @rdname generateComponentsCAMERA
setMethod("generateComponentsP", c("workflow", "ComponentsCAMERAParam"),
          \(obj, param = NULL, ...) doWfCompon(obj, algo = "CAMERA", param = param, ...))

#' @rdname generateComponentsCAMERA
setMethod("generateComponentsPCAMERA", "workflow",
          \(obj, param = NULL, ...) doWfCompon(obj, algo = "CAMERA", param = param, ...))

#' @rdname generateComponentsIntClust
setMethod("generateComponentsP", c("workflow", "ComponentsIntClustParam"),
          \(obj, param = NULL, ...) doWfCompon(obj, algo = "IntClust", param = param, ...))

#' @rdname generateComponentsIntClust
setMethod("generateComponentsPIntClust", "workflow",
          \(obj, param = NULL, ...) doWfCompon(obj, algo = "IntClust", param = param, ...))

#' @rdname generateComponentsOpenMS
setMethod("generateComponentsP", c("workflow", "ComponentsOpenMSParam"),
          \(obj, param = NULL, ...) doWfCompon(obj, algo = "OpenMS", param = param, ...))

#' @rdname generateComponentsOpenMS
setMethod("generateComponentsPOpenMS", "workflow",
          \(obj, param = NULL, ...) doWfCompon(obj, algo = "OpenMS", param = param, ...))

#' @rdname generateComponentsCliqueMS
setMethod("generateComponentsP", c("workflow", "ComponentsCliqueMSParam"),
          \(obj, param = NULL, ...) doWfCompon(obj, algo = "CliqueMS", param = param, ...))

#' @rdname generateComponentsCliqueMS
setMethod("generateComponentsPCliqueMS", "workflow",
          \(obj, param = NULL, ...) doWfCompon(obj, algo = "CliqueMS", param = param, ...))

#' @rdname generateComponentsSpecClust
setMethod("generateComponentsPSpecClust", "workflow",
          \(obj, param = NULL, ...) doWfCompon(obj, algo = "SpecClust", param = param, ...))

#' @rdname generateComponentsTPs
setMethod("generateComponentsPTPs", "workflow",
          \(obj, param = NULL, ...) doWfCompon(obj, algo = "TPs", param = param, ...))

#' @rdname generateComponentsNontarget
setMethod("generateComponentsP", c("workflow", "ComponentsNontargetParam"),
          \(obj, param = NULL, ...) doWfCompon(obj, algo = "Nontarget", param = param, ...))

#' @rdname generateComponentsNontarget
setMethod("generateComponentsPNontarget", "workflow",
          \(obj, param = NULL, ...) doWfCompon(obj, algo = "Nontarget", param = param, ...))
