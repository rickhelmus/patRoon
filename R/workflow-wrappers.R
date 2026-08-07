# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include workflow.R
NULL

# Wrapper method functions.
# NOTE: these are defined here to avoid circular dependencies between the workflow and the other classes.

doWfStep <- function(obj, func, slotNameIn, slotNameOut, param, paramClass, ...)
{
    if (is.null(param))
        param <- new(paramClass, template = templateDir(obj))
    args <- c(sapply(slotNameIn, slot, object = obj, simplify = FALSE), list(param = param, ...))
    names(args)[1] <- "obj" # HACK: first should always be obj, any other slot inputs should remain
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

#' @rdname generateCompoundsLibrary
setMethod("generateCompoundsPLibrary", "workflow",
          \(obj, param = NULL, ...) doWfCompounds(obj, algo = "Library", param = param, ...))

#' @rdname generateCompoundsLibrary
setMethod("generateCompoundsP", c("workflow", "CompoundsLibraryParam"),
          \(obj, param = NULL, ...) doWfCompounds(obj, algo = "Library", param = param, ...))


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

doWfTPs <- function(..., algo)
{
    doWfStep(func = paste0("generateTPsP", algo), slotNameIn = "fGroups", slotNameOut = "TPs",
             paramClass = paste0("TPs", algo, "Param"), ...)
}

#' @rdname generateTPsBioTransformer
setMethod("generateTPsP", c("workflow", "TPsBioTransformerParam"),
          \(obj, param = NULL, ...) doWfTPs(obj, algo = "BioTransformer", param = param, ...))

#' @rdname generateTPsBioTransformer
setMethod("generateTPsPBioTransformer", "workflow",
          \(obj, param = NULL, ...) doWfTPs(obj, algo = "BioTransformer", param = param, ...))


#' @rdname generateTPsCTS
setMethod("generateTPsP", c("workflow", "TPsCTSParam"),
          \(obj, param = NULL, ...) doWfTPs(obj, algo = "CTS", param = param, ...))

#' @rdname generateTPsCTS
setMethod("generateTPsPCTS", "workflow",
          \(obj, param = NULL, ...) doWfTPs(obj, algo = "CTS", param = param, ...))

#' @rdname generateTPsLibrary
setMethod("generateTPsP", c("workflow", "TPsLibraryParam"),
          \(obj, param = NULL, ...) doWfTPs(obj, algo = "Library", param = param, ...))

#' @rdname generateTPsLibrary
setMethod("generateTPsPLibrary", "workflow",
          \(obj, param = NULL, ...) doWfTPs(obj, algo = "Library", param = param, ...))


#' @rdname generateTPsAnnForm
setMethod("generateTPsP", c("workflow", "TPsAnnFormParam"),
          \(obj, param = NULL, ...) doWfTPs(obj, algo = "AnnForm", param = param, ...))

#' @rdname generateTPsAnnForm
setMethod("generateTPsPAnnForm", "workflow",
          \(obj, param = NULL, ...) doWfTPs(obj, algo = "AnnForm", param = param, ...))

#' @rdname generateTPsAnnComp
setMethod("generateTPsP", c("workflow", "TPsAnnCompParam"),
          \(obj, param = NULL, ...) doWfTPs(obj, algo = "AnnComp", param = param, ...))

#' @rdname generateTPsAnnComp
setMethod("generateTPsPAnnComp", "workflow",
          \(obj, param = NULL, ...) doWfTPs(obj, algo = "AnnComp", param = param, ...))

#' @rdname generateTPsLibraryFormula
setMethod("generateTPsP", c("workflow", "TPsLibraryFormulaParam"),
          \(obj, param = NULL, ...) doWfTPs(obj, algo = "LibraryFormula", param = param, ...))

#' @rdname generateTPsLibraryFormula
setMethod("generateTPsPLibraryFormula", "workflow",
          \(obj, param = NULL, ...) doWfTPs(obj, algo = "LibraryFormula", param = param, ...))

#' @rdname generateTPsLogic
setMethod("generateTPsP", c("workflow", "TPsLogicParam"),
          \(obj, param = NULL, ...) doWfTPs(obj, algo = "Logic", param = param, ...))

#' @rdname generateTPsLogic
setMethod("generateTPsPLogic", "workflow",
          \(obj, param = NULL, ...) doWfTPs(obj, algo = "Logic", param = param, ...))


doWfFilter <- function(obj, slotName, paramClass, param = NULL, ...)
{
    doWfStep(func = "filterP", slotNameIn = slotName, slotNameOut = slotName, paramClass = paramClass, obj = obj,
             param = param, ...)
}

#' @rdname feature-filtering
setMethod("filterP", c("workflow", "FilterFeatGroupsParam"),
          \(obj, param = NULL, ...) doWfFilter(obj, slotName = "fGroups", paramClass = "FilterFeatGroupsParam", param = param, ...))

#' @rdname MSPeakLists-class
setMethod("filterP", c("workflow", "FilterMSPeakListsParam"),
          \(obj, param = NULL, ...) doWfFilter(obj, slotName = "MSPeakLists", paramClass = "FilterMSPeakListsParam", param = param, ...))


setMethod("wfWrap", "workflow", function(obj, expr)
{
    slotMap <- names(obj)
    names(slotMap) <- paste0(".", slotMap)
    names(slotMap)[slotMap == "MSPeakLists"] <- ".mslists"

    expr <- substitute(expr)

    if (!is.call(expr))
        stop("Expression must be an assignment (e.g. .fGroups <- ...")

    op <- as.character(expr[[1]])
    if (!(op %in% c("<-", "=")))
        stop("Expression must be an assignment using '<-' or '='")

    lhs <- expr[[2]]
    # Require the left-hand side to be a dot-prefixed symbol naming a workflow slot
    if (!is.symbol(lhs))
        stop("Left-hand side must be a dot-prefixed symbol naming a workflow slot (e.g. .fGroups)")

    lhsName <- as.character(lhs)
    if (!(lhsName %in% names(slotMap)))
        stop("Left-hand side must be a dot-prefixed symbol naming a workflow slot (e.g. .fGroups)")

    # Replace dot-prefixed symbols in the expression with slot(obj, "<name>") calls
    replaceSlots <- function(e, isFuncPos = FALSE)
    {
        if (is.symbol(e))
        {
            symName <- as.character(e)
            if (!isFuncPos && nchar(symName) > 1 && substr(symName, 1, 1) == ".")
            {
                if (symName %in% names(slotMap))
                {
                    slotName <- slotMap[[symName]]
                    return(substitute(slot(obj, NAME), list(NAME = slotName)))
                }
            }
            return(e)
        }
        else if (is.call(e))
        {
            # Process the call: recurse into the function position (which we never replace) and all argument positions
            # (which we may replace).
            funcPos <- e[[1]]
            args <- if (length(e) > 1) as.list(e[-1]) else list()
            newArgs <- lapply(args, function(a) replaceSlots(a, isFuncPos = FALSE))
            as.call(c(list(replaceSlots(funcPos, isFuncPos = TRUE)), newArgs))
        }
        else if (is.pairlist(e))
            as.pairlist(lapply(e, function(x) replaceSlots(x, isFuncPos = FALSE)))
        else if (is.expression(e))
            as.expression(lapply(e, function(x) replaceSlots(x, isFuncPos = FALSE)))
        else
            e # Atomic vectors, NULL, etc. – leave as-is
    }

    exprReplaced <- replaceSlots(expr)

    # Evaluate the replaced assignment in a temporary environment so that slot<- updates the 'obj' binding
    evalEnv <- new.env(parent = parent.frame())
    evalEnv$obj <- obj
    eval(exprReplaced, envir = evalEnv, enclos = parent.frame())

    # Retrieve the (possibly) modified workflow object and return it
    obj <- evalEnv$obj

    return(obj)
})
