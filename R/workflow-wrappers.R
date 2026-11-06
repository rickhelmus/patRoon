# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include workflow.R
#' @include workflow-set.R
NULL

# Wrapper method functions.
# NOTE: these are defined here to avoid circular dependencies between the workflow and the other classes.

doWfStepFunc <- function(func, inputObjs, param, inputs = list(), setFirstArgObj = TRUE)
{
    args <- c(inputObjs, inputs, list(param = param))
    if (setFirstArgObj)
        names(args)[1] <- "obj" # HACK: set first arg to "obj", any other slot inputs should remain
    do.call(func, args)
}

makeSlotNameInNamed <- function(slotNameIn)
{
    names(slotNameIn) <- if (!is.null(names(slotNameIn)))
        fifelse(nzchar(names(slotNameIn)), names(slotNameIn), slotNameIn)
    else
        slotNameIn
    return(slotNameIn)
}

resolveWfStepArgs <- function(param, paramClass, args)
{
    # Route arguments matching a parameter field into the parameter object; keep method-specific arguments as inputs.
    argNames <- names(args)
    paramNames <- names(param)
    isParamArg <- !is.null(argNames) & nzchar(argNames) & argNames %in% paramNames
    isParamArg[is.na(isParamArg)] <- FALSE

    paramOverrides <- args[isParamArg]
    functionInputs <- args[!isParamArg]
    if (length(paramOverrides) > 0)
        param <- do.call(resolveParamOverrides, c(list(param = param), paramOverrides))

    return(list(param = param, functionInputs = functionInputs))
}

doWfStep <- function(object, func, slotNameIn, slotNameOut, paramClass, param = NULL, inputs = list(),
                     setFirstArgObj = TRUE, ...)
{
    slotNameIn <- makeSlotNameInNamed(slotNameIn)

    if (is.null(param))
        param <- new(paramClass, template = templateDir(object))
    else
        checkmate::assertClass(param, paramClass)
    
    resolved <- resolveWfStepArgs(param, paramClass, list(...))
    param <- resolved$param
    inputs <- c(inputs, resolved$functionInputs)

    slot(object, slotNameOut) <- doWfStepFunc(func, sapply(slotNameIn, slot, object = object, simplify = FALSE), param,
                                              inputs = inputs, setFirstArgObj = setFirstArgObj)
    
    step <- data.table(func = func, paramClass = paramClass, param = list(param),
                       slotNameIn = list(slotNameIn), slotNameOut = slotNameOut,
                       set = NA_character_)
    
    object@steps <- rbind(object@steps, step)
    return(object)
}

doWfFeatSet <- function(object, func, slotNameIn, slotNameOut, paramClass, param = NULL, inputs = list(), ..., sets = NULL)
{
    assertSets(object, sets, multiple = TRUE, null.ok = TRUE)
    
    if (is.null(param))
        param <- new(paramClass, template = templateDir(object))
    else
        checkmate::assertClass(param, paramClass)
    
    mySets <- sets(object)
    if (length(mySets) == 0)
    {
        warning("No sets available in this workflowSet object", call. = FALSE)
        return(object)
    }
    
    if (is.null(sets))
        sets <- mySets
    
    slotNameIn <- makeSlotNameInNamed(slotNameIn)

    resolved <- resolveWfStepArgs(param, paramClass, list(...))
    param <- resolved$param
    inputs <- c(inputs, resolved$functionInputs)
    
    stepRows <- list()
    object@setObjects[sets] <- Map(sets, setObjects(object)[sets], f = function(sn, so)
    {
        if (is.null(so[[slotNameIn]]))
            warning(sprintf("No %s available in set '%s'", slotNameIn, sn), call. = FALSE)
        else
        {
            so[[slotNameOut]] <- doWfStepFunc(func, setNames(so[slotNameIn], names(slotNameIn)), param,
                                              inputs = c(inputs, resolved$inputs))
            stepRows[[length(stepRows) + 1L]] <<- data.table(func = func, paramClass = paramClass,
                                                             param = list(param),
                                                             slotNameIn = list(slotNameIn),
                                                             slotNameOut = slotNameOut, set = sn)
        }
        return(so)
    })

    if (length(stepRows) > 0)
        object@steps <- rbindlist(c(list(object@steps), stepRows))
    
    return(object)
}


#' @rdname MSConversion
setMethod("convertMSFilesP", c("workflow", "ANY"), function(obj, param, ...)
{
    convertMSFilesP(analysisInfo(obj), param, ...)
    return(obj)
})


setMethod("doWfFeat", "workflow", function(obj, ..., algo)
{
    doWfStep(obj, func = paste0("findFeaturesP", algo), slotNameIn = "analysisInfo", slotNameOut = "features",
             paramClass = paste0("Features", algo, "Param"), ...)
})

setMethod("doWfFeat", "workflowSet", function(obj, ..., algo)
{
    doWfFeatSet(obj, func = paste0("findFeaturesP", algo), slotNameIn = "analysisInfo", slotNameOut = "features",
                paramClass = paste0("Features", algo, "Param"), ...)
})

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

#' @rdname findFeaturesSIRIUS
setMethod("findFeaturesP", c("workflow", "FeaturesSIRIUSParam"),
          \(obj, param = NULL, ...) doWfFeat(obj, algo = "SIRIUS", param = param, ...))

#' @rdname findFeaturesSIRIUS
setMethod("findFeaturesPSIRIUS", "workflow",
          \(obj, param = NULL, ...) doWfFeat(obj, algo = "SIRIUS", param = param, ...))

#' @rdname findFeaturesPiek
setMethod("findFeaturesP", c("workflow", "FeaturesPiekParam"),
          \(obj, param = NULL, ...) doWfFeat(obj, algo = "Piek", param = param, ...))

#' @rdname findFeaturesPiek
setMethod("findFeaturesPPiek", "workflow",
          \(obj, param = NULL, ...) doWfFeat(obj, algo = "Piek", param = param, ...))


setMethod("doWfGroupFeat", "workflow", function(obj, ..., algo, slotIn = "features")
{
    doWfStep(obj, func = paste0("groupFeaturesP", algo), slotNameIn = slotIn, slotNameOut = "fGroups",
             paramClass = paste0("FeatureGroups", algo, "Param"), ...)
})

setMethod("doWfGroupFeat", "workflowSet", function(obj, ..., algo, slotIn = "features")
{
    # if there already is a features object, makeSet was called and we can call groupFeatures() as usual. Otherwise,
    # call it per set.
    # UNDONE: make sure that this is documented
    f <- if (is.null(obj@features)) doWfFeatSet else doWfStep
    f(obj, func = paste0("groupFeaturesP", algo), slotNameIn = slotIn, slotNameOut = "fGroups",
      paramClass = paste0("FeatureGroups", algo, "Param"), ...)
})

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

#' @rdname groupFeaturesSIRIUS
setMethod("groupFeaturesP", c("workflow", "FeatureGroupsSIRIUSParam"),
          \(obj, param = NULL, ...) doWfGroupFeat(obj, algo = "SIRIUS", param = param, ..., slotIn = "analysisInfo"))

#' @rdname groupFeaturesSIRIUS
setMethod("groupFeaturesPSIRIUS", "workflow",
          \(obj, param = NULL, ...) doWfGroupFeat(obj, algo = "SIRIUS", param = param, ..., slotIn = "analysisInfo"))

#' @rdname groupFeaturesGreedy
setMethod("groupFeaturesP", c("workflow", "FeatureGroupsGreedyParam"),
          \(obj, param = NULL, ...) doWfGroupFeat(obj, algo = "Greedy", param = param, ...))

#' @rdname groupFeaturesGreedy
setMethod("groupFeaturesPGreedy", "workflow",
          \(obj, param = NULL, ...) doWfGroupFeat(obj, algo = "Greedy", param = param, ...))


doWfMSPL <- function(obj, ...)
{
    doWfStep(obj, func = "generateMSPeakListsP", slotNameIn = "fGroups", slotNameOut = "MSPeakLists",
             paramClass = "MSPeakListsParam", ...)
}

#' @rdname generateMSPeakLists
setMethod("generateMSPeakListsP", "workflow", \(obj, param = NULL, ...) doWfMSPL(obj, param = param, ...))


doWfFormulas <- function(obj, ..., algo)
{
    doWfStep(obj, func = paste0("generateFormulasP", algo), slotNameIn = c("fGroups", "MSPeakLists"),
             slotNameOut = "formulas", paramClass = paste0("Formulas", algo, "Param"), ...)
}

#' @rdname generateFormulasGenForm
setMethod("generateFormulasP", c("workflow", "FormulasGenFormParam"),
          \(obj, param = NULL, ...) doWfFormulas(obj, algo = "GenForm", param = param, ...))

#' @rdname generateFormulasGenForm
setMethod("generateFormulasPGenForm", "workflow",
          \(obj, param = NULL, ...) doWfFormulas(obj, algo = "GenForm", param = param, ...))

#' @rdname generateFormulasSIRIUS
setMethod("generateFormulasPSIRIUS", "workflow",
          \(obj, param = NULL, ...) doWfFormulas(obj, algo = "SIRIUS", param = param, ...))

#' @rdname generateFormulasSIRIUS
setMethod("generateFormulasP", c("workflow", "FormulasSIRIUSParam"),
          \(obj, param = NULL, ...) doWfFormulas(obj, algo = "SIRIUS", param = param, ...))


doWfCompounds <- function(obj, ..., algo)
{
    doWfStep(obj, func = paste0("generateCompoundsP", algo), slotNameIn = c("fGroups", "MSPeakLists"),
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

#' @rdname generateCompoundsSIRIUS
setMethod("generateCompoundsPSIRIUS", "workflow",
          \(obj, param = NULL, ...) doWfCompounds(obj, algo = "SIRIUS", param = param, ...))

#' @rdname generateCompoundsSIRIUS
setMethod("generateCompoundsP", c("workflow", "CompoundsSIRIUSParam"),
          \(obj, param = NULL, ...) doWfCompounds(obj, algo = "SIRIUS", param = param, ...))


doWfCompon <- function(obj, ..., algo, slotNameIn = "fGroups", inputs = list())
{
    doWfStep(obj, func = paste0("generateComponentsP", algo), slotNameIn = slotNameIn, slotNameOut = "components",
             paramClass = paste0("Components", algo, "Param"), inputs = inputs, ...)
}

#' @rdname generateComponentsRAMClustR
setMethod("generateComponentsP", c("workflow", "ComponentsRAMClustRParam"),
          \(obj, param = NULL, ..., ionization = NULL) doWfCompon(obj, algo = "RAMClustR", param = param,
                                                                  inputs = list(ionization = ionization), ...))

#' @rdname generateComponentsRAMClustR
setMethod("generateComponentsPRAMClustR", "workflow",
          \(obj, param = NULL, ..., ionization = NULL) doWfCompon(obj, algo = "RAMClustR", param = param,
                                                                  inputs = list(ionization = ionization), ...))

#' @rdname generateComponentsCAMERA
setMethod("generateComponentsP", c("workflow", "ComponentsCAMERAParam"),
          \(obj, param = NULL, ..., ionization = NULL) doWfCompon(obj, algo = "CAMERA", param = param,
                                                                  inputs = list(ionization = ionization), ...))

#' @rdname generateComponentsCAMERA
setMethod("generateComponentsPCAMERA", "workflow",
          \(obj, param = NULL, ..., ionization = NULL) doWfCompon(obj, algo = "CAMERA", param = param,
                                                                  inputs = list(ionization = ionization), ...))

#' @rdname generateComponentsIntClust
setMethod("generateComponentsP", c("workflow", "ComponentsIntClustParam"),
          \(obj, param = NULL, ...) doWfCompon(obj, algo = "IntClust", param = param, ...))

#' @rdname generateComponentsIntClust
setMethod("generateComponentsPIntClust", "workflow",
          \(obj, param = NULL, ...) doWfCompon(obj, algo = "IntClust", param = param, ...))

#' @rdname generateComponentsOpenMS
setMethod("generateComponentsP", c("workflow", "ComponentsOpenMSParam"),
          \(obj, param = NULL, ..., ionization = NULL) doWfCompon(obj, algo = "OpenMS", param = param,
                                                                  inputs = list(ionization = ionization), ...))

#' @rdname generateComponentsOpenMS
setMethod("generateComponentsPOpenMS", "workflow",
          \(obj, param = NULL, ..., ionization = NULL) doWfCompon(obj, algo = "OpenMS", param = param,
                                                                  inputs = list(ionization = ionization), ...))

#' @rdname generateComponentsCliqueMS
setMethod("generateComponentsP", c("workflow", "ComponentsCliqueMSParam"),
          \(obj, param = NULL, ..., ionization = NULL) doWfCompon(obj, algo = "CliqueMS", param = param,
                                                                  inputs = list(ionization = ionization), ...))

#' @rdname generateComponentsCliqueMS
setMethod("generateComponentsPCliqueMS", "workflow",
          \(obj, param = NULL, ..., ionization = NULL) doWfCompon(obj, algo = "CliqueMS", param = param,
                                                                  inputs = list(ionization = ionization), ...))

#' @rdname generateComponentsNet
setMethod("generateComponentsP", c("workflow", "ComponentsNetParam"),
          \(obj, param = NULL, ..., ionization = NULL) doWfStep(func = "generateComponentsPNet",
                                                               slotNameIn = c("fGroups", "MSPeakLists"), slotNameOut = "components",
                                                               paramClass = "ComponentsNetParam", object = obj, param = param,
                                                               inputs = list(ionization = ionization), ...))

#' @rdname generateComponentsNet
setMethod("generateComponentsPNet", "workflow",
          \(obj, param = NULL, ..., ionization = NULL) doWfStep(func = "generateComponentsPNet",
                                                               slotNameIn = c("fGroups", "MSPeakLists"), slotNameOut = "components",
                                                               paramClass = "ComponentsNetParam", object = obj, param = param,
                                                               inputs = list(ionization = ionization), ...))

#' @rdname generateComponentsSpecClust
setMethod("generateComponentsP", c("workflow", "ComponentsSpecClustParam"),
          \(obj, param = NULL, ...) doWfCompon(obj, algo = "SpecClust", param = param, ...,
                                               slotNameIn = c("fGroups", "MSPeakLists")))

#' @rdname generateComponentsSpecClust
setMethod("generateComponentsPSpecClust", "workflow",
          \(obj, param = NULL, ...) doWfCompon(obj, algo = "SpecClust", param = param, ...,
                                               slotNameIn = c("fGroups", "MSPeakLists")))

#' @rdname generateComponentsTPs
setMethod("generateComponentsP", c("workflow", "ComponentsTPsParam"),
          \(obj, param = NULL, ...) doWfCompon(obj, algo = "TPs", param = param, ...,
                                               slotNameIn = c("fGroups", "TPs", "MSPeakLists", "formulas", "compounds")))

#' @rdname generateComponentsTPs
setMethod("generateComponentsPTPs", "workflow",
          \(obj, param = NULL, ...) doWfCompon(obj, algo = "TPs", param = param, ...,
                                               slotNameIn = c("fGroups", "TPs", "MSPeakLists", "formulas", "compounds")))

#' @rdname generateComponentsNontarget
setMethod("generateComponentsP", c("workflow", "ComponentsNontargetParam"),
          \(obj, param = NULL, ..., ionization = NULL) doWfCompon(obj, algo = "Nontarget", param = param,
                                                                  inputs = list(ionization = ionization), ...))

#' @rdname generateComponentsNontarget
setMethod("generateComponentsPNontarget", "workflow",
          \(obj, param = NULL, ..., ionization = NULL) doWfCompon(obj, algo = "Nontarget", param = param,
                                                                  inputs = list(ionization = ionization), ...))


doWfTPs <- function(obj, ..., parents, algo, slotNameIn = character(), inputs = list())
{
    if (identical(parents, "fGroups"))
    {
        if (is.null(obj@fGroups) || !isScreening(obj@fGroups))
            stop("No suspect information in this workflow object! You may need to perform a prior suspect screening or adapt the parents argument.", call. = FALSE)
        parents <- obj@fGroups
    }
    else if (identical(parents, "compounds"))
    {
        if (is.null(obj@compounds))
            stop("No compound annotations in this workflow object!", call. = FALSE)
        parents <- obj@compounds
    }
    doWfStep(func = paste0("generateTPsP", algo), slotNameIn = slotNameIn, slotNameOut = "TPs",
             paramClass = paste0("TPs", algo, "Param"), object = obj, inputs = c(list(obj = parents), inputs),
             setFirstArgObj = FALSE, ...)
}

doWfTPsLogic <- function(obj, ..., algo, inputs = list())
{
    doWfStep(obj, func = paste0("generateTPsP", algo), slotNameIn = "fGroups", slotNameOut = "TPs",
             paramClass = paste0("TPs", algo, "Param"), inputs = inputs, ...)
}


#' @rdname generateTPsBioTransformer
setMethod("generateTPsP", c("workflow", "TPsBioTransformerParam"),
          \(obj, param = NULL, parents = "fGroups", ...) doWfTPs(obj, algo = "BioTransformer", param = param, ..., parents = parents))

#' @rdname generateTPsBioTransformer
setMethod("generateTPsPBioTransformer", "workflow",
          \(obj, param = NULL, parents = "fGroups", ...) doWfTPs(obj, algo = "BioTransformer", param = param, ..., parents = parents))


#' @rdname generateTPsCTS
setMethod("generateTPsP", c("workflow", "TPsCTSParam"),
          \(obj, param = NULL, parents = "fGroups", ...) doWfTPs(obj, algo = "CTS", param = param, ..., parents = parents))

#' @rdname generateTPsCTS
setMethod("generateTPsPCTS", "workflow",
          \(obj, param = NULL, parents = "fGroups", ...) doWfTPs(obj, algo = "CTS", param = param, ..., parents = parents))

#' @rdname generateTPsLibrary
setMethod("generateTPsP", c("workflow", "TPsLibraryParam"),
          \(obj, param = NULL, parents = "fGroups", ...) doWfTPs(obj, algo = "Library", param = param, ..., parents = parents))

#' @rdname generateTPsLibrary
setMethod("generateTPsPLibrary", "workflow",
          \(obj, param = NULL, parents = "fGroups", ...) doWfTPs(obj, algo = "Library", param = param, ..., parents = parents))


#' @rdname generateTPsAnnForm
setMethod("generateTPsP", c("workflow", "TPsAnnFormParam"),
          \(obj, param = NULL, parents = "fGroups", ...) doWfTPs(obj, algo = "AnnForm", param = param, ..., parents = parents, slotNameIn = "formulas"))

#' @rdname generateTPsAnnForm
setMethod("generateTPsPAnnForm", "workflow",
          \(obj, param = NULL, parents = "fGroups", ...) doWfTPs(obj, algo = "AnnForm", param = param, ..., parents = parents, slotNameIn = "formulas"))

#' @rdname generateTPsAnnComp
setMethod("generateTPsP", c("workflow", "TPsAnnCompParam"),
          \(obj, param = NULL, parents = "fGroups", ..., TPsRef = NULL) doWfTPs(obj, algo = "AnnComp", param = param,
                                                                                ..., parents = parents,
                                                                                slotNameIn = c("compounds", fGroupsComps = "fGroups"),
                                                                                inputs = list(TPsRef = TPsRef)))

#' @rdname generateTPsAnnComp
setMethod("generateTPsPAnnComp", "workflow",
          \(obj, param = NULL, parents = "fGroups", ..., TPsRef = NULL) doWfTPs(obj, algo = "AnnComp", param = param,
                                                                                ..., parents = parents,
                                                                                slotNameIn = c("compounds", fGroupsComps = "fGroups"),
                                                                                inputs = list(TPsRef = TPsRef)))

#' @rdname generateTPsLibraryFormula
setMethod("generateTPsP", c("workflow", "TPsLibraryFormulaParam"),
          \(obj, param = NULL, parents = "fGroups", ..., TPLibrary) doWfTPs(obj, algo = "LibraryFormula", param = param,
                                                                            ..., parents = parents,
                                                                            inputs = list(TPLibrary = TPLibrary)))

#' @rdname generateTPsLibraryFormula
setMethod("generateTPsPLibraryFormula", "workflow",
          \(obj, param = NULL, parents = "fGroups", ..., TPLibrary) doWfTPs(obj, algo = "LibraryFormula", param = param,
                                                                            ..., parents = parents,
                                                                            inputs = list(TPLibrary = TPLibrary)))

#' @rdname generateTPsLogic
setMethod("generateTPsP", c("workflow", "TPsLogicParam"),
          \(obj, param = NULL, ..., adduct = NULL, transformations = NULL) doWfTPsLogic(obj, algo = "Logic",
                                                                                        param = param, ...,
                                                                                        inputs = list(adduct = adduct, transformations = transformations)))

#' @rdname generateTPsLogic
setMethod("generateTPsPLogic", "workflow",
          \(obj, param = NULL, ..., adduct = NULL, transformations = NULL) doWfTPsLogic(obj, algo = "Logic",
                                                                                        param = param, ...,
                                                                                        inputs = list(adduct = adduct, transformations = transformations)))


doWfFilter <- function(obj, slotName, paramClass, param = NULL, ...)
{
    doWfStep(obj, func = "filterP", slotNameIn = slotName, slotNameOut = slotName, paramClass = paramClass,
             param = param, ...)
}

#' @rdname feature-filtering
setMethod("filterP", c("workflow", "FilterFeatGroupsParam"),
          \(obj, param = NULL, ...) doWfFilter(obj, slotName = "fGroups", paramClass = "FilterFeatGroupsParam", param = param, ...))

#' @rdname MSPeakLists-class
setMethod("filterP", c("workflow", "FilterMSPeakListsParam"),
          \(obj, param = NULL, ...) doWfFilter(obj, slotName = "MSPeakLists", paramClass = "FilterMSPeakListsParam", param = param, ...))


#' @rdname suspect-screening
setMethod("screenSuspectsP", "workflow",
          \(obj, param = NULL, ..., suspects, adduct = NULL) doWfStep(obj, func = "screenSuspectsP",
                                                                      slotNameIn = "fGroups", slotNameOut = "fGroups",
                                                                      paramClass = "ScreenSuspectsParam", param = param,
                                                                      inputs = list(suspects = suspects, adduct = adduct), ...))

#' @rdname featureGroups-class
setMethod("selectIonsP", "workflow",
          \(obj, param = NULL, ..., prefAdduct) doWfStep(obj, func = "selectIonsP",
                                                         slotNameIn = c("fGroups", "components"),
                                                         slotNameOut = "fGroups", paramClass = "SelectIonsParam",
                                                         param = param, inputs = list(prefAdduct = prefAdduct), ...))

#' @rdname featureGroups-class
setMethod("normIntsP", "workflow",
          \(obj, param = NULL, ..., standards = NULL) doWfStep(obj, func = "normIntsP", slotNameIn = "fGroups",
                                                               slotNameOut = "fGroups", paramClass = "NormIntsParam",
                                                               param = param, inputs = list(standards = standards), ...))

#' @rdname featureGroups-class
setMethod("calculatePeakQualitiesP", "workflow",
          \(obj, param = NULL, ...) doWfStep(obj, func = "calculatePeakQualitiesP", slotNameIn = "fGroups",
                                             slotNameOut = "fGroups", paramClass = "CalculatePeakQualitiesParam",
                                             param = param, ...))

#' @rdname id-conf
setMethod("estimateIDConfidenceP", "workflow", function(obj, param = NULL, ...,
                                                        IDFile = system.file("misc", "IDLevelRules.yml", package = "patRoon"),
                                                        logPath = NULL)
{
    inputs <- list(IDFile = IDFile, logPath = logPath)

    if (!is.null(obj@formulas))
    {
        obj <- doWfStep(obj, func = "estimateIDConfidenceP", slotNameIn = "formulas", slotNameOut = "formulas",
                        paramClass = "EstimateIDConfidenceParam", param = param, inputs = inputs, ...)
    }
    if (!is.null(obj@compounds))
    {
        obj <- doWfStep(obj, func = "estimateIDConfidenceP", slotNameIn = c("compounds", "formulas"),
                        slotNameOut = "compounds", paramClass = "EstimateIDConfidenceParam", param = param,
                        inputs = inputs, ...)
    }
    
    if (!is.null(obj@fGroups) && isScreening(obj@fGroups))
    {
        obj <- doWfStep(obj, func = "estimateIDConfidenceP", slotNameIn = c("fGroups", "formulas", "compounds"),
                        slotNameOut = "fGroups", paramClass = "EstimateIDConfidenceParam",  param = param,
                        inputs = inputs, ...)
    }
    
    # UNDONE: warn if nothing happened
    
    return(obj)
})

#' @rdname pred-tox
setMethod("predictToxP", "workflow", function(obj, param = NULL, ...)
{
    # UNDONE: enable this if MS2Tox will ever support SIRIUS 6
    # if (!is.null(obj@formulas) && inherits(obj@formulas, "formulasSIRIUS"))
    # {
    #     obj <- doWfStep(obj, func = "predictToxP", slotNameIn = "formulas", slotNameOut = "formulas",
    #                     paramClass = "PredictToxParam", param = param, ...)
    # }
    if (!is.null(obj@compounds))
    {
        obj <- doWfStep(obj, func = "predictToxP", slotNameIn = "compounds", slotNameOut = "compounds",
                        paramClass = "PredictToxParam", param = param, ...)
    }
    if (!is.null(obj@fGroups) && isScreening(obj@fGroups))
    {
        obj <- doWfStep(obj, func = "predictToxP", slotNameIn = "fGroups", slotNameOut = "fGroups",
                        paramClass = "PredictToxParam", param = param, ...)
    }
    return(obj)
})

#' @rdname pred-quant
setMethod("predictRespFactorsP", "workflow", function(obj, param = NULL, ..., calibrants)
{
    # UNDONE: enable this if MS2Quant will ever support SIRIUS 6
    # if (!is.null(obj@formulas) && inherits(obj@formulas, "formulasSIRIUS"))
    # {
    #     obj <- doWfStep(obj, func = "predictRespFactorsP", slotNameIn = c("formulas", "fGroups"),
    #                     slotNameOut = "formulas", paramClass = "PredictRespFactorsParam",
    #                     param = param, ...)
    # }
    if (!is.null(obj@compounds))
    {
        obj <- doWfStep(obj, func = "predictRespFactorsP", slotNameIn = c("compounds", "fGroups"),
                        slotNameOut = "compounds", paramClass = "PredictRespFactorsParam", param = param,
                        inputs = list(calibrants = calibrants), ...)
    }
    if (!is.null(obj@fGroups) && isScreening(obj@fGroups))
    {
        obj <- doWfStep(obj, func = "predictRespFactorsP", slotNameIn = "fGroups", slotNameOut = "fGroups",
                        paramClass = "PredictRespFactorsParam", param = param,
                        inputs = list(calibrants = calibrants), ...)
    }
    return(obj)
})

setMethod("assignMobilitiesP", c("workflow", "character"), function(obj, param, ...)
{
    checkmate::assertChoice(param, c("fGroups", "compounds"))
    obj <- doWfStep(obj, func = "assignMobilitiesP",
                    slotNameIn = if (param == "fGroups") "fGroups" else c("compounds", "fGroups"), slotNameOut = param,
                    paramClass = if (param == "fGroups") "AssignMobilitiesFeatureGroupsParam" else "AssignMobilitiesCompoundsParam",
                    param = NULL, ...)
    return(obj)
})
setMethod("assignMobilitiesP", c("workflow", "AssignMobilitiesFeatureGroupsParam"), function(obj, param, ...)
{
    obj <- doWfStep(obj, func = "assignMobilitiesP", slotNameIn = "fGroups", slotNameOut = "fGroups",
                    paramClass = "AssignMobilitiesFeatureGroupsParam", param = param, ...)
    return(obj)
})
setMethod("assignMobilitiesP", c("workflow", "AssignMobilitiesCompoundsParam"), function(obj, param, ...)
{
    obj <- doWfStep(obj, func = "assignMobilitiesP", slotNameIn = c("compounds", "fGroups"), slotNameOut = "compounds",
                    paramClass = "AssignMobilitiesCompoundsParam", param = param, ...)
    return(obj)
})

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
