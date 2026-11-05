# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

verifyWFHasFGroups <- function(obj)
{
    if (is.null(obj@fGroups))
        stop("No feature groups available in this workflow object", call. = FALSE)
}

getWFFGroupsOrFeats <- function(obj) 
{
    if (!is.null(obj@fGroups))
        return(obj@fGroups)
    if (!is.null(obj@features))
        return(obj@features)
    stop("No feature groups or features available in this workflow object", call. = FALSE)
}

getWFFGroups <- function(obj) 
{
    verifyWFHasFGroups(obj)
    return(obj@fGroups)
}

getWFObject <- function(obj, n)
{
    n <- substitute(n)
    if (!is.character(n))
        n <- deparse(n)
    ret <- slot(obj, n)
    if (is.null(ret))
        stop(sprintf("No %s available in this workflow object", n), call. = FALSE)
    return(ret)
}

setOptClass <- \(baseName) setClassUnion(paste0(baseName, "Opt"), c(baseName, "NULL"))
setOptClass("features")
setOptClass("featureGroups")
setOptClass("MSPeakLists")
setOptClass("formulas")
setOptClass("compounds")
setOptClass("components") # UNDONE: this will become a list probably
setOptClass("transformationProducts")

workflow <- setClass("workflow", slots = c(analysisInfo = "data.table", features = "featuresOpt",
                                           fGroups = "featureGroupsOpt", MSPeakLists = "MSPeakListsOpt",
                                           formulas = "formulasOpt", compounds = "compoundsOpt",
                                           components = "componentsOpt", TPs = "transformationProductsOpt",
                                           templateDir = "character"))

setMethod("initialize", "workflow", function(.Object, analysisInfo, ...)
{
    # UNDONE: make anaInfo optional? --> at least if features are given
    # UNDONE: don't store anaInfo and only keep it in features?
    # UNDONE(?): if only fGroups are given, derive features/anaInfo from that. If only features, derive anaInfo from that
    
    # catch anaInfo assignment here for (1) validation and (2) data.table conversion
    args <- list(...)
    if (!missing(analysisInfo))
        args[["analysisInfo"]] <- assertAndPrepareAnaInfo(analysisInfo)
    do.call(callNextMethod, c(list(.Object), args))
})

setValidity("workflow", function(object)
{
    # UNDONE: ensure that objects are in sync (same features, anaInfo)?
    
    if (is.null(object@fGroups))
    {
        for (slotName in c("MSPeakLists", "formulas", "compounds", "components"))
        {
            if (!is.null(slot(object, slotName)))
                return(sprintf("%s is set, but feature groups are NULL", slotName))
        }
    }
    if (is.null(object@MSPeakLists))
    {
        for (slotName in c("formulas", "compounds"))
        {
            if (!is.null(slot(object, slotName)))
                return(sprintf("%s is set, but MSPeakLists are NULL", slotName))
        }
    }
    return(TRUE)
})

setMethod("analysisInfo", "workflow", function(obj, df = FALSE)
{
    checkmate::assertFlag(df)
    return(if (df) as.data.frame(obj@analysisInfo) else obj@analysisInfo)
})

#' @export
setMethod("analyses", "featureGroups", function(obj) analysisInfo(obj)$analysis)

setMethod("templateDir", "workflow", function(obj) obj@templateDir)

#' @export
setMethod("getFeatures", "workflow", function(obj) getFeatures(getWFFGroupsOrFeats(obj)))

#' @export
setMethod("featureTable", "workflow", function(obj, ...) featureTable(getWFFGroupsOrFeats(obj), ...))

#' @export
setMethod("getFeatureQualityNames", "workflow", function(obj, ...) getFeatureQualityNames(getWFFGroupsOrFeats(obj), ...))

#' @export
setMethod("fromIMS", "workflow", function(obj) fromIMS(getWFFGroupsOrFeats(obj)))

#' @export
setMethod("hasIMS", "workflow", function(obj) hasIMS(getWFFGroupsOrFeats(obj)))

#' @export
setMethod("groupNames", "workflow", function(obj) groupNames(getWFObj(obj, fGroups)))

#' @export
setMethod("updateGroups", "workflow", function(obj, ...) updateGroups(getWFObj(obj, fGroups), ...))

#' @export
setMethod("groupTable", "workflow", function(obj, ...) groupTable(getWFObj(obj, fGroups), ...))

#' @export
setMethod("groupFeatIndex", "workflow", function(obj) groupFeatIndex(getWFObj(obj, fGroups)))

#' @export
setMethod("groupInfo", "workflow", function(obj) groupInfo(getWFObj(obj, fGroups)))

#' @export
setMethod("groupQualities", "workflow", function(obj) groupQualities(getWFObj(obj, fGroups)))

#' @export
setMethod("groupScores", "workflow", function(obj) groupScores(getWFObj(obj, fGroups)))

#' @export
setMethod("internalStandards", "workflow", function(obj) internalStandards(getWFObj(obj, fGroups)))

#' @export
setMethod("internalStandardAssignments", "workflow", function(obj, ...) internalStandardAssignments(getWFObj(obj, fGroups), ...))

#' @export
setMethod("checkFeatures", "workflow", function(obj, ...) checkFeatures(getWFObj(obj, fGroups), ...))

#' @export
setMethod("peakLists", "workflow", function(obj, ...) peakLists(getWFObj(obj, MSPeakLists), ...))

#' @export
setMethod("averagedPeakLists", "workflow", function(obj, ...) averagedPeakLists(getWFObj(obj, MSPeakLists), ...))

#' @export
setMethod("addFormulaScoring", "workflow", function(obj, ...) addFormulaScoring(getWFObject(obj, compounds), getWFObject(obj, formulas), ...))

#' @export
setMethod("annotatedPeakList", "workflow", function(obj, index, groupName, onlyAnnotated = FALSE)
{
    if (!is.null(obj@compounds))
        return(annotatedPeakList(obj@compounds, index = index, groupName = groupName,
                                 MSPeakLists = obj@MSPeakLists, formulas = obj@formulas,
                                 onlyAnnotated = onlyAnnotated))
    if (!is.null(obj@formulas))
        return(annotatedPeakList(obj@formulas, index = index, groupName = groupName,
                                 MSPeakLists = obj@MSPeakLists, onlyAnnotated = onlyAnnotated))
    stop("No compounds or formulas available in this workflow object", call. = FALSE)
})

#' @export
setMethod("componentTable", "workflow", function(obj) componentTable(getWFObj(obj, components)))

#' @export
setMethod("componentInfo", "workflow", function(obj) componentInfo(getWFObj(obj, components)))

#' @export
setMethod("checkComponents", "workflow", function(obj, ...)
{
    checkComponents(getWFObj(obj, components), getWFObj(obj, fGroups), ...)
})

#' @export
setMethod("show", "workflow", function(object)
{
    printf("A workflow object.\n")
    showObjectSize(object)
    
    printf("Template directory: %s\n", if (length(object@templateDir) > 0) object@templateDir else "none")

    dataSlots <- names(object)[!vapply(names(object), function(n) is.null(slot(object, n)), logical(1))]
    printf("Available workflow data: %s\n", if (length(dataSlots) > 0) getStrListWithMax(dataSlots, 8, ", ") else "none")

    slotLabels <- c(features = "Features", fGroups = "Feature groups", MSPeakLists = "MS peak lists",
                    formulas = "Formulas", compounds = "Compounds", components = "Components",
                    TPs = "Transformation products")
    for (slotName in names(slotLabels))
    {
        slotObject <- slot(object, slotName)
        if (!is.null(slotObject))
            printf("%s: %d\n", slotLabels[[slotName]], length(slotObject))
    }

    showAnaInfo(analysisInfo(object))
    invisible(NULL)
})

#' @export
setMethod("getBPCs", "workflow", function(obj, ...) getBPCs(getWFFGroupsOrFeats(obj), ...))

#' @export
setMethod("getEICs", "workflow", function(obj, ...) getEICs(getWFFGroupsOrFeats(obj), ...))

#' @export
setMethod("getTICs", "workflow", function(obj, ...) getTICs(getWFFGroupsOrFeats(obj), ...))

#' @export
setMethod("names", "workflow", function(x)
{
    c("analysisInfo", "features", "fGroups", "MSPeakLists", "formulas", "compounds", "components", "TPs")
})

#' @export
setMethod("replicates", "workflow", function(obj) unique(analysisInfo(obj)$replicate))

#' @export
setMethod("[", c("workflow", "ANY", "ANY", "missing"), function(x, i, j, ..., ni, replicates, IMS, results,
                                                                reorder = FALSE, drop = TRUE)
{
    checkmate::assertFlag(reorder)
    
    doFilterByResults <- !missing(results)
    if (doFilterByResults)
    {
        if (isTRUE(results))
            results <- c("formulas", "compounds", "components")
        else
            checkmate::assertChoice(results, c("formulas", "compounds", "components"))
    }

    if (!missing(i) && !missing(ni))
        stop("Cannot simultaneously specify i and ni arguments.", call. = FALSE)

    if (!is.null(x@fGroups))
    {
        # Resolve ni here so that the same analysis selection can be used for all workflow slots.
        if (!missing(ni))
        {
            ai <- x@analysisInfo
            env <- parent.frame()
            ai <- eval(substitute(ai[.ni], list(.ni = substitute(ni))), envir = env)
            i <- ai$analysis
        }

        fgSubSetArgs <- list(x@fGroups, ..., reorder = reorder)
        if (!missing(i))
            fgSubSetArgs$i <- i
        if (!missing(j))
            fgSubSetArgs$j <- j
        if (!missing(replicates))
            fgSubSetArgs$replicates <- replicates
        if (!missing(IMS))
            fgSubSetArgs$IMS <- IMS
        x@fGroups <- do.call(`[`, fgSubSetArgs)
        
        selectedAnalyses <- analyses(x@fGroups)
        selectedGroups <- groupNames(x@fGroups)

        if (doFilterByResults)
        {
            resultGroups <- unique(unlist(lapply(results, function(slotName)
            {
                resultObj <- slot(x, slotName)
                if (is.null(resultObj))
                    return(character())
                return(groupNames(resultObj))
            }), use.names = FALSE))
            selectedGroups <- intersect(selectedGroups, resultGroups)
            x@fGroups <- x@fGroups[, selectedGroups]
        }
    }
    else
    {
        if ((!missing(replicates) && !is.null(replicates)) ||
            (!missing(IMS) && !is.null(IMS)) || doFilterByResults)
            stop("replicates, IMS, and results filters require feature groups in this workflow", call. = FALSE)

        if (!missing(ni))
        {
            ai <- x@analysisInfo
            env <- parent.frame()
            ai <- eval(substitute(ai[.ni], list(.ni = substitute(ni))), envir = env)
            i <- ai$analysis
        }
        else if (missing(i))
            i <- x@analysisInfo$analysis

        selectedAnalyses <- assertSubsetArgAndToChr(i, x@analysisInfo$analysis)
        selectedGroups <- NULL
    }

    if (reorder)
        x@analysisInfo <- x@analysisInfo[match(selectedAnalyses, analysis)]
    else
        x@analysisInfo <- x@analysisInfo[analysis %in% selectedAnalyses]

    if (!is.null(x@features))
        x@features <- x@features[selectedAnalyses, reorder = reorder]

    if (!is.null(x@MSPeakLists))
    {
        if (is.null(selectedGroups))
            x@MSPeakLists <- x@MSPeakLists[selectedAnalyses, ]
        else
            x@MSPeakLists <- x@MSPeakLists[selectedAnalyses, selectedGroups]

        if (reorder)
            x@MSPeakLists <- reorderAnalyses(x@MSPeakLists, intersect(selectedAnalyses, analyses(x@MSPeakLists)))
    }

    if (!is.null(selectedGroups))
    {
        if (!is.null(x@formulas))
            x@formulas <- x@formulas[selectedGroups]
        if (!is.null(x@compounds))
            x@compounds <- x@compounds[selectedGroups]
        if (!is.null(x@components))
            x@components <- x@components[, selectedGroups]
    }

    return(x)
})

#' @export
setMethod("[[", c("workflow", "ANY", "missing"), function(x, i, j)
{
    checkmate::assertChoice(i, names(x))
    return(slot(x, i))
})

#' @export
setReplaceMethod("[[", c("workflow", "ANY", "missing"), function(x, i, j, value)
{
    checkmate::assertChoice(i, names(x))
    slot(x, i) <- value
    validObject(x)
    return(x)
})

#' @export
setMethod("$", "workflow", function(x, name)
{
    eval(substitute(slot(x, NAME_ARG), list(NAME_ARG = name)))
})

#' @export
setReplaceMethod("$", "workflow", function(x, name, value)
{
    checkmate::assertChoice(name, names(x))
    eval(substitute(slot(x, NAME_ARG) <- VALUE_ARG, list(NAME_ARG = name, VALUE_ARG = value)))
    validObject(x)
    return(x)
})

#' @export
setMethod("report", "workflow", function(obj, MSPeakLists = obj@MSPeakLists, formulas = obj@formulas,
                                         compounds = obj@compounds, components = obj@components, TPs = obj@TPs, ...)
{
    verifyWFHasFGroups(obj)
    report(obj = obj@fGroups, MSPeakLists = MSPeakLists, formulas = formulas, compounds = compounds,
           components = components, TPs = TPs, ...)
})
