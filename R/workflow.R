# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

verifyWFHasFGroups <- function(obj)
{
    if (is.null(obj@fGroups))
        stop("No feature groups available in this workflow object", call. = FALSE)
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
    
    hasFGroups <- !is.null(object@fGroups)
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

setMethod("templateDir", "workflow", function(obj) obj@templateDir)

#' @describeIn workflow Obtain feature group names. Requires feature groups to be available.
#' @export
setMethod("groupNames", "workflow", function(obj)
{
    verifyWFHasFGroups(obj)
    groupNames(obj@fGroups)
})

#' @describeIn workflow Returns names of the workflow data slots that contain non-NULL objects.
#' @export
setMethod("names", "workflow", function(x)
{
    c("analysisInfo", "features", "fGroups", "MSPeakLists", "formulas", "compounds", "components", "TPs")
})

#' @export
setMethod("replicates", "workflow", function(obj) unique(analysisInfo(obj)$replicate))

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
