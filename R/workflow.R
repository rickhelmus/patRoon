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
    # UNDONE: more checks for ...
    #   - ensure that all bojects are given (eg if only compounds are given, then also MSPL, fGroups etc must be there)
    #   - ensure that objects are in sync (same features, anaInfo)
    #   - if only fGroups are given, derive features/anaInfo from that. If only features, derive anaInfo from that
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo)
    callNextMethod(.Object, analysisInfo = analysisInfo, ...)
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
setMethod("report", "workflow", function(obj, MSPeakLists = obj@MSPeakLists, formulas = obj@formulas,
                                         compounds = obj@compounds, components = obj@components, TPs = obj@TPs, ...)
{
    verifyWFHasFGroups(obj)
    report(obj = obj@fGroups, MSPeakLists = MSPeakLists, formulas = formulas, compounds = compounds,
           components = components, TPs = TPs, ...)
})
