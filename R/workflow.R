# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

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
    # UNDONE: make anaInfo optional?
    # UNDONE: don't store anaInfo and only keep it in features?
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo)
    callNextMethod(.Object, analysisInfo = analysisInfo, ...)
})

setMethod("analysisInfo", "workflow", function(obj, df = FALSE)
{
    checkmate::assertFlag(df)
    return(if (df) as.data.frame(obj@analysisInfo) else obj@analysisInfo)
})

setMethod("templateDir", "workflow", function(obj) obj@templateDir)
