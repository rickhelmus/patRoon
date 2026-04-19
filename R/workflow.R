setOptClass <- \(baseName) setClassUnion(paste0(baseName, "Opt"), c(baseName, "NULL"))
setOptClass("features")
setOptClass("featureGroups")
setOptClass("MSPeakLists")
setOptClass("formulas")
setOptClass("compounds")

workflow <- setClass("workflow", slots = c(analysisInfo = "data.table", features = "featuresOpt",
                                           fGroups = "featureGroupsOpt", MSPeakLists = "MSPeakListsOpt",
                                           formulas = "formulasOpt", compounds = "compoundsOpt"))

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

setMethod("findFeatures", "workflow", function(obj, ...)
{
    obj@features <- findFeatures(analysisInfo(obj), ...)
    return(obj)
})

setMethod("findFeaturesP", c("workflow", "ANY"), function(obj, param, ...)
{
    findFeaturesP(analysisInfo(obj), param, ...)
})
