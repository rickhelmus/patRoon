#' @include feature_groups.R
NULL

#' @rdname importFeatureGroupsTable
featureGroupsTable <- setClass("featureGroupsTable", contains = "featureGroups")

setMethod("initialize", "featureGroupsTable",
          function(.Object, ...) callNextMethod(.Object, algorithm = "table", ...))

importFeatureGroupsTable <- function(analysisInfo, input)
{
    # UNDONE: doc that additional columns go in the features table
    # UNDONE: verify that groups are unique per ana
    
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo)
    checkmate::assert(
        checkmate::checkDataFrame(input),
        checkmate::checkFileExists(input, access = "r"),
        .var.name = "input"
    )
    if (is.character(input))
        input <- fread(input)
    
    gInfoCols <- c("ret", "mz", "mobility", "CCS")
    gInfoColsPrefix <- paste0("group_", gInfoCols)
    
    ac <- checkmate::makeAssertCollection()
    assertListVal(input, "group", checkmate::assertCharacter, any.missing = FALSE, min.chars = 1, add = ac)
    for (col in gInfoColsPrefix)
    {
        assertListVal(input, col, checkmate::assertNumeric, any.missing = FALSE, finite = TRUE, add = ac,
                      mustExist = FALSE)
    }
    checkmate::reportAssertions(ac)
    
    inputFeat <- copy(input)
    inputFeat <- removeDTColumnsIfPresent(inputFeat, c("group", gInfoColsPrefix, "group_mobility_collapsed",
                                                       "group_CCS_collapsed"))
    importedFeat <- importFeaturesTable(analysisInfo, inputFeat)
    
    for (col in gInfoCols)
    {
        gcol <- paste0("group_", col)
        if (is.null(input[[gcol]]) && !is.null(input[[col]]))
            input[, (gcol) := mean(get(col)), by = "group"]
    }
    
    gInfo <- unique(input, by = "group")
    gInfo <- subsetDTColumnsIfPresent(gInfo, c("group", gInfoCols, "ims_parent_group"))
    setnames(gInfo, gInfoColsPrefix, gInfoCols, skip_absent = TRUE)
    
    if (hasMobilities(inputFeat) && is.null(gInfo[["ims_parent_group"]]))
        gInfo[, ims_parent_group := NA_character_]
    
    gTable <- data.table(matrix(0, nrow = nrow(analysisInfo), ncol = nrow(gInfo)))
    setnames(gTable, gInfo$group)
    ftindex <- data.table(matrix(0L, nrow = nrow(analysisInfo), ncol = nrow(gInfo)))
    setnames(ftindex, gInfo$group)
    for (grp in gInfo$group)
    {
        ftg <- input[group == grp, c("ID", "analysis", "intensity")]
        anai <- match(ftg$analysis, analysisInfo$analysis) # align analysis order
        
        set(gTable, anai, grp, ftg$intensity)
        
        finds <- sapply(seq_len(nrow(ftg)), \(i) importedFeat[[ftg$analysis[i]]][ID == ftg$ID[i], which = TRUE])
        set(ftindex, anai, grp, finds)
    }
    
    # UNDONE: sets: fill in annotations slot (calc ion_mz column or neutralMass?), handle that there is no actual grouping algo...
    # --> make separate function for sets, with groupAlgo/groupArgs arguments
    # --> make args empty by default, and check in adducts<-()/selectIons() if set?
    ret <- featureGroupsTable(groups = gTable, groupInfo = gInfo, ftindex = ftindex, features = importedFeat)
    
    return(ret)
}
