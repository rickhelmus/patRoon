#' @include feature_groups.R
NULL

#' @rdname importFeatureGroupsTable
featureGroupsTable <- setClass("featureGroupsTable", contains = "featureGroups")

setMethod("initialize", "featureGroupsTable",
          function(.Object, ...) callNextMethod(.Object, algorithm = "table", ...))

importFeatureGroupsTable <- function(analysisInfo, input, addCols = NULL, groupAlgo, groupArgs = NULL)
{
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo)
    checkmate::assert(
        checkmate::checkDataFrame(input),
        checkmate::checkFileExists(input, access = "r"),
        .var.name = "input"
    )
    input <- if (is.character(input)) fread(input) else makeDT(input)
    
    gInfoCols <- c("ret", "mz", "mobility", "CCS")
    gInfoColsPrefix <- paste0("group_", gInfoCols)
    
    hasSets <- !is.null(input[["set"]])
    sufSets <- \(x) paste0(x, "-", setsNames)
    setsCol <- \(x, s) paste0(x, "-", s) 
    setsNames <- setsAnnAddCols <- setsAnnNumCols <- character()
    if (hasSets)
    {
        setsNames <- unique(input$set)
        setsAnnNumCols <- c(sufSets("group_ion_mz"), "group_neutralMass")
        setsAnnAddCols <- sufSets("group_adduct")
    }
    
    ac <- checkmate::makeAssertCollection()
    # NOTE: groupAlgo is checked below and only if sets data is being imported
    checkmate::assertList(groupArgs, null.ok = TRUE, add = ac)
    
    assertListVal(input, "group", checkmate::assertCharacter, any.missing = FALSE, min.chars = 1, add = ac)
    assertUniqueDTBy(input, c("group", "analysis"), add = ac)
    for (col in gInfoColsPrefix)
    {
        assertListVal(input, col, checkmate::assertNumeric, any.missing = FALSE, finite = TRUE, add = ac,
                      mustExist = FALSE)
    }
    if (hasSets)
    {
        for (col in setsAnnNumCols)
        {
            assertListVal(input, col, checkmate::assertNumeric, any.missing = col != "group_neutralMass", finite = TRUE,
                          add = ac, mustExist = FALSE)
        }
        for (col in setsAnnAddCols)
        {
            assertListVal(input, col, checkmate::assertCharacter, any.missing = TRUE, min.chars = 1, add = ac,
                          mustExist = FALSE)
        }
    }
    checkmate::reportAssertions(ac)
    
    setsAnn <- NULL
    if (hasSets)
    {
        # fill in annotations table: either from group columns (as exported by as.data.table()) or from feature columns
        # adduct: only reported as group_adduct-<set>
        # neutralMass: as mz column per feature, or group_neutralMass per group
        # ion_mz: per feature or group_ion_mz-<set>
        # --> make sure that both feature and feature group variations are present, this is needed for importing the
        # features and making the sets annotations table
        
        # adduct: not reported by ADT, feats to groups: take from unique by group table, group to feats: copy from corresponding set col
        # neutralMass: feats to groups: mean average by group; group to feats: not needed
        # ion_mz: feats to groups: mean average, group to feats: copy

        ensureCols <- function(featCol, setsCols)
        {
            if (!featCol %in% names(input) && !all(setsCols %in% names(input)))
                stop("Please specify at least one of the following columns: ",
                     paste0(c(featCol, setsCols), collapse = ", "), call. = FALSE)
        }
        
        # ensure either feature or fGroup column is present
        ensureCols("adduct", sufSets("group_adduct"))
        ensureCols("mz", "group_neutralMass")
        ensureCols("ion_mz", sufSets("group_ion_mz"))
        
        # add missing feature columns from group data (mz from neutralMass not needed: mz should always be there)
        if (is.null(input[["adduct"]]))
            input[, adduct := get(paste0("group_adduct-", set)), by = "set"]
        if (is.null(input[["ion_mz"]]))
            input[, ion_mz := get(paste0("group_ion_mz-", set)), by = "set"]
        
        for (s in setsNames)
        {
            cn <- setsCol("group_adduct", s)
            if (is.null(input[[cn]]))
                input[, (cn) := adduct]
            cn <- setsCol("group_ion_mz", s)
            if (is.null(input[[cn]]))
                input[, (cn) := mean(ion_mz), by = c("group", "set")]
        }
        if (is.null(input[["neutralMass"]]))
            input[, neutralMass := mean(mz), by = "group"]
    }
    
    inputFeat <- copy(input)
    inputFeat <- removeDTColumnsIfPresent(inputFeat, c("group", gInfoColsPrefix, "group_mobility_collapsed",
                                                       "group_CCS_collapsed", setsAnnAddCols, setsAnnNumCols))
    importedFeat <- importFeaturesTable(analysisInfo, inputFeat, addCols = addCols)
    
    for (col in gInfoCols)
    {
        gcol <- paste0("group_", col)
        if (is.null(input[[gcol]]) && !is.null(input[[col]]))
            input[, (gcol) := mean(get(col)), by = "group"]
    }
    
    gInfo <- unique(input, by = "group")
    gInfo <- subsetDTColumnsIfPresent(gInfo, c("group", gInfoCols, "ims_parent_group"))
    setnames(gInfo, gInfoColsPrefix, gInfoCols, skip_absent = TRUE)
    
    if (hasMobilities(importedFeat) && is.null(gInfo[["ims_parent_group"]]))
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
    
    constArgs <- list(groups = gTable, groupInfo = gInfo, ftindex = ftindex, features = importedFeat)
    ret <- if (hasSets)
    {
        assertGroupFeatAlgo(groupAlgo)
        
        ann <- unique(input, by = c("set", "group"))[, c("set", "group", setsAnnAddCols, setsAnnNumCols), with = FALSE]
        for (s in setsNames)
        {
            cols <- c(setsCol("group_adduct", s), setsCol("group_ion_mz", s))
            ann[set == s, c("adduct", "ion_mz") := mget(cols)]
            ann[, (cols) := NULL]
        }
        setnames(ann, "group_neutralMass", "neutralMass")
        setcolorder(ann, c("set", "group", "adduct", "neutralMass", "ion_mz"))
        do.call(featureGroupsSet, c(constArgs, list(annotations = ann, algorithm = "table-set", groupAlgo = groupAlgo,
                                                    groupArgs = if (is.null(groupArgs)) list() else groupArgs)))
    }
    else
        do.call(featureGroupsTable, constArgs)
    
    printf("Imported %d feature groups.\n", nrow(gInfo))
    
    return(ret)
}
