# SPDX-FileCopyrightText: 2016 - 2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
#' @include features.R
NULL

#' @rdname importFeaturesTable
featuresTable <- setClass("featuresTable", contains = "features")

setMethod("initialize", "featuresTable", function(.Object, ...) callNextMethod(.Object, algorithm = "table", ...))

#' @export
importFeaturesTable <- function(analysisInfo, input)
{
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo)
    checkmate::assert(
        checkmate::checkDataFrame(input),
        checkmate::checkFileExists(input, access = "r"),
        .var.name = "input"
    )
    
    if (is.character(input))
        input <- fread(input)
 
    ac <- checkmate::makeAssertCollection()
    checkmate::assert(
        checkmate::checkNumeric(input$ID, any.missing = FALSE),
        checkmate::checkCharacter(input$ID, any.missing = FALSE, min.chars = 1),
        .var.name = "input$ID", add = ac
    )
    for (col in c("ret", "mz", "intensity"))
        assertListVal(input, col, checkmate::assertNumeric, any.missing = FALSE, finite = TRUE, add = ac)
    
    maybeAddCol <- function(col, ...)
    {
        if (is.null(input[[col]]))
            eval(substitute(input[, (col) := .v], list(.v = substitute(...))))
        else
            assertListVal(input, col, checkmate::assertNumeric, any.missing = FALSE, finite = TRUE, add = ac)
        return(input)
    }
    
    maybeAddCol("area", 2.5 * intensity)
    maybeAddCol("retmin", ret - defaultLim("retention", "narrow"))
    maybeAddCol("retmax", ret + defaultLim("retention", "narrow"))
    maybeAddCol("mzmin", mz - defaultLim("mz", "narrow"))
    maybeAddCol("mzmax", mz + defaultLim("retention", "narrow"))
    
    hasMob <- !is.null(input[["mobility"]])
    if (hasMob)
    {
        assertListVal(input, "mobility", checkmate::assertNumeric, any.missing = TRUE, finite = TRUE, add = ac)
        maybeAddCol("mobmin", mobility - defaultLim("mobility", "narrow"))
        maybeAddCol("mobmax", mobility + defaultLim("mobility", "narrow"))
        maybeAddCol("mob_area", NA_real_)
        maybeAddCol("mob_intensity", NA_real_)
        if (is.null(input[["ims_parent_ID"]]))
            input[, ims_parent_ID := NA_character_]
        if (is.null(input[["ims_parent_ID"]]))
            input[, mob_assign_method := NA_character_]
    }
    else
    {
        mobCols <- intersect(c("mobmin", "mobmax", "mob_area", "mob_intensity", "ims_parent_ID", "mob_assign_method"),
                             names(input))
        if (length(mobCols) > 0)
        {
            warning("The following columns are invalid in non-IMS workflows and will be removed: ",
                    paste(mobCols, collapse = ", "), call. = FALSE)
            input[, (mobCols) := NULL]
        }
    }
    
    setsCols <- c("set", "adduct", "ion_mz")
    hasSets <- all(setsCols %chin% names(input))
    if (hasSets)
    {
        assertListVal(input, "set", checkmate::assertCharacter, any.missing = FALSE, min.chars = 1, add = ac)
        assertListVal(input, "adduct", checkmate::assertCharacter, any.missing = FALSE, min.chars = 1, add = ac)
        assertListVal(input, "ion_mz", checkmate::assertNumeric, any.missing = FALSE, finite = TRUE, add = ac)
    }
    else
    {
        setsColsPresent <- intersect(setsCols, names(input))
        if (length(setsColsPresent) > 0)
        {
            warning(sprintf("In sets workflows the %s column(s) must be present, but only %s is found. ", setsCols,
                            setsColsPresent),
                    "Assuming this is not a sets workflow and removing interfering columns...", call. = FALSE)
            input[, (setsColsPresent) := NULL]
        }
    }
    
    checkmate::reportAssertions(ac)

    input[, ID := as.character(ID)]
    
    # sync common analyses between analysisInfo
    anasCommon <- intersect(analysisInfo$analysis, input$analysis)
    if (length(anasCommon) == 0)
        stop("No common analyses found between analysisInfo and input data.", call. = FALSE)
    
    anasAIOnly <- setdiff(analysisInfo$analysis, input$analysis)
    if (length(anasAIOnly) > 0)
        warning("The following analyses in analysisInfo are not present in the input data: ",
                paste(anasAIOnly, collapse = ", "), call. = FALSE)
    
    anasInputOnly <- setdiff(input$analysis, analysisInfo$analysis)
    if (length(anasInputOnly) > 0)
        warning("The following analyses in the input data are not present in analysisInfo: ",
                paste(anasInputOnly, collapse = ", "), call. = FALSE)
    
    analysisInfo <- analysisInfo[analysis %in% anasCommon]

    if (!is.null(analysisInfo[["set"]]))
    {
        warning("The 'set' column in the analysisInfo is not used and will be removed")
        analysisInfo[, set := NULL]
    }
    
    if (hasSets)
    {
        # move set assignments to anaInfo
        analysisInfo[, set := input$set[match(analysis, input$analysis)]]
        input[, set := NULL]
    }
    
    if (!is.null(input[["group"]]))
    {
        warning("The 'group' column conflicts internally and will be removed.", call. = FALSE)
        input[, group := NULL]
    }
    normCols <- intersect(c("intensity_rel", "area_rel"), names(input))
    if (length(normCols) > 0)
    {
        warning("The following columns are used internally for normalization and will be removed.",
                paste(normCols, collapse = ", "), call. = FALSE)
        input[, (normCols) := NULL]
    }
    
    setcolorder(input, c("ID", "ret", "retmin", "retmax", "mz", "mzmin", "mzmax", "mobility", "mobmin", "mobmax",
                         "intensity", "area"), skip_absent = TRUE)
    
    fTable <- split(input, by = "analysis", keep.by = FALSE)
    fTable <- fTable[anasCommon] # subset and sync to anaInfo order
    
    constArgs <- list(analysisInfo = analysisInfo, features = fTable, hasMobilities = hasMob)
    if (hasSets)
        return(do.call(featuresSet, c(constArgs, list(algorithm = "table-set"))))
    return(do.call(featuresTable, constArgs))
}
