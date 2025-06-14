# SPDX-FileCopyrightText: 2016 - 2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

checkHasNames <- function(x, n, subset = FALSE, type = "unique")
{
    if (subset)
        return(checkmate::checkNames(names(x), subset.of = n, type = type))
    return(checkmate::checkNames(names(x), must.include = n, type = type))
}
assertHasNames <- checkmate::makeAssertionFunction(checkHasNames)

checkRange <- function(x, null.ok = FALSE)
{
    ret <- checkmate::checkNumeric(x, any.missing = FALSE, lower = 0, len = 2, null.ok = null.ok)
    if (isTRUE(ret) && !is.null(x) && x[1] > x[2])
        ret <- paste0("lower range (", x[1], ") higher than upper (", x[2], ")")
    return(ret)
}
assertRange <- checkmate::makeAssertionFunction(checkRange)

assertScoreRange <- function(x, scNames, .var.name = checkmate::vname(x), add = NULL)
{
    checkmate::assertList(x, null.ok = TRUE, types = "numeric", .var.name = .var.name, add = add)
    if (!is.null(x))
    {
        checkmate::assertNames(names(x), type = "unique", subset.of = scNames, .var.name = .var.name, add = add)
        checkmate::qassertr(x, "N2", .var.name = .var.name)
    }
}

checkS4 <- function(x, null.ok = FALSE)
{
    if (is.null(x))
    {
        if (null.ok)
            return(TRUE)
        return("object is NULL")
    }
    if (!isS4(x))
        return("object is not an S4 object")
    return(TRUE)
}
assertS4 <- checkmate::makeAssertionFunction(checkS4)

checkChoiceSilent <- function(x, ch)
{
    ret <- checkmate::checkString(x, min.chars = 1)
    if (isTRUE(ret) && !x %in% ch)
        ret <- paste("Must be element of", getStrListWithMax(ch, 6, ", "))
    return(ret)
}
assertChoiceSilent <- checkmate::makeAssertionFunction(checkChoiceSilent)

checkUniqueDTBy <- function(x, by)
{
    ret <- checkmate::checkDataTable(x)
    if (isTRUE(ret))
        ret <- checkHasNames(x, by)
    if (isTRUE(ret) && anyDuplicated(x, by = by))
        ret <- sprintf("Found non-unique value pairs in columns %s", paste0(by, collapse = "+"))
    return(ret)
}
assertUniqueDTBy <- checkmate::makeAssertionFunction(checkUniqueDTBy)

assertListVal <- function(x, field, assertFunc, mustExist = TRUE, ..., .var.name = checkmate::vname(x))
{
    if (!field %in% names(x))
    {
        if (mustExist)
            stop(sprintf("Field '%s' is missing from %s.", field, .var.name), call. = FALSE)
        return(invisible(NULL))
    }
    assertFunc(x[[field]], ..., .var.name = sprintf("%s[[\"%s\"]]", .var.name, field))
}

assertCharOrFactor <- function(x, empty.ok = FALSE, null.ok = FALSE, ..., .var.name = .var.name)
{
    checkmate::assert(
        checkmate::checkFactor(x, empty.levels.ok = empty.ok, any.missing = empty.ok, null.ok = null.ok),
        checkmate::checkCharacter(x, min.chars = if (empty.ok) 0 else 1, any.missing = empty.ok, null.ok = null.ok),
        .var.name = .var.name
    )
}

assertMSFileType <- function(x, null.ok = FALSE, .var.name = checkmate::vname(x), add = NULL)
{
    checkmate::assertChoice(x, getMSFileTypes(), null.ok = null.ok, .var.name = .var.name, add = add)
}

assertMSConvAlgo <- function(x, .var.name = checkmate::vname(x), add = NULL)
{
    checkmate::assertChoice(x, c("pwiz", "openms", "bruker", "im_collapse", "timsconvert"), .var.name = .var.name,
                            add = add)
}

assertAnalysisInfo <- function(x, fileTypes = NULL, allowedFormats = NULL, null.ok = FALSE,
                               .var.name = checkmate::vname(x), add = NULL)
{
    if (is.null(x) && null.ok)
        return(TRUE)

    if (!is.null(add))
        mc <- length(add$getMessages())

    checkmate::assertDataFrame(x, min.rows = 1, .var.name = .var.name, add = add)
    assertHasNames(x, c("analysis", "replicate", "blank"), .var.name = .var.name, add = add)
    
    assertListVal(x, "analysis", checkmate::assertCharacter, any.missing = FALSE, .var.name = .var.name, add = add)
    assertListVal(x, "replicate", checkmate::assertCharacter, any.missing = FALSE, .var.name = .var.name, add = add)
    assertListVal(x, "blank", checkmate::assertCharacter, any.missing = TRUE, .var.name = .var.name, add = add)
    
    pathCols <- paste0("path_", getMSFileTypes())
    if (!any(pathCols %in% names(x)))
        stop(sprintf("%s does not contain any file paths! Please include at least one of: %s", .var.name,
                     paste0(pathCols, collapse = ", ")), call. = FALSE)
    for (col in pathCols)
    {
        if (!is.null(x[[col]]))
            assertListVal(x, col, checkmate::assertCharacter, any.missing = TRUE, .var.name = .var.name, add = add)
    }

    checkmate::assert(
        checkmate::checkNull(x[["norm_conc"]]),
        checkmate::checkCharacter(x[["norm_conc"]]),
        checkmate::checkNumeric(x[["norm_conc"]]),
        .var.name = sprintf("%s[[\"norm_conc\"]]", .var.name)
    )
    
    # only continue if previous assertions didn't fail: x needs to be used as list which otherwise gives error
    # NOTE: this is only applicable if add != NULL, otherwise previous assertions will throw errors
    if (is.null(add) || length(add$getMessages()) == mc)
    {
        if (!is.null(fileTypes))
        {
            if (is.null(allowedFormats))
                allowedFormats <- lapply(fileTypes, getMSFileFormats)
            
            # stops if files are missing
            getMSFilesFromAnaInfo(x, fileTypes, allowedFormats, mustExist = TRUE)
        }

        checkmate::assertVector(x$analysis, unique = TRUE, .var.name = paste0(.var.name, "$analysis"), add = add)
    }

    invisible(NULL)
}

assertAndPrepareAnaInfo <- function(x, ..., null.ok = FALSE, .var.name = checkmate::vname(x), add = NULL)
{
    if (is.null(x) && null.ok)
        return(NULL)
    
    if (!is.null(add))
        mc <- length(add$getMessages())

    checkmate::assertDataFrame(x, .var.name = .var.name)

    x <- makeDT(x) # convert to DT or make a unique copy
    x <- unFactorDT(x)

    if (is.null(x[["path_centroid"]]) && !is.null(x[["path"]]))
    {
        warning("The usage of the 'path' column in the analysis information is deprecated. ",
                "The column will be renamed to 'path_centroid', but this may not be what you want...", call. = FALSE)
        setnames(x, "path", "path_centroid")
    }
    if (is.null(x[["replicate"]]) && !is.null(x[["group"]]))
    {
        warning("The usage of the 'group' column in the analysis information is deprecated. Please rename this column to 'replicate'.",
                call. = FALSE)
        setnames(x, "group", "replicate")
    }
    
    assertAnalysisInfo(x, ..., .var.name = .var.name, add = add)

    if ((is.null(add) || length(add$getMessages()) == mc) && !is.null(x))
    {
        chrCols <- c("analysis", "replicate", "blank", getAnaInfoPathCols(x))
        for (col in chrCols)
            x[, (col) := as.character(get(col))]
        if (!is.null(x[["norm_conc"]]))
            x[, norm_conc := as.numeric(norm_conc)]
        x[is.na(blank), blank := ""]
    }
    
    return(x)
}

assertGroupFeatAlgo <- function(x, .var.name = checkmate::vname(x), add = NULL)
{
    checkmate::assertChoice(x, c("openms", "xcms", "xcms3", "kpic2"), .var.name = .var.name, add = add)
}

assertGroupFeatVerbose <- function(x, .var.name = checkmate::vname(x), add = NULL)
{
    checkmate::assert(
        checkmate::checkFlag(x),
        checkmate::checkChoice(x, "full"),
        .var.name = .var.name, add = add
    )
}

assertAnaInfoBy <- function(x, anaInfo, withFGroups, null.ok = FALSE, any.missing = FALSE,
                            .var.name = checkmate::vname(x), add = NULL)
{
    if (is.null(x) && null.ok)
        return(NULL)
    
    ch <- names(anaInfo)
    if (withFGroups)
        ch <- c(ch, "fGroups")
    checkmate::assertChoice(x, ch, .var.name = .var.name) # no add: let it fail
    if ((!withFGroups || x != "fGroups") && !any.missing)
        checkmate::assertAtomicVector(anaInfo[[x]], .var.name = sprintf("anaInfo[[\"%s\"]]", x),
                                      any.missing = any.missing, add = add)
    
}

assertAndPrepareAnaInfoBy <- function(x, anaInfo, withFGroups, null.ok = FALSE, any.missing = FALSE,
                                      .var.name = checkmate::vname(x), add = NULL)
{
    if (is.null(x) && null.ok)
        return(NULL)

    if (isTRUE(x))
        x <- "replicate"
    else if (isFALSE(x))
        x <- "analysis"
    
    assertAnaInfoBy(x, anaInfo, withFGroups, null.ok, any.missing, .var.name = .var.name, add = add)
    
    return(x)
}

assertSuspectList <- function(x, needsAdduct, skipInvalid, null.ok = FALSE, .var.name = checkmate::vname(x), add = NULL)
{
    # NOTE: we don't check mobility/CCS columns here, since we don't know yet which adduct specific column to take,
    # things would become complex. Instead, all checks will be done in expandSuspMobilities()
    
    if (is.null(x) && null.ok)
        return(invisible(NULL))
    
    mzCols <- c("mz", "neutralMass", "SMILES", "InChI", "formula")
    
    # this seems necessary for proper naming in subsequent assertions (why??)
    .var.name <- force(.var.name)
    
    checkmate::assertDataFrame(x, any.missing = TRUE, min.rows = 1, .var.name = .var.name, add = add)
    assertHasNames(x, "name", .var.name = .var.name, add = add)

    checkmate::assertNames(intersect(names(x), mzCols), subset.of = mzCols,
                           .var.name = paste0("names(", .var.name, ")"), add = add)

    needsAdduct <- needsAdduct && (is.null(x[["mz"]]) || any(is.na(x$mz)))
    if (needsAdduct)
    {
        msg <- "Adduct information is required to calculate ionized suspect masses. "
        
        if (is.null(x[["adduct"]]))
            stop(msg, "Please either set the adduct argument or add an adduct column in the suspect list.")
        if (any(is.na(x[["adduct"]]) & (is.null(x[["mz"]]) | is.na(x[["mz"]]))))
            stop(msg, "Please either set the adduct argument or make sure that suspects without mz information have data in the adduct column.")
    }
    
    for (col in c("name", "SMILES", "InChI", "formula", "InChIKey", "adduct", "fragments_mz", "fragments_formula"))
    {
        emptyOK <- col != "name" && (col != "adduct" || !needsAdduct)
        assertListVal(x, col, assertCharOrFactor, empty.ok = emptyOK, mustExist = !emptyOK, .var.name = .var.name,
                      add = add)
    }

    for (col in c("mz", "neutralMass", "rt"))
        assertListVal(x, col, checkmate::assertNumeric, any.missing = TRUE, mustExist = FALSE,
                      lower = if (col != "rt") 0 else -Inf, finite = TRUE, .var.name = .var.name, add = add)

    if (!skipInvalid)
    {
        cx <- if (is.data.table(x)) copy(x) else as.data.table(x)
        cx[, OK := any(!sapply(.SD, is.na)), by = seq_len(nrow(cx)), .SDcols = intersect(names(x), mzCols)]
        if (all(!cx$OK))
            stop("Suspect list does not contain any (data to calculate) suspect masses", call. = FALSE)
        else if (any(!cx$OK))
            stop("Suspect list does not contain any (data to calculate) suspect masses for row(s): ",
                 paste0(which(!cx$OK), collapse = ", "), call. = FALSE)
    }
    
    invisible(NULL)
}

assertAndPrepareSuspectsSets <- function(x, sets, skipInvalid, .var.name = checkmate::vname(x))
{
    if (checkmate::testDataFrame(x))
    {
        assertSuspectList(x, FALSE, skipInvalid, .var.name = .var.name)
        if (length(sets) > 1)
        {
            cols <- c("mz", "adduct", "fragments_mz")
            for (cl in cols)
            {
                if (!is.null(x[[cl]]) && !all(is.na(x[[cl]])))
                {
                    warning("The suspect list seems to contain an mz, adduct or fragments_mz column, ",
                            "which are generally specific to the ionization mode used. ",
                            "These columns most likely need to removed since the same suspect list will be used for all sets.",
                            call. = FALSE)
                    break
                }
            }
        }
        x <- sapply(sets, function(s) x, simplify = FALSE) # same for all sets
    }
    else
    {
        checkmate::assertList(x, "data.frame", any.missing = FALSE, all.missing = FALSE,
                              len = length(sets))
        checkmate::assert(
            checkmate::checkNames(names(x), "unnamed"),
            checkmate::checkNames(names(x), "unique", must.include = sets),
            .var.name = .var.name
        )
        if (checkmate::testNames(names(x), "unnamed"))
            names(x) <- sets
    }
    
    # sync order
    x <- x[sets]
    
    return(x)
}

assertLogicTransformations <- function(x, null.ok = FALSE, .var.name = checkmate::vname(x), add = NULL)
{
    if (null.ok && is.null(x))
        return(NULL)
    
    checkmate::assertDataFrame(x, min.rows = 1, any.missing = FALSE, col.names = "unique", .var.name = .var.name,
                               add = add)
    checkmate::assertNames(colnames(x), permutation.of = c("transformation", "add", "sub", "retDir"), what = "colnames",
                           add = add)

    assertListVal(x, "transformation", checkmate::assertCharacter, min.chars = 1, any.missing = FALSE, unique = TRUE,
                  .var.name = .var.name, add = add)
    assertListVal(x, "add", checkmate::assertCharacter, .var.name = .var.name, add = add)
    assertListVal(x, "sub", checkmate::assertCharacter, .var.name = .var.name, add = add)
    assertListVal(x, "retDir", checkmate::assertSubset, choices = c(-1, 0, 1), .var.name = .var.name, add = add)
}

assertCanCreateDir <- function(x, .var.name = checkmate::vname(x), add = NULL)
{
    if (!is.null(add))
        mc <- length(add$getMessages())

    checkmate::assertString(x, min.chars = 1, .var.name = .var.name, add = add)

    # only continue if previous assertions didn't fail: x needs to be a valid path for next assertions
    # NOTE: this is only applicable if add != NULL, otherwise previous assertions will throw errors
    if (is.null(add) || length(add$getMessages()) == mc)
    {
        # find first existing directory and see if it's writable
        x <- normalizePath(x, mustWork = FALSE)
        repeat
        {
            if (file.exists(x))
            {
                checkmate::assertDirectoryExists(x, "w", .var.name = .var.name, add = add)
                break
            }
            x <- normalizePath(dirname(x), mustWork = FALSE)
        }
    }
    invisible(NULL)
}

assertCanCreateDirs <- function(x, .var.name = checkmate::vname(x), add = NULL)
{
    for (ana in x)
        assertCanCreateDir(ana, .var.name, add)
}

assertDirExists <- function(x, access, null.ok = FALSE, .var.name = checkmate::vname(x), add = NULL)
{
    if (is.null(x) && null.ok)
        return(invisible(NULL))
    checkmate::assertDirectoryExists(x, access, .var.name = .var.name, add = add)
}

assertDACloseSaveArgs <- function(x, save, .var.name = checkmate::vname(x), add = NULL)
{
    checkmate::assertFlag(x, .var.name = .var.name, add = add)
    checkmate::assertFlag(save, .var.name = "save", add = add)
}

assertXYLim <- function(x, ylim, .var.name = checkmate::vname(x), add = NULL)
{
    checkmate::assertNumeric(x, finite = TRUE, .var.name = .var.name, len = 2, null.ok = TRUE, add = add)
    checkmate::assertNumeric(ylim, finite = TRUE, .var.name = "ylim", len = 2, null.ok = TRUE, add = add)
}

assertConsCommonArgs <- function(absMinAbundance, relMinAbundance, uniqueFrom, uniqueOuter, objNames, add = NULL)
{
    checkmate::assertNumber(absMinAbundance, .var.name = "absMinAbundance", null.ok = TRUE, add = add)
    checkmate::assertNumber(relMinAbundance, .var.name = "relMinAbundance", null.ok = TRUE, add = add)
    checkmate::assert(checkmate::checkLogical(uniqueFrom, min.len = 1, max.len = length(objNames), any.missing = FALSE, null.ok = TRUE),
                      checkmate::checkIntegerish(uniqueFrom, lower = 1, upper = length(objNames), any.missing = FALSE),
                      checkmate::checkSubset(uniqueFrom, objNames, empty.ok = FALSE),
                      .var.name = "uniqueFrom")
    checkmate::assertFlag(uniqueOuter, .var.name = "uniqueOuter", add = add)

    if (!is.null(uniqueFrom) && (!is.null(absMinAbundance) || !is.null(relMinAbundance)))
        stop("Cannot apply both unique and abundance filters simultaneously.")
}

checkCSVFile <- function(x, cols)
{
    ret <- checkmate::checkFileExists(x, "r")
    if (isTRUE(ret))
    {
        t <- fread(x, nrows = 1) # nrows=0 doesn't always bug (may give internal error)
        missingc <- setdiff(cols, names(t))
        if (length(missingc) > 0)
            ret <- paste0("Missing columns: ", paste0(missingc, collapse = ", "))
    }
    return(ret)
}
assertCSVFile <- checkmate::makeAssertionFunction(checkCSVFile)

# used for "[" methods
checkSubsetArg <- function(x)
{
    ret <- checkmate::checkIntegerish(x)
    if (!isTRUE(ret))
        ret <- checkmate::checkCharacter(x)
    if (!isTRUE(ret))
        ret <- checkmate::checkLogical(x)
    if (!isTRUE(ret))
        ret <- "Should be valid numeric, character or logical"
    return(ret)
}
assertSubsetArg <- checkmate::makeAssertionFunction(checkSubsetArg)

assertSubsetArgAndToChr <- function(x, choices, .var.name = checkmate::vname(x), add = NULL)
{
    assertSubsetArg(x, .var.name = .var.name, add = add)
    if (!is.character(x))
        x <- choices[x]
    x <- intersect(x, choices)
    return(x)
}

# used for "[[" methods
checkExtractArg <- function(x)
{
    ret <- checkmate::checkInt(x, lower = 0)
    if (!isTRUE(ret))
        ret <- checkmate::checkString(x)
    if (!isTRUE(ret))
        ret <- "Should be valid numeric or character scalar"
    return(ret)
}
assertExtractArg <- checkmate::makeAssertionFunction(checkExtractArg)

checkDeleteArg <- function(x)
{
    ret <- checkmate::checkNull(x)
    if (!isTRUE(ret))
        ret <- checkmate::checkIntegerish(x, any.missing = FALSE)
    if (!isTRUE(ret))
        ret <- checkmate::checkCharacter(x, any.missing = FALSE)
    if (!isTRUE(ret))
        ret <- checkmate::checkLogical(x, any.missing = FALSE)
    if (!isTRUE(ret))
        ret <- "Should be NULL, valid numeric, character or logical"
    return(ret)
}
assertDeleteArg <- checkmate::makeAssertionFunction(checkDeleteArg)

assertDeleteArgAndToChr <- function(x, choices, .var.name = checkmate::vname(x), add = NULL)
{
    if (!is.null(add))
        mc <- length(add$getMessages())
    
    assertDeleteArg(x, .var.name = .var.name, add = add)
    
    if (!is.null(add) && length(add$getMessages()) != mc)
        return(x) # assert failed

    if (is.null(x))
        x <- choices
    else
    {
        if (!is.character(x))
            x <- choices[x]
        x <- intersect(x, choices)
    }
    
    return(x)
}

assertIMSArg <- function(x, null.ok = FALSE, .var.name = checkmate::vname(x), add = NULL)
{
    if (is.null(x) && null.ok)
        return(invisible(NULL))

    checkmate::assert(
        checkmate::checkFlag(x),
        checkmate::checkChoice(x, c("both", "maybe")),
        .var.name = .var.name, add = add
    )
}

assertApplyIMSArg <- function(x, .var.name = checkmate::vname(x), add = NULL)
{
    checkmate::assert(
        checkmate::checkFlag(x),
        checkmate::checkChoice(x, "both"),
        .var.name = .var.name, add = add
    )
}

assertFGAsDataTableArgs <- function(fGroups, areas, features, qualities, regression, regressionBy, averageFunc,
                                    normalized, FCParams, concAggrParams, toxAggrParams, normConcToTox, anaInfoCols,
                                    IMS, collapseSuspects, onlyHits)
{
    anaInfo <- analysisInfo(fGroups)
    
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertFlag, . ~ areas + features + normalized + normConcToTox,
           fixed = list(add = ac))
    checkmate::assert(
        checkmate::checkFALSE(qualities),
        checkmate::checkChoice(qualities, c("quality", "score", "both")),
        .var.name = "qualities", add = ac
    )
    if (!isFALSE(regression))
    {
        assertAnaInfoBy(regression, anaInfo, FALSE, any.missing = TRUE, add = ac)
        checkmate::assertNumeric(anaInfo[[regression]],
                                 .var.name = sprintf("analysisInfo(fGroups)[[\"%s\"]]", regression), add = ac)
    }
    checkmate::assertChoice(regressionBy, names(anaInfo), null.ok = TRUE, add = ac)
    checkmate::assertFunction(averageFunc, add = ac)
    assertFCParams(FCParams, fGroups, null.ok = TRUE, add = ac)
    aapply(assertPredAggrParams, . ~ concAggrParams + toxAggrParams, null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertSubset(anaInfoCols, names(anaInfo), empty.ok = TRUE, add = ac)
    assertIMSArg(IMS, add = ac)
    checkmate::assertString(collapseSuspects, null.ok = TRUE, add = ac)
    checkmate::assertFlag(onlyHits, add = ac)
    checkmate::reportAssertions(ac)
}

assertNormalizationMethod <- function(x, withNone = TRUE, .var.name = checkmate::vname(x), add = NULL)
{
    ch <- c("max", "minmax")
    if (withNone)
        ch <- c(ch, "none")
    checkmate::assertChoice(x, ch, .var.name = .var.name, add = add)
}

assertClusterMethod <- function(x, .var.name = checkmate::vname(x), add = NULL)
{
    checkmate::assertChoice(x, c("bin", "distance", "hclust"), .var.name = .var.name, add = add)
}

assertEIXParams <- function(x, .var.name = checkmate::vname(x), add = NULL)
{
    checkmate::assertList(x, names = "unique", .var.name = .var.name) # no add: should fail
    
    assertListVal(x, "window", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name, add = add)
    assertListVal(x, "mzExpWindow", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name, add = add)
    assertListVal(x, "mzExpIMSWindow", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name, add = add)
    assertListVal(x, "minIntensityIMS", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name, add = add)
    assertListVal(x, "mobExpWindow", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name, add = add)
    assertListVal(x, "topMost", checkmate::assertCount, positive = TRUE, null.ok = TRUE, .var.name = .var.name,
                  add = add)
    assertListVal(x, "topMostByReplicate", checkmate::assertFlag, .var.name = .var.name, add = add)
    assertListVal(x, "onlyPresent", checkmate::assertFlag, .var.name = .var.name, add = add)
    assertListVal(x, "setsAdductPos", checkAndToAdduct, .var.name = .var.name)
    assertListVal(x, "setsAdductNeg", checkAndToAdduct, .var.name = .var.name)
    
    if (!is.null(x[["topMost"]]) && !isTRUE(x$onlyPresent))
        stop("onlyPresent must be TRUE if topMost is set", call. = FALSE)
    
    invisible(NULL)
}

assertEICParams <- function(x, .var.name = checkmate::vname(x), add = NULL)
{
    assertEIXParams(x, .var.name = .var.name, add = add)
    invisible(NULL)
}

assertEIMParams <- function(x, .var.name = checkmate::vname(x), add = NULL)
{
    assertEIXParams(x, .var.name = .var.name, add = add)
    
    assertListVal(x, "maxRTWindow", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name, add = add)
    assertListVal(x, "IMSWindow", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name, add = add)
    assertListVal(x, "clusterMethod", assertClusterMethod, .var.name = .var.name, add = add)
    
    invisible(NULL)
}

assertFCParams <- function(x, fGroups, null.ok = FALSE, .var.name = checkmate::vname(x), add = NULL)
{
    if (null.ok && is.null(x))
        return(NULL)
    
    checkmate::assertList(x, names = "unique", .var.name = .var.name) # no add: should fail
    
    assertListVal(x, "replicates", checkmate::assertCharacter, any.missing = FALSE, len = 2, .var.name = .var.name,
                  add = add)
    assertListVal(x, "replicates", checkmate::assertSubset, choices = replicates(fGroups), .var.name = .var.name,
                  add = add)
    assertListVal(x, "thresholdFC", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name, add = add)
    assertListVal(x, "thresholdPV", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name, add = add)
    assertListVal(x, "zeroValue", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name, add = add)
    assertListVal(x, "zeroMethod", checkmate::assertChoice, choices = c("add", "fixed", "omit"), .var.name = .var.name,
                  add = add)
    assertListVal(x, "PVTestFunc", checkmate::assertFunction, .var.name = .var.name, add = add)
    assertListVal(x, "PVAdjFunc", checkmate::assertFunction, .var.name = .var.name, add = add)
}

assertPredAggrParams <- function(x, null.ok = FALSE, .var.name = checkmate::vname(x), add = NULL)
{
    if (null.ok && is.null(x))
        return(NULL)
    
    checkmate::assertList(x, names = "unique", .var.name = .var.name) # no add: should fail
    
    assertListVal(x, "typeFunc", checkmate::assertFunction, .var.name = .var.name, add = add)
    assertListVal(x, "groupFunc", checkmate::assertFunction, .var.name = .var.name, add = add)
    assertListVal(x, "candidateFunc", checkmate::assertFunction, .var.name = .var.name, add = add)
    assertListVal(x, "preferType", checkmate::assertChoice, choices = c("suspect", "compound", "SIRIUS_FP", "none"),
                  .var.name = .var.name, add = add)
}

assertAvgPListParams <- function(x, .var.name = checkmate::vname(x), add = NULL)
{
    checkmate::assertList(x, names = "unique", .var.name = .var.name) # no add: should fail

    assertListVal(x, "clusterMzWindow", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name,
                  add = add)
    assertListVal(x, "minAbundanceRel", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name,
                  add = add)
    assertListVal(x, "minAbundanceAbs", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name,
                  add = add)
    assertListVal(x, "minAbundanceIMSRel", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name,
                  add = add)
    assertListVal(x, "minAbundanceIMSAbs", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name,
                  add = add)
    assertListVal(x, "topMost", checkmate::assertCount, positive = TRUE, .var.name = .var.name, add = add)
    assertListVal(x, "minIntensityPre", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name,
                  add = add)
    assertListVal(x, "minIntensityPost", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name,
                  add = add)
    assertListVal(x, "minIntensityIMS", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name,
                  add = add)
    assertListVal(x, "method", assertClusterMethod, .var.name = .var.name, add = add)
    assertListVal(x, "withPrecursorMS", checkmate::assertFlag, .var.name = .var.name, add = add)
    assertListVal(x, "pruneMissingPrecursorMS", checkmate::assertFlag, .var.name = .var.name, add = add)
    assertListVal(x, "retainPrecursorMSMS", checkmate::assertFlag, .var.name = .var.name, add = add)
}

assertPListIsolatePrecParams <- function(x, .var.name = checkmate::vname(x), add = NULL)
{
    if (is.null(x))
        return(NULL)

    checkmate::assertList(x, names = "unique", .var.name = .var.name) # no add: should fail

    assertListVal(x, "maxIsotopes", checkmate::assertCount, .var.name = .var.name, add = add)
    assertListVal(x, "mzDefectRange", checkmate::assertNumeric, any.missing = FALSE, len = 2, finite = TRUE,
                  .var.name = .var.name, add = add)
    assertListVal(x, "intRange", checkmate::assertNumeric, any.missing = FALSE, len = 2, finite = TRUE,
                  .var.name = .var.name, add = add)
    assertListVal(x, "z", checkmate::assertCount, positive = TRUE, .var.name = .var.name, add = add)
    assertListVal(x, "maxGap", checkmate::assertCount, positive = TRUE, .var.name = .var.name, add = add)
}

assertSpecSimParams <- function(x, .var.name = checkmate::vname(x), add = NULL)
{
    checkmate::assertList(x, names = "unique", .var.name = .var.name) # no add: should fail

    assertListVal(x, "method", checkmate::assertChoice, choices = c("cosine", "jaccard"), .var.name = .var.name,
                  add = add)
    assertListVal(x, "removePrecursor", checkmate::assertFlag, .var.name = .var.name, add = add)
    assertListVal(x, "mzWeight", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name, add = add)
    assertListVal(x, "intWeight", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name, add = add)
    assertListVal(x, "absMzDev", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name, add = add)
    assertListVal(x, "relMinIntensity", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name,
                  add = add)
    assertListVal(x, "minPeaks", checkmate::assertCount, positive = TRUE, .var.name = .var.name, add = add)
    assertListVal(x, "shift", checkmate::assertChoice, choices = c("none", "precursor", "both"), .var.name = .var.name,
                  add = add)
    assertListVal(x, "setCombineMethod", checkmate::assertChoice, choices = c("mean", "min", "max"),
                  .var.name = .var.name, add = add)
}

assertTPStructParams <- function(x, .var.name = checkmate::vname(x), add = NULL)
{
    checkmate::assertList(x, names = "unique", .var.name = .var.name) # no add: should fail

    assertListVal(x, "calcLogP", checkmate::assertChoice, choices = c("rcdk", "obabel", "none"), .var.name = .var.name,
                  add = add)
    assertListVal(x, "forceCalcLogP", checkmate::assertFlag, .var.name = .var.name, add = add)
    assertListVal(x, "forceCalcRetDir", checkmate::assertFlag, .var.name = .var.name, add = add)
    assertListVal(x, "minLogPDiff", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name, add = add)
    assertListVal(x, "calcSims", checkmate::assertFlag, .var.name = .var.name, add = add)
    assertListVal(x, "fpType", checkmate::assertString, min.chars = 1, .var.name = .var.name, add = add)
    assertListVal(x, "fpSimMethod", checkmate::assertString, min.chars = 1, .var.name = .var.name, add = add)
}

assertCheckSession <- function(x, mustExist, null.ok = FALSE, .var.name = checkmate::vname(x), add = NULL)
{
    if (null.ok && is.null(x))
        return(NULL)
    
    checkmate::assertString(x, min.chars = 1, .var.name = .var.name, add = add)
    if (mustExist)
        checkmate::assertFileExists(x, "r", .var.name = .var.name, add = add)
    else
        checkmate::assertPathForOutput(x, overwrite = TRUE, .var.name = .var.name, add = add)
    
    # UNDONE: validate YAML?
}

checkSetLabels <- function(x, len)
{
    ret <- checkmate::checkCharacter(x, min.chars = 1, any.missing = FALSE, len = len, unique = TRUE)
    if (isTRUE(ret) && any(grepl("\\W", x)))
        ret <- "Set labels should only contain alphanumeric characters and underscores"
    if (isTRUE(ret) && any(grepl("genform|sirius|bruker|metfrag|library", x)))
        ret <- "Set labels cannot contain annotation algorithm names"
    return(ret)
}
assertSetLabels <- checkmate::makeAssertionFunction(checkSetLabels)

assertSets <- function(obj, s, multiple, null.ok = multiple, .var.name = checkmate::vname(s), add = NULL)
{
    if (multiple)
        checkmate::assertSubset(s, sets(obj), empty.ok = null.ok, .var.name = .var.name, add = add)
    else
        checkmate::assertChoice(s, sets(obj), null.ok = null.ok, .var.name = .var.name, add = add)
}

assertMakeSetArgs <- function(objects, class, adducts, adductNullOK, labels, add = NULL)
{
    checkmate::assertList(objects, types = class, any.missing = FALSE,
                          unique = TRUE, .var.name = "obj and ...", min.len = 1,
                          add = add)
    if (!adductNullOK || !is.null(adducts))
        checkmate::assert(checkmate::checkCharacter(adducts, any.missing = FALSE, min.len = 1,
                                                    max.len = length(objects)),
                          checkmate::checkList(adducts, types = c("adduct", "character"), any.missing = FALSE,
                                               min.len = 1, max.len = length(objects)),
                          .var.name = "adducts")
    checkmate::assertCharacter(labels, len = length(objects), min.chars = 1, unique = TRUE,
                               null.ok = !is.null(adducts), add = add)
}

assertDynamicTreeCutArgs <- function(maxTreeHeight, deepSplit, minModuleSize, add = NULL)
{
    checkmate::assertNumber(maxTreeHeight, 0, finite = TRUE, add = add)
    checkmate::assertFlag(deepSplit, add = add)
    checkmate::assertCount(minModuleSize, positive = TRUE, add = add)
}

assertSIRIUSLogin <- function(x, .var.name = checkmate::vname(x), add = NULL)
{
    checkmate::assert(
        checkmate::checkFALSE(x),
        checkmate::checkChoice(x, c("check", "interactive")),
        checkmate::checkCharacter(x, any.missing = FALSE, min.chars = 1, len = 2),
        .var.name = .var.name, add = add
    )
    if (is.character(x) && length(x) == 2)
        assertHasNames(x, c("username", "password"), .var.name = x, add = add)
}

assertAndPrepareReportSettings <- function(settings, setAggr = TRUE)
{
    emptyListToVec <- function(val, evec)
    {
        # yaml package always returns empty sequences as lists --> convert to empty vector
        return(if (is.list(val) && length(val) == 0) evec else val)
    }
    assertAndToFunc <- function(x, .var.name = checkmate::vname(val))
    {
        checkmate::assertString(x, .var.name = .var.name)
        x <- get(x)
        checkmate::assertFunction(x, .var.name = .var.name)
        return(x)
    }
    
    checkmate::assertList(settings, any.missing = FALSE)
    assertHasNames(settings, c("general", "summary", "features", "MSPeakLists", "formulas", "compounds", "TPs",
                               "internalStandards"))

    ac <- checkmate::makeAssertCollection()
    
    checkmate::assertList(settings$general)
    # NOTE: version 1 wasn't specified, so might be absent
    checkmate::assertCount(settings$general[["version"]], positive = TRUE, null.ok = TRUE) # don't add!
    
    # check version first: if file is older we silently update it so subsequent checks won't fail
    defSettings <- readYAML(system.file("report", "settings.yml", package = "patRoon"))
    if (is.null(settings$general[["version"]]) || settings$general$version < defSettings$general$version)
    {
        warning("Report settings file is older than current and might be incomplete. ",
                "Use genReportSettingsFile() to update the file.", call. = FALSE)
        settings <- adjustReportSettings(defSettings, settings)
    }
    else if (settings$general$version > defSettings$general$version)
        warning("Report settings file is newer than current! ",
                "Update patRoon to support all settings", call. = FALSE)
    
    checkmate::assertSubset(settings$general$format, choices = "html", add = ac)
    checkmate::assertPathForOutput(settings$general$path, overwrite = TRUE, add = ac)
    checkmate::assertCount(settings$general$keepUnusedPlots, positive = FALSE, add = ac)
    checkmate::assertFlag(settings$general$selfContained, add = ac)
    checkmate::assertFlag(settings$general$noDate, add = ac)
    
    checkmate::assertSubset(settings$summary, c("chord", "venn", "upset"), add = ac)
    
    checkmate::assertList(settings$features)
    checkmate::assertFlag(settings$features$retMin, add = ac)
    checkmate::assertList(settings$features$chromatograms)
    checkmate::assertFlag(settings$features$chromatograms$large, add = ac)
    checkmate::assertFlag(settings$features$chromatograms$small, add = ac)
    checkmate::assert(
        checkmate::checkFlag(settings$features$chromatograms$features),
        checkmate::checkChoice(settings$features$chromatograms$features, "all"),
        .var.name = "settings$features$chromatograms$features",
        add = ac
    )
    checkmate::assertChoice(settings$features$chromatograms$intMax, c("eic", "feature"), add = ac)
    checkmate::assertFlag(settings$features$mobilograms$large, add = ac)
    checkmate::assertFlag(settings$features$mobilograms$small, add = ac)
    checkmate::assert(
        checkmate::checkFlag(settings$features$mobilograms$features),
        checkmate::checkChoice(settings$features$mobilograms$features, "all"),
        .var.name = "settings$features$mobilograms$features",
        add = ac
    )
    checkmate::assertFlag(settings$features$intensityPlots, add = ac)
    if (setAggr)
    {
        settings$features$aggregateConcs <- getDefPredAggrParams(assertAndToFunc(settings$features$aggregateConcs))
        settings$features$aggregateTox <- getDefPredAggrParams(assertAndToFunc(settings$features$aggregateTox))
    }
    
    checkmate::assertList(settings$MSPeakLists)
    checkmate::assertFlag(settings$MSPeakLists$spectra, add = ac)
    
    checkmate::assertList(settings$formulas)
    checkmate::assertFlag(settings$formulas$include, add = ac)
    assertNormalizationMethod(settings$formulas$normalizeScores, add = ac)
    settings$formulas$exclNormScores <- emptyListToVec(settings$formulas$exclNormScores, character())
    checkmate::assertCharacter(settings$formulas$exclNormScores, min.chars = 1, any.missing = FALSE, add = ac)
    checkmate::assertCount(settings$formulas$topMost, positive = TRUE, add = ac)
    
    checkmate::assertList(settings$compounds)
    assertNormalizationMethod(settings$compounds$normalizeScores, add = ac)
    settings$compounds$exclNormScores <- emptyListToVec(settings$compounds$exclNormScores, character())
    checkmate::assertCharacter(settings$compounds$exclNormScores, min.chars = 1, any.missing = FALSE, add = ac)
    checkmate::assertCount(settings$compounds$topMost, positive = TRUE, add = ac)
    
    checkmate::assertList(settings$TPs)
    checkmate::assertFlag(settings$TPs$graphs, add = ac)
    checkmate::assertCount(settings$TPs$graphStructuresMax, add = ac)
    
    checkmate::assertList(settings$internalStandards)
    checkmate::assertFlag(settings$internalStandards$graph, add = ac)
    
    checkmate::reportAssertions(ac)
    
    return(settings)
}

assertAndPrepareLimits <- function(limits)
{
    checkmate::assertList(limits, any.missing = FALSE)
    assertHasNames(limits, c("general", "retention", "mz", "mobility_bruker", "mobility_agilent", "CCS"))
    
    assertListVal(limits, "general", checkmate::assertList, any.missing = FALSE)
    assertListVal(limits$general, "version", checkmate::assertCount, positive = TRUE)
    
    # check version first: if file is older we silently update it so subsequent checks won't fail
    # NOTE: this for the future
    
    ac <- checkmate::makeAssertCollection()
    assertListVal(limits$general, "IMS", checkmate::assertChoice, choices = c("bruker", "agilent"), add = ac)

    
    levels <- c("very_narrow", "narrow", "medium", "wide")
    levels <- c(levels, paste0(levels, "_rel"))
    for (categ in setdiff(names(limits), "general"))
    {
        assertListVal(limits, categ, checkmate::assertList, any.missing = FALSE, add = ac)
        assertListVal(limits, categ, assertHasNames, n = levels, subset = TRUE, add = ac)
        for (limit in names(limits[[categ]]))
            assertListVal(limits[[categ]], limit, checkmate::assertNumber, lower = 0, finite = TRUE,
                          .var.name = paste0("limits$", categ), add = ac)
        
    }
    
    checkmate::reportAssertions(ac)
    
    if (limits$general$IMS == "bruker")
    {
        limits$mobility_agilent <- NULL
        names(limits)[which(names(limits) == "mobility_bruker")] <- "mobility"
    }
    else
    {
        limits$mobility_bruker <- NULL
        names(limits)[which(names(limits) == "mobility_agilent")] <- "mobility"
    }
    
    return(limits)
}

checkConcUnit <- function(x)
{
    bases <- c("n", "u", "m", "")
    units <- c(paste0(bases, "gL"), paste0(bases, "M"))
    units <- c(units, paste("log", units), paste("log10", units), paste("log2", units))
    checkmate::checkChoice(x, units)
}
assertConcUnit <- checkmate::makeAssertionFunction(checkConcUnit)

assertAndPrepareQuantCalib <- function(calibration, concUnit)
{
    checkmate::assertDataFrame(calibration, any.missing = FALSE)
    
    calibration <- if (is.data.table(calibration)) copy(calibration) else as.data.table(calibration)
    
    ac <- checkmate::makeAssertCollection()
    
    coln <- names(calibration)
    
    maybeTakeMS2QCol <- function(mcol, pcol)
    {
        if (!pcol %in% coln && mcol %in% coln)
        {
            printf("NOTE: using MS2Quant column '%s' for '%s'\n", mcol, pcol)
            setnames(calibration, mcol, pcol)
        }
    }

    maybeTakeMS2QCol("identifier", "name")
    maybeTakeMS2QCol("retention_time", "rt")
    maybeTakeMS2QCol("area", "intensity")
    maybeTakeMS2QCol("conc_M", "conc")

    for (col in c("name", "SMILES"))
        assertListVal(calibration, col, checkmate::assertCharacter, any.missing = FALSE, add = add)
    for (col in c("rt", "intensity", "conc"))
    {
        assertListVal(calibration, col, checkmate::assertNumeric, mustExist = TRUE,
                      lower = if (col != "rt") 0 else -Inf, finite = TRUE, add = add)
    }
    
    checkmate::reportAssertions(ac)
    
    calibration[, conc := mapply(convertConc, conc, babelConvert(SMILES, "smi", "MW", mustWork = TRUE),
                                 MoreArgs = list(unitFrom = concUnit, unitTo = "M"))]
    
    return(calibration[])
}

checkQuantEluent <- function(x, fGroups)
{
    ret <- checkmate::checkDataFrame(x, any.missing = FALSE, ncols = 2, types = "numeric")
    if (isTRUE(ret))
        ret <- checkHasNames(x, c("time", "B"))
    if (isTRUE(ret))
        ret <- checkmate::checkNumeric(x$time, lower = 0, finite = TRUE)
    if (isTRUE(ret))
        ret <- checkmate::checkNumeric(x$B, lower = 0, upper = 100, finite = TRUE)
    
    if (isTRUE(ret) && length(fGroups) > 0 && max(x$time) < max(groupInfo(fGroups)$ret))
        ret <- paste("The highest retention time in the eluent table is less than the highest feature retention time.",
                     "Make sure retention times are specified in seconds")
        
    return(ret)
}
assertQuantEluent <- checkmate::makeAssertionFunction(checkQuantEluent)

assertConvertMSFilesArgs <- function(formatFrom, formatTo, overWrite, algorithm, add)
{
    # no adding: should fail first
    assertMSConvAlgo(algorithm)
    
    checkmate::assertChoice(formatTo, getMSConversionFormats(algorithm, "output"), add = add)
    checkmate::assertFlag(overWrite, add = add)
    
    validFormatsFrom <- switch(algorithm,
                               pwiz = getMSFileFormats(),
                               openms = getMSFileFormats("centroid"),
                               bruker = "bruker",
                               # NOTE: for im_collapse more specific format checks are done later
                               im_collapse = c(getMSFileFormats("raw"), getMSFileFormats("ims")),
                               timsconvert = "bruker_ims")
    checkmate::assertChoice(formatFrom, validFormatsFrom, add = add)
}

assertPlotEIXArgs <- function(obj, analysis, groupName, showPeakArea, showFGroupRect, title, groupBy, showLegend,
                              annotate, showProgress, xlim, ylim, add)
{
    aapply(checkmate::assertSubset, . ~ analysis + groupName, list(analyses(obj), names(obj)), empty.ok = TRUE,
           fixed = list(add = add))
    aapply(checkmate::assertFlag, . ~ showPeakArea + showFGroupRect + showLegend + showProgress,
           fixed = list(add = add))
    checkmate::assertString(title, null.ok = TRUE, add = add)
    assertAnaInfoBy(groupBy, analysisInfo(obj), TRUE, null.ok = TRUE, add = add)
    annotate <- checkmate::matchArg(annotate, c("none", "ret", "mz", "mob"), several.ok = TRUE, add = add)
    assertXYLim(xlim, ylim, add = add)
}

assertFindPeakParams <- function(x, null.ok = FALSE, .var.name = checkmate::vname(x), add = NULL)
{
    if (null.ok && is.null(x))
        return(invisible(NULL))
    
    assertListVal(x, "algorithm", checkmate::assertChoice, choices = c("xcms3", "envipick", "openms", "piek"),
                  .var.name = .var.name) # NOTE: don't add, must fail for next line
    if (x$algorithm == "openms")
    {
        assertListVal(x, "minPeakWidth", checkmate::assertNumber, finite = TRUE, .var.name = .var.name, add = add)
        assertListVal(x, "backgroundSubtraction", checkmate::assertChoice, choices = c("none", "original", "exact"),
                      .var.name = .var.name, add = add)
        assertListVal(x, "SGolayFrameLength", checkmate::assertCount, positive = TRUE, .var.name = .var.name, add = add)
        assertListVal(x, "SGolayPolyOrder", checkmate::assertCount, positive = TRUE, .var.name = .var.name, add = add)
        assertListVal(x, "useGauss", checkmate::assertFlag, .var.name = .var.name, add = add)
        assertListVal(x, "gaussWidth", checkmate::assertNumber, finite = TRUE, .var.name = .var.name, add = add)
        assertListVal(x, "SN", checkmate::assertNumber, finite = TRUE, .var.name = .var.name, add = add)
        assertListVal(x, "SNWinLen", checkmate::assertNumber, finite = TRUE, .var.name = .var.name, add = add)
        assertListVal(x, "SNBinCount", checkmate::assertCount, positive = TRUE, .var.name = .var.name, add = add)
        assertListVal(x, "method", checkmate::assertChoice, choices = c("legacy", "corrected", "crawdad"),
                      .var.name = .var.name, add = add)
        assertListVal(x, "integrationType", checkmate::assertChoice,
                      choices = c("intensity_sum", "simpson", "trapezoid"), .var.name = .var.name, add = add)
        assertListVal(x, "baselineType", checkmate::assertChoice,
                      choices = c("base_to_base", "vertical_division", "vertical_division_min", "vertical_division_max"),
                      .var.name = .var.name, add = add)
        assertListVal(x, "fitEMG", checkmate::assertFlag, .var.name = .var.name, add = add)
        assertListVal(x, "extraOpts", checkmate::assertList, .var.name = .var.name, add = add)
    }
    else if (x$algorithm == "piek")
    {
        assertListVal(x, "minIntensity", checkmate::assertNumber, finite = TRUE, .var.name = .var.name, add = add)
        assertListVal(x, "SN", checkmate::assertCount, .var.name = .var.name, add = add)
        assertListVal(x, "peakWidth", checkmate::assertNumeric, finite = TRUE, len = 2, any.missing = FALSE,
                      .var.name = .var.name, add = add)
        assertListVal(x, "RTRange", checkmate::assertNumeric, len = 2, any.missing = FALSE, .var.name = .var.name,
                      add = add)
        assertListVal(x, "maxPeaksPerSignal", checkmate::assertCount, .var.name = .var.name, add = add)
    }
    # NOTE: for XCMS/enviPick we just let the package functions throw an error...
    
    assertListVal(x, "forcePeakRange", assertRange, .var.name = .var.name, add = add)
    assertListVal(x, "relMinIntensity", checkmate::assertNumber, .var.name = .var.name, add = add)
    
    invisible(NULL)
}

assertFindMobilitiesArgs <- function(mobPeakParams, chromPeakParams, EIMParams, EICParams, peakRTWindow, fallbackEIC,
                                     calcArea, IMSWindow, CCSParams, parallel, add = NULL)
{
    assertFindPeakParams(mobPeakParams, null.ok = TRUE, add = add)
    assertFindPeakParams(chromPeakParams, null.ok = TRUE, add = add)
    assertEIMParams(EIMParams, add = add)
    assertEICParams(EICParams, add = add)
    aapply(checkmate::assertNumber, . ~ peakRTWindow + IMSWindow, finite = TRUE, fixed = list(add = add))
    aapply(checkmate::assertFlag, . ~ fallbackEIC + parallel, fixed = list(add = add))
    checkmate::assertChoice(calcArea, c("integrate", "sum"), add = add)
    assertCCSParams(CCSParams, null.ok = TRUE, add = add)
    invisible(NULL)
}

assertCCSParams <- function(x, null.ok = FALSE, .var.name = checkmate::vname(x), add = NULL)
{
    if (null.ok && is.null(x))
        return(invisible(NULL))
    
    assertListVal(x, "method", checkmate::assertChoice,
                  choices = c("bruker", "mason-schamp_k", "mason-schamp_1/k", "agilent"),# "waters"),
                  .var.name = .var.name, add = add)
    
    if (x$method == "bruker" && !backendAvailable("opentims"))
        stop("Cannot calculate mobility data with Bruker library: not supported on this platform/OS", call. = FALSE)
    
    if (x$method == "agilent" && is.null(x[["calibrant"]]))
        stop("Please set the calibrant CCS parameter", call. = FALSE)
    
    assertListVal(x, "defaultCharge", checkmate::assertInt, .var.name = .var.name, add = add)
    
    assertListVal(x, "calibrant", function(..., .var.name)
    {
        checkmate::assert(
                  checkmate::checkNull(...),
                  checkmate::checkDirectoryExists(...),
                  checkmate::checkFileExists(..., "r"),
                  checkmate::checkList(..., any.missing = FALSE),
                  .var.name = .var.name
        )
    }, .var.name = .var.name)
    
    if (is.list(x$calibrant))
    {
        vn <- paste0(.var.name, "$calibrant")
        assertListVal(x$calibrant, "massGas", checkmate::assertNumber, finite = TRUE, mustExist = FALSE,
                      .var.name = vn, add = add)
        assertListVal(x$calibrant, "TFix", checkmate::assertNumber, finite = TRUE, .var.name = vn, add = add)
        assertListVal(x$calibrant, "beta", checkmate::assertNumber, finite = TRUE, .var.name = vn, add = add)
    }
}

assertMobilityConversionArgs <- function(mobility, mz, CCSParams, charge, add = NULL)
{
    aapply(checkmate::assertNumeric, . ~ mobility + mz, finite = TRUE, lower = 0, fixed = list(add = add))
    checkmate::assertIntegerish(charge, null.ok = TRUE, add = add)
    assertCCSParams(CCSParams, add = add)
    
    if (length(mobility) != length(mz))
        stop("The length of the input mobility and mz data should be equal", call. = FALSE)
}

assertIMSRangeParams <- function(x, null.ok = FALSE, .var.name = checkmate::vname(x), add = NULL)
{
    if (null.ok && is.null(x))
        return(invisible(NULL))
    
    assertListVal(x, "param", checkmate::assertChoice, choices = c("mobility", "CCS"), .var.name = .var.name, add = add)
    assertListVal(x, "lower", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name, add = add)
    assertListVal(x, "upper", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name, add = add)
    assertListVal(x, "mzRelative", checkmate::assertFlag, .var.name = .var.name, add = add)
    
}

assertIMSMatchParams <- function(x, null.ok = FALSE, .var.name = checkmate::vname(x), add = NULL)
{
    if (null.ok && is.null(x))
        return(invisible(NULL))
    
    assertListVal(x, "param", checkmate::assertChoice, choices = c("mobility", "CCS"), .var.name = .var.name, add = add)
    assertListVal(x, "relative", checkmate::assertFlag, .var.name = .var.name, add = add)
    assertListVal(x, "window", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name, add = add)
    assertListVal(x, "minMatches", checkmate::assertCount, positive = FALSE, .var.name = .var.name, add = add)
}

assertPiekGenEICParams <- function(x, .var.name = checkmate::vname(x), add = NULL)
{
    checkmate::assertList(x, .var.name = .var.name)
    assertListVal(x, "methodMZ", checkmate::assertChoice, choices = c("bins", "suspects", "ms2"), .var.name = .var.name)
    assertListVal(x, "methodIMS", checkmate::assertChoice, choices = c("bins", "suspects", "ms2"), null.ok = TRUE,
                  .var.name = .var.name)
    assertListVal(x, "retRange", assertRange, null.ok = TRUE, .var.name = .var.name, add = add)
    assertListVal(x, "minEICIntensity", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name, add = add)
    assertListVal(x, "minEICAdjTime", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name, add = add)
    assertListVal(x, "minEICAdjPoints", checkmate::assertCount, .var.name = .var.name, add = add)
    assertListVal(x, "minEICAdjIntensity", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name, add = add)
    assertListVal(x, "topMostEIC", checkmate::assertCount, positive = FALSE, .var.name = .var.name, add = add)
    assertListVal(x, "topMostEICPre", checkmate::assertCount, positive = FALSE, .var.name = .var.name, add = add)
    
    if (x$methodMZ == "bins")
    {
        assertListVal(x, "mzRange", assertRange, .var.name = .var.name, add = add)
        assertListVal(x, "mzStep", checkmate::assertNumber, lower = 0.000001, finite = TRUE,
                      .var.name = .var.name, add = add)
    }
    else if (x$methodMZ == "suspects")
    {
        assertListVal(x, "rtWindow", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name, add = add)
        assertListVal(x, "mzWindow", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name, add = add)
        
        # UNDONE these to separate param and also use elsewhere?
        assertListVal(x, "skipInvalid", checkmate::assertFlag, .var.name = .var.name, add = add)
        assertListVal(x, "prefCalcChemProps", checkmate::assertFlag, .var.name = .var.name, add = add)
        assertListVal(x, "neutralChemProps", checkmate::assertFlag, .var.name = .var.name, add = add)
    }
    else if (x$methodMZ == "ms2")
    {
        assertListVal(x, "rtWindow", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name, add = add)
        assertListVal(x, "mzWindow", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name, add = add)
        assertListVal(x, "minTIC", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name, add = add)
        assertListVal(x, "clusterMethod", assertClusterMethod, .var.name = .var.name, add = add)
    }

    if (!is.null(x[["methodIMS"]]))
    {
        if (x$methodIMS != "bins" && x$methodIMS != x$methodMZ)
            stop("'methodIMS' should be set to 'bins' or match 'methodMZ'.", call. = FALSE)
        
        if (x$methodIMS == "bins")
        {
            assertListVal(x, "mobRange", assertRange, .var.name = .var.name, add = add)
            assertListVal(x, "mobStep", checkmate::assertNumber, lower = 0.000001, finite = TRUE,
                          .var.name = .var.name, add = add)
        }
        else if (x$methodIMS == "suspects")
        {
            assertListVal(x, "IMSWindow", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name,
                          add = add)
        }
        else if (x$methodIMS == "ms2")
        {
            assertListVal(x, "IMSWindow", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name,
                          add = add)
        }
    }
}

# from https://github.com/mllg/checkmate/issues/115
aapply = function(fun, formula, ..., fixed = list())
{
    fun = match.fun(fun)
    terms = terms(formula)
    vnames = attr(terms, "term.labels")
    ee = attr(terms, ".Environment")

    dots = list(...)
    dots$.var.name = vnames
    dots$x = unname(mget(vnames, envir = ee))
    .mapply(fun, dots, MoreArgs = fixed)

    invisible(NULL)
}
