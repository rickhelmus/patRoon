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

assertAnalysisInfo <- function(x, allowedFormats = NULL, verifyCentroided = FALSE, null.ok = FALSE,
                               .var.name = checkmate::vname(x), add = NULL)
{
    if (is.null(x) && null.ok)
        return(TRUE)

    if (!is.null(add))
        mc <- length(add$getMessages())

    checkmate::assertDataFrame(x, min.rows = 1, .var.name = .var.name, add = add)
    assertHasNames(x, c("path", "analysis", "group", "blank"), .var.name = .var.name, add = add)

    assertListVal(x, "path", checkmate::assertCharacter, any.missing = FALSE, .var.name = .var.name, add = add)
    assertListVal(x, "analysis", checkmate::assertCharacter, any.missing = FALSE, .var.name = .var.name, add = add)
    assertListVal(x, "group", checkmate::assertCharacter, any.missing = FALSE, .var.name = .var.name, add = add)
    assertListVal(x, "blank", checkmate::assertCharacter, any.missing = TRUE, .var.name = .var.name, add = add)

    checkmate::assert(
        checkmate::checkNull(x[["conc"]]),
        checkmate::checkCharacter(x[["conc"]]),
        checkmate::checkNumeric(x[["conc"]]),
        .var.name = sprintf("%s[[\"conc\"]]", .var.name)
    )

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
        checkmate::assertDirectoryExists(x$path, .var.name = .var.name, add = add)

        # UNDONE: more extensions? (e.g. mzData)
        if (is.null(allowedFormats))
            allowedFormats <- MSFileFormats()
        
        # stops if files are missing
        getMSFilePaths(x$analysis, x$path, allowedFormats, mustExist = TRUE)

        if (verifyCentroided)
            verifyDataCentroided(x)
        
        checkmate::assertVector(x$analysis, unique = TRUE, .var.name = paste0(.var.name, "$analysis"), add = add)
    }

    invisible(NULL)
}

assertAndPrepareAnaInfo <- function(x, ..., add = NULL)
{
    if (!is.null(add))
        mc <- length(add$getMessages())

    assertAnalysisInfo(x, ..., add = add)

    if ((is.null(add) || length(add$getMessages()) == mc) && !is.null(x))
    {
        x <- makeDT(x) # convert to DT or make a unique copy
        x <- unFactorDT(x)

        if (is.null(x[["blank"]]) && !is.null(x[["ref"]]))
        {
            warning("The usage of a 'ref' column in the analysis information is deprecated. Please re-name this column to 'blank'.")
            setnames(x, "ref", "blank")
        }
        
        if (!is.null(x[["conc"]]))
            x[, conc := as.numeric(conc)]
        if (!is.null(x[["norm_conc"]]))
            x[, norm_conc := as.numeric(norm_conc)]
        x[is.na(blank), blank := ""]
    }
    
    return(x)
}

assertAndPrepareAnaInfoAverage <- function(x, anaInfo, .var.name = checkmate::vname(x), add = NULL)
{
    checkmate::assert(
        checkmate::checkFlag(x),
        checkmate::checkChoice(x, names(anaInfo)),
        .var.name = .var.name, add = add
    )
    if (isTRUE(x))
        x <- "group"
    else if (isFALSE(x))
        x <- "analysis"
    return(x)
}

assertSuspectList <- function(x, needsAdduct, skipInvalid, .var.name = checkmate::vname(x), add = NULL)
{
    mzCols <- c("mz", "neutralMass", "SMILES", "InChI", "formula")
    allCols <- c("name", "adduct", "rt", mzCols)
    
    # this seems necessary for proper naming in subsequent assertions (why??)
    .var.name <- force(.var.name)
    
    # subset with relevant columns: avoid checking others in subsequent assertDataFrame call
    if (checkmate::testDataFrame(x))
    {
        if (is.data.table(x))
            x <- x[, intersect(names(x), allCols), with = FALSE]
        else
            x <- x[, intersect(names(x), allCols), drop = FALSE]
    }
    
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

assertFGAsDataTableArgs <- function(fGroups, areas, features, qualities, regression, regressionBy, averageFunc,
                                    normalized, FCParams, concAggrParams, toxAggrParams, normConcToTox, collapseSuspects,
                                    onlyHits)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertFlag, . ~ areas + features + regression + normalized + normConcToTox,
           fixed = list(add = ac))
    checkmate::assertString(regressionBy, na.ok = FALSE, min.chars = 1, null.ok = TRUE, add = ac)
    checkmate::assertFunction(averageFunc, add = ac)
    assertFCParams(FCParams, fGroups, null.ok = TRUE, add = ac)
    aapply(assertPredAggrParams, . ~ concAggrParams + toxAggrParams, null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertString(collapseSuspects, null.ok = TRUE, add = ac)
    checkmate::assertFlag(onlyHits, add = ac)
    checkmate::reportAssertions(ac)
    
    checkmate::assert(checkmate::checkFALSE(qualities),
                      checkmate::checkChoice(qualities, c("quality", "score", "both")),
                      .var.name = "qualities")
}

assertNormalizationMethod <- function(x, withNone = TRUE, .var.name = checkmate::vname(x), add = NULL)
{
    ch <- c("max", "minmax")
    if (withNone)
        ch <- c(ch, "none")
    checkmate::assertChoice(x, ch, .var.name = .var.name, add = add)
}

assertEICParams <- function(x, .var.name = checkmate::vname(x), add = NULL)
{
    checkmate::assertList(x, names = "unique", .var.name = .var.name) # no add: should fail
    
    assertListVal(x, "rtWindow", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name, add = add)
    assertListVal(x, "mzExpWindow", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name, add = add)
    assertListVal(x, "topMost", checkmate::assertCount, positive = TRUE, null.ok = TRUE, .var.name = .var.name,
                  add = add)
    assertListVal(x, "topMostByRGroup", checkmate::assertFlag, .var.name = .var.name, add = add)
    assertListVal(x, "onlyPresent", checkmate::assertFlag, .var.name = .var.name, add = add)
    assertListVal(x, "setsAdductPos", checkAndToAdduct, .var.name = .var.name)
    assertListVal(x, "setsAdductNeg", checkAndToAdduct, .var.name = .var.name)
    
    if (!is.null(x[["topMost"]]) && !isTRUE(x$onlyPresent))
        stop("onlyPresent must be TRUE if topMost is set", call. = FALSE)
    
    invisible(NULL)
}

assertFCParams <- function(x, fGroups, null.ok = FALSE, .var.name = checkmate::vname(x), add = NULL)
{
    if (null.ok && is.null(x))
        return(NULL)
    
    checkmate::assertList(x, names = "unique", .var.name = .var.name) # no add: should fail
    
    assertListVal(x, "rGroups", checkmate::assertCharacter, any.missing = FALSE, len = 2, .var.name = .var.name,
                  add = add)
    assertListVal(x, "rGroups", checkmate::assertSubset, choices = replicateGroups(fGroups), .var.name = .var.name,
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
    assertListVal(x, "topMost", checkmate::assertCount, positive = TRUE, .var.name = .var.name, add = add)
    assertListVal(x, "minIntensityPre", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name,
                  add = add)
    assertListVal(x, "minIntensityPost", checkmate::assertNumber, lower = 0, finite = TRUE, .var.name = .var.name,
                  add = add)
    assertListVal(x, "avgFun", checkmate::assertFunction, .var.name = .var.name, add = add)
    assertListVal(x, "method", checkmate::assertChoice, choices = c("distance", "hclust"), .var.name = .var.name,
                  add = add)
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
    if (isTRUE(ret) && any(grepl(",", x, fixed = TRUE)))
        ret <- "Set labels cannot contain commas"
    if (isTRUE(ret) && any(grepl("-", x, fixed = TRUE)))
        ret <- "Set labels cannot contain minus signs (-)"
    if (isTRUE(ret) && any(grepl("genform|sirius|bruker|metfrag", x)))
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
    
    if (isTRUE(ret) && length(fGroups) > 0 && max(x$time) < max(groupInfo(fGroups)$rts))
        ret <- paste("The highest retention time in the eluent table is less than the highest feature retention time.",
                     "Make sure retention times are specified in seconds")
        
    return(ret)
}
assertQuantEluent <- checkmate::makeAssertionFunction(checkQuantEluent)

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
