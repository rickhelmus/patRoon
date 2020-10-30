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

assertAnalysisInfo <- function(x, allowedFormats = NULL, null.ok = FALSE, .var.name = checkmate::vname(x), add = NULL)
{
    if (is.null(x) && null.ok)
        return(TRUE)

    if (!is.null(add))
        mc <- length(add$getMessages())

    checkmate::assertDataFrame(x, min.rows = 1, .var.name = .var.name, add = add)
    assertHasNames(x, c("path", "analysis", "group", "blank"), .var.name = .var.name, add = add)

    assertCharField <- function(f) checkmate::assertCharacter(x[[f]], .var.name = sprintf("%s[\"%s\"]", .var.name, f), add = add)
    assertCharField("path")
    assertCharField("analysis")
    assertCharField("group")
    assertCharField("blank")

    checkmate::assert(
        checkmate::checkNull(x[["conc"]]),
        checkmate::checkCharacter(x[["conc"]]),
        checkmate::checkNumeric(x[["conc"]]),
        .var.name = sprintf("%s[[\"conc\"", .var.name)
    )

    # only continue if previous assertions didn't fail: x needs to be used as list which otherwise gives error
    # NOTE: this is only applicable if add != NULL, otherwise previous assertions will throw errors
    if (is.null(add) || length(add$getMessages()) == mc)
    {
        checkmate::assertDirectoryExists(x$path, .var.name = .var.name, add = add)

        # UNDONE: more extensions? (e.g. mzData)
        if (is.null(allowedFormats))
            allowedFormats <- MSFileFormats()

        exts <- unique(unlist(MSFileExtensions()[allowedFormats]))

        existFiles <- mapply(x$path, x$analysis, FUN = function(path, ana)
        {
            for (f in allowedFormats)
            {
                exts <- MSFileExtensions()[[f]]
                for (e in exts)
                {
                    p <- file.path(path, paste0(ana, ".", e))
                    if (file.exists(p) && file.info(p, extra_cols = FALSE)$isdir == MSFileFormatIsDir(f, e))
                        return(TRUE)
                }
            }
            message(sprintf("Analysis does not exist: %s (in %s)", ana, path))
            return(FALSE)
        })
        
        if (any(!existFiles))
            checkmate::makeAssertion(x, sprintf("No analyses found with correct data format (valid: %s)",
                                                paste0(allowedFormats, collapse = ", ")),
                                     var.name = .var.name, collection = add)
    }

    invisible(NULL)
}

assertAndPrepareAnaInfo <- function(x, ..., add = NULL)
{
    if (!is.null(x))
        x <- unFactorDF(x)

    if (!is.null(add))
        mc <- length(add$getMessages())

    if (!is.null(x) && checkmate::checkDataFrame(x) && is.null(x[["blank"]]) && !is.null(x[["ref"]]))
    {
        warning("The usage of a 'ref' column in the analysis information is deprecated. Please re-name this column to 'blank'.")
        setnames(x, "ref", "blank")
    }

    assertAnalysisInfo(x, ..., add = add)

    if ((is.null(add) || length(add$getMessages()) == mc) &&
        (!is.null(x) && !is.null(x[["conc"]])))
        x[["conc"]] <- as.numeric(x[["conc"]])

    return(x)
}

assertSuspectList <- function(x, adduct, skipInvalid, .var.name = checkmate::vname(x), add = NULL)
{
    mzCols <- c("mz", "neutralMass", "SMILES", "InChI", "formula")
    allCols <- c("name", "adduct", "rt", mzCols)
    
    # this seems necessary for proper naming in subsequent assertions (why??)
    .var.name <- force(.var.name)
    
    # subset with relevant columns: avoid checking others in subsequent assetDataFrame call
    if (checkmate::testDataFrame(x))
    {
        if (is.data.table(x))
            x <- x[, intersect(names(x), allCols), with = FALSE]
        else
            x <- x[, intersect(names(x), allCols)]
    }
    
    checkmate::assertDataFrame(x, any.missing = skipInvalid, min.rows = 1, .var.name = .var.name, add = add)
    assertHasNames(x, "name", .var.name = .var.name, add = add)
    

    checkmate::assertNames(intersect(names(x), mzCols), subset.of = mzCols,
                           .var.name = paste0("names(", .var.name, ")"), add = add)

    assertCharField <- function(f, null.ok = TRUE)
    {
        checkmate::assert(
            checkmate::checkFactor(x[[f]], empty.levels.ok = skipInvalid,
                                   any.missing = skipInvalid, null.ok = null.ok),
            checkmate::checkCharacter(x[[f]], min.chars = if (skipInvalid) 0 else 1,
                                      null.ok = null.ok),
            .var.name = sprintf("%s[\"%s\"]", .var.name, f)    
        )
        
    }
    assertCharField("SMILES"); assertCharField("InChI"); assertCharField("formula");
    assertCharField("adduct", null.ok = !is.null(adduct) || !is.null(x[["mz"]]))

    assertNumField <- function(f) checkmate::assertNumeric(x[[f]], .var.name = sprintf("%s[\"%s\"]", .var.name, f),
                                                           lower = 0, finite = TRUE,
                                                           null.ok = TRUE, add = add)
    assertNumField("mz"); assertNumField("neutralMass"); assertNumField("rt")

    invisible(NULL)
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
    checkmate::assertNumber(absMinAbundance, .var.name = "absMinAbundance", null.ok = TRUE, add = ac)
    checkmate::assertNumber(relMinAbundance, .var.name = "relMinAbundance", null.ok = TRUE, add = ac)
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

assertNormalizationMethod <- function(x, .var.name = checkmate::vname(x), add = NULL)
{
    checkmate::assertChoice(x, c("none", "max", "minmax"), .var.name = .var.name, add = add)
}

assertAvgPListParams <- function(x, .var.name = checkmate::vname(x), add = NULL)
{
    checkmate::assertList(x, names = "unique", .var.name = .var.name) # no add: should fail

    assertVal <- function(f, v, ...) f(x[[v]], ..., .var.name = paste0(.var.name, "$", v), add = add)

    assertVal(checkmate::assertNumber, "clusterMzWindow", lower = 0, finite = TRUE)
    assertVal(checkmate::assertCount, "topMost", positive = TRUE)
    assertVal(checkmate::assertNumber, "minIntensityPre", lower = 0, finite = TRUE)
    assertVal(checkmate::assertNumber, "minIntensityPost", lower = 0, finite = TRUE)
    assertVal(checkmate::assertFunction, "avgFun")
    assertVal(checkmate::assertChoice, "method", choices = c("distance", "hclust"))
    assertVal(checkmate::assertFlag, "retainPrecursorMSMS")
}

assertPListIsolatePrecParams <- function(x, .var.name = checkmate::vname(x), add = NULL)
{
    if (is.null(x))
        return(NULL)

    checkmate::assertList(x, names = "unique", .var.name = .var.name) # no add: should fail

    assertVal <- function(f, v, ...) f(x[[v]], ..., .var.name = paste0(.var.name, "$", v), add = add)

    assertVal(checkmate::assertCount, "maxIsotopes")
    assertVal(checkmate::assertNumeric, "mzDefectRange", any.missing = FALSE, len = 2, finite = TRUE)
    assertVal(checkmate::assertNumeric, "intRange", any.missing = FALSE, len = 2, finite = TRUE)
    assertVal(checkmate::assertCount, "z", positive = TRUE)
    assertVal(checkmate::assertCount, "maxGap", positive = TRUE)
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
