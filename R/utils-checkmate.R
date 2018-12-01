checkHasNames <- function(x, n, type = "unique") checkmate::checkNames(names(x), must.include = n, type = type)
assertHasNames <- checkmate::makeAssertionFunction(checkHasNames)

checkRange <- function(x, null.ok = FALSE)
{
    ret <- checkmate::checkNumeric(x, any.missing = FALSE, finite = TRUE, len = 2, null.ok = null.ok)
    if (isTRUE(ret) && !is.null(x) && x[2] != -1 && x[1] > x[2])
        ret <- paste0("lower range (", x[1], ") higher than upper (", x[2], ")")
    return(ret)
}
assertRange <- checkmate::makeAssertionFunction(checkRange)

checkChoiceSilent <- function(x, ch)
{
    ret <- checkmate::checkString(x, min.chars = 1)
    if (isTRUE(ret) && !x %in% ch)
        ret <- paste("Must be element of", getStrListWithMax(ch, 6, ", "))
    return(ret)
}
assertChoiceSilent <- checkmate::makeAssertionFunction(checkChoiceSilent)

assertAnalysisInfo <- function(x, allowedFormats = NULL, .var.name = checkmate::vname(x), add = NULL)
{
    if (!is.null(add))
        mc <- length(add$getMessages())
    
    checkmate::assertDataFrame(x, types = "character", any.missing = FALSE, min.rows = 1, .var.name = .var.name, add = add)
    assertHasNames(x, c("path", "analysis", "group", "ref"), .var.name = .var.name, add = add)
    
    # only continue if previous assertions didn't fail: x needs to be used as list which otherwise gives error
    # NOTE: this is only applicable if add != NULL, otherwise previous assertions will throw errors
    if (is.null(add) || length(add$getMessages()) == mc)
    {
        checkmate::assertDirectoryExists(x$path, .var.name = .var.name, add = add)
        
        # UNDONE: more extensions? (e.g. mzData)
        if (is.null(allowedFormats))
            allowedFormats <- c("mzML", "mzXML", "d")
        
        res <- FALSE
        for (f in allowedFormats)
        {
            p <- file.path(x$path, paste0(x$analysis, ".", f))
            if (f == "d")
                res <- checkmate::checkDirectoryExists(p)
            else
                res <- checkmate::checkFileExists(p)
            
            if (isTRUE(res))
                break
        }
        
        allowedFormats <- sub("d", "bruker", allowedFormats) # user friendly name
        if (!isTRUE(res))
            checkmate::makeAssertion(x, sprintf("No analyses found with correct data format (valid: %s)",
                                                paste0(allowedFormats, collapse = ", ")),
                                     var.name = .var.name, collection = add)
    }
    
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

assertMultiProcArgs <- function(x, maxProcAmount, .var.name = checkmate::vname(x), add = NULL)
{
    if (!is.null(x))
        assertCanCreateDir(x, .var.name = .var.name, add = add)
    # HACK: fix var name
    checkmate::assertCount(maxProcAmount, positive = TRUE, .var.name = "maxProcAmount", add = add)
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
    ret <- checkmate::checkIntegerish(x, lower = 0)
    if (!isTRUE(ret))
        ret <- checkmate::checkCharacter(x)
    if (!isTRUE(ret))
        ret <- checkmate::checkLogical(x)
    if (!isTRUE(ret))
        ret <- "Should be valid numeric, character or logical"
    return(ret)
}
assertSubsetArg <- checkmate::makeAssertionFunction(checkSubsetArg)

# used for "[[" methods
checkExtractArg <- function(x)
{
    ret <- checkmate::checkInt(x, lower = 0)
    if (!isTRUE(ret))
        ret <- checkmate::checkString(x)
    if (!isTRUE(ret))
        ret <- "Should be valid numeric or character skalar"
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
    assertVal(checkmate::assertNumber, "minIntensity", lower = 0, finite = TRUE)
    assertVal(checkmate::assertFunction, "avgFun")
    assertVal(checkmate::assertChoice, "method", choices = c("distance", "hclust"))
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
