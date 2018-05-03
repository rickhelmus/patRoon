checkHasNames <- function(x, n, type = "unique") checkmate::checkNames(names(x), must.include = n, type = type)
assertHasNames <- checkmate::makeAssertionFunction(checkHasNames)

assertAnalysisInfo <- function(x, requiredFormat = NULL, allowedFormats = NULL, .var.name = checkmate::vname(x), add = NULL)
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
        if (!is.null(requiredFormat))
            checkmate::assertFileExists(file.path(x$path, paste0(x$analysis, ".", requiredFormat)),
                                        .var.name = .var.name, add = add)
        else
        {
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
            
            if (!isTRUE(res))
                checkmate::makeAssertion(x, sprintf("No valid data format found (valid: %s)",
                                                    paste0(allowedFormats, collapse = ", ")),
                                         var.name = .var.name, collection = add)
        }
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

assertMultiProcArgs <- function(x, maxProcAmount, .var.name = checkmate::vname(x), add = NULL)
{
    if (!is.null(x))
        assertCanCreateDir(x, .var.name = .var.name, add = add)
    # HACK: fix var name
    checkmate::assertCount(maxProcAmount, positive = TRUE, .var.name = "maxProcAmount", add = add)
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
