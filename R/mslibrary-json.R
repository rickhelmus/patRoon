#' @include mslibrary.R
NULL

loadMSLibraryMoNAJSON <- function(file, potAdducts = TRUE, potAdductsLib = TRUE, absMzDev = 0.002, calcSPLASH = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFileExists(file, "r", add = ac)
    checkmate::assertNumber(absMzDev, lower = 0, finite = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ potAdductsLib + calcSPLASH, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    hash <- makeHash(makeFileHash(file), potAdducts, potAdductsLib, absMzDev, calcSPLASH)
    cd <- loadCacheData("MSLibraryJSON", hash)
    if (!is.null(cd))
        return(cd)
    
    lib <- readMoNAJSON(normalizePath(file))
    lib <- sanitizeMSLibrary(lib, potAdducts, potAdductsLib, absMzDev, calcSPLASH)
    
    ret <- MSLibrary(records = lib$records[], spectra = lib$spectra, algorithm = "json")
    
    saveCacheData("MSLibraryJSON", ret, hash)
    
    return(ret)
}
