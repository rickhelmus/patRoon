#' @include mslibrary.R
NULL

loadMSLibraryMSP <- function(file, parseComments = TRUE, potAdducts = NULL, potAdductsLib = TRUE, absMzDev = 0.002,
                             calcSPLASH = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFileExists(file, "r", add = ac)
    aapply(checkmate::assertFlag, . ~ parseComments + potAdductsLib + calcSPLASH, fixed = list(add = ac))
    checkmate::assert(checkmate::checkNull(potAdducts),
                      checkmate::checkFALSE(potAdducts),
                      checkmate::checkCharacter(potAdducts, any.missing = FALSE, min.chars = 1),
                      checkmate::checkList(potAdducts, types = c("adduct", "character"), any.missing = FALSE),
                      .var.name = "potAdducts")
    checkmate::assertNumber(absMzDev, lower = 0, finite = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    hash <- makeHash(makeFileHash(file), parseComments, potAdducts, potAdductsLib, absMzDev, calcSPLASH)
    cd <- loadCacheData("MSLibraryMSP", hash)
    if (!is.null(cd))
        return(cd)
    
    lib <- readMSP(normalizePath(file), parseComments)
    lib <- sanitizeMSLibrary(lib, potAdducts, potAdductsLib, absMzDev, calcSPLASH)
    
    ret <- MSLibrary(records = lib$records[], spectra = lib$spectra, algorithm = "msp")
    
    saveCacheData("MSLibraryMSP", ret, hash)
    
    return(ret)
}
