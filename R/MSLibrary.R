#' @include main.R
#' @include workflow-step.R
NULL

#' @export
MSLibrary <- setClass("MSLibrary", slots = c(records = "data.table", spectra = "list"),
                      contains = "workflowStep")

#' @export
setMethod("records", "MSLibrary", function(obj) obj@records)

#' @export
setMethod("spectra", "MSLibrary", function(obj) obj@spectra)

#' @export
setMethod("length", "MSLibrary", function(x) nrow(records(x)))

#' @export
setMethod("names", "MSLibrary", function(x) names(spectra(x)))

#' @export
setMethod("show", "MSLibrary", function(object)
{
    callNextMethod()
    
    # UNDONE: more?
    printf("Total records: %d\n", length(object))
})

#' @export
setMethod("[", c("MSLibrary", "ANY", "missing", "missing"), function(x, i, ...)
{
    if (!missing(i))
    {
        i <- assertSubsetArgAndToChr(i, names(x))
        x@records <- x@records[DB_ID %in% i]
        x@spectra <- x@spectra[i]
    }
    
    return(x)
})

#' @export
setMethod("[[", c("MSLibrary", "ANY", "missing"), function(x, i, j)
{
    assertExtractArg(i)
    return(x@spectra[[i]])
})

#' @export
setMethod("$", "MSLibrary", function(x, name)
{
    eval(substitute(x@spectra$NAME_ARG, list(NAME_ARG = name)))
})

#' @export
setMethod("as.data.table", "MSLibrary", function(x)
{
    allSpecs <- rbindlist(spectra(x), idcol = "DB_ID")
    return(merge(records(x), allSpecs, by = "DB_ID"))
})

#' @export
setMethod("convertToSuspects", "MSLibrary", function(obj)
{
    if (length(obj) == 0)
        stop("Cannot create suspect list: no data", call. = FALSE)
    
    ret <- copy(records(obj))
    
    mapCols <- c(Name = "name",
                 SMILES = "SMILES",
                 InChI = "InChI",
                 InChIKey = "InChIKey",
                 Formula = "formula",
                 Precursor_Type = "adduct",
                 ExactMass = "neutralMass")
    mapCols <- mapCols[names(mapCols) %in% names(ret)]
    setnames(ret, names(mapCols), mapCols)
    ret <- unique(ret, by = "InChIKey")
    ret <- ret[, mapCols, with = FALSE]
    ret <- prepareSuspectList(ret, NULL, FALSE, FALSE)
    
    return(ret)
})


loadMSPLibrary <- function(file)
{
    checkmate::assertFileExists(file, "r")
    
    lib <- readMSP(file)
    lib$records <- as.data.table(lib$records)
    lib$spectra <- lapply(lib$spectra, as.data.table)
    
    return(MSLibrary(records = lib$records, spectra = lib$spectra, algorithm = "msp"))
}
