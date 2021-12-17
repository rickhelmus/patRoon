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
                 Precursor_type = "adduct",
                 ExactMass = "neutralMass")
    mapCols <- mapCols[names(mapCols) %in% names(ret)]
    setnames(ret, names(mapCols), mapCols)
    ret <- unique(ret, by = "InChIKey")
    ret <- ret[, mapCols, with = FALSE]
    ret <- prepareSuspectList(ret, NULL, FALSE, FALSE)
    
    return(ret)
})

setMethod("export", "MSLibrary", function(obj, type, out)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertChoice(type, "msp", add = ac)
    checkmate::assertPathForOutput(out, overwrite = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    # convert records to character matrix to simplfy Rcpp processing
    recs <- copy(records(obj))
    recs[, (names(recs)) := lapply(.SD, as.character)]
    recs <- as.matrix(recs)
    writeMSPLibrary(recs, spectra(obj), normalizePath(out))
})



loadMSPLibrary <- function(file, parseComments = TRUE, potAdducts = NULL, absMzDev = 0.002)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFileExists(file, "r", add = ac)
    aapply(checkmate::assertFlag, . ~ parseComments, fixed = list(add = ac))
    checkmate::assert(checkmate::checkNull(potAdducts),
                      checkmate::checkFALSE(potAdducts),
                      checkmate::checkCharacter(potAdducts, any.missing = FALSE, min.chars = 1),
                      checkmate::checkList(potAdducts, types = c("adduct", "character"), any.missing = FALSE),
                      .var.name = "potAdducts")
    checkmate::assertNumber(absMzDev, lower = 0, finite = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    lib <- readMSP(normalizePath(file), parseComments)
    
    printf("Converting to tables\n")
    lib$records <- as.data.table(lib$records)
    # UNDONE: not for now, takes a lot of time for large amounts
    # lib$spectra <- withProg(length(lib$spectra), FALSE, lapply(lib$spectra, function(s)
    # {
    #     ret <- as.data.table(s)
    #     doProgress()
    #     return(ret)
    # }))
    
    # C++ code sets "NA" as string, convert to NA. Similarly, library may have 'n/a' markers...
    for (j in seq_along(lib$records))
        set(lib$records, which(lib$records[[j]] %chin% c("NA", "n/a", "N/A")), j, NA_character_)
    
    # Ensure case of column names used by patRoon are consistent
    chCols <- c("Name", "SMILES", "InChI", "InChIKey", "Formula", "Precursor_type", "Ion_mode", "Splash")
    numCols <- c("ExactMass", "PrecursorMZ", "MW")
    allCols <- c(chCols, numCols)
    change <- match(tolower(allCols), tolower(names(lib$records)), nomatch = integer())
    setnames(lib$records, change, allCols)
    
    if (!all(c("Name", "DB_ID") %in% names(lib$records)))
        stop("MSP file misses mandatory Name and/or DB# data.")
    
    # make sure columns at least exist, which makes future checks easier
    for (col in allCols)
    {
        if (is.null(lib$records[[col]]))
            lib$records[, (col) := if (col %in% chCols) NA_character_ else NA_real_]
    }

    # cleanup character data: trim whitespace and convert empties to NA
    lib$records[, (chCols) := lapply(.SD, function(x)
    {
        x <- trimws(x)
        return(fifelse(nzchar(x), x, NA_character_))
    }), .SDcols = chCols]

    # ensure mass data is numeric, but ignore conversion warnings
    suppressWarnings(lib$records[, (numCols) := lapply(.SD, as.numeric), .SDcols = numCols])
    
    # UNDONE: this is quite slow, skip for now?
    # printf("Sanitizing SMILES/InChI values... ")
    # lib$records[!is.na(SMILES), SMILES := babelConvert(SMILES, "smi", "smi", FALSE)]
    # lib$records[!is.na(InChI), InChI := babelConvert(InChI, "inchi", "inchi", FALSE)]
    # printf("Done!\n")
    
    lib$records <- convertChemDataIfNeeded(lib$records, destFormat = "smi", destCol = "SMILES",
                                           fromFormats = "inchi", fromCols = "InChI")
    lib$records <- convertChemDataIfNeeded(lib$records, destFormat = "inchi", destCol = "InChI",
                                           fromFormats = "smi", fromCols = "SMILES")
    lib$records <- convertChemDataIfNeeded(lib$records, destFormat = "inchikey", destCol = "InChIKey",
                                           fromFormats = c("smi", "inchi"), fromCols = c("SMILES", "InChI"))
    lib$records <- convertChemDataIfNeeded(lib$records, destFormat = "formula", destCol = "Formula",
                                           fromFormats = c("smi", "inchi"), fromCols = c("SMILES", "InChI"))
    
    # normalize polarity: ensure uppercase, sometimes shortened as P/N
    lib$records[, Ion_mode := toupper("POSITIVE")]
    lib$records[Ion_mode == "P", Ion_mode := "POSITIVE"]
    lib$records[Ion_mode == "N", Ion_mode := "NEGATIVE"]
    
    # cleanup adducts: it's quite a mess in some DBs, for now convert a few common 'flavors'
    adductMapping <- c("^M\\+H$" = "[M+H]+",
                       "^\\(M\\+H\\)\\+$" = "[M+H]+",
                       "^M\\-H$" = "[M-H]-",
                       "^M\\+Na$" = "[M+Na]+",
                       "^M\\+K$" = "[M+K]+",
                       "^M\\+NH4$" = "[M+HN4]+",
                       "^M\\+Cl$" = "[M+Cl]-", # UNDONE: MoNA records say POSITIVE...
                       "^M\\+Na\\-2H$" = "[M+Na-2H]-",
                       "^M\\+2Na$" = "[M+2Na]2+",
                       "^M\\+$" = "[M]+",
                       "^M\\-$" = "[M]-",
                       "^M\\+H\\-H2O$" = "[M+H-H2O]+",
                       "\\]\\+\\+$" = "\\]2+", # ++ --> 2+
                       "\\]\\-\\-$" = "\\]2-", # -- --> 2-
                       "\\]\\+\\*$" = "\\]\\+", # +* (radical) --> +
                       "\\]\\-\\*$" = "\\]\\-", # -* (radical) --> -
                       "[\\-\\+]{1}ACN" = "C2H3N", # ACN
                       "[\\-\\+]{1}FA" = "CH2O2", # formic acid
                       "[\\-\\+]{1}Hac" = "C2H4O2", # acetic acid
                       "[\\-\\+]{1}DMSO" = "C2H6OS" # DMSO
    )
    for (i in seq_along(adductMapping))
        lib$records[, Precursor_type := sub(names(adductMapping)[i], adductMapping[i], Precursor_type)]

    # UNDONE: more checks (e.g. doesn't detect nonexistent elements)
    printf("Verify/Standardize adducts\n")
    lib$records[!is.na(Precursor_type), Precursor_type := normalizeAdducts(Precursor_type, err = FALSE)]
    
    if (!isFALSE(potAdducts))
    {
        printf("Guessing missing adducts\n")
        # UNDONE: include lib adducts optionally
        if (is.null(potAdducts))
        {
            potAdducts <- unique(c(GenFormAdducts()$adduct_generic, MetFragAdducts()$adduct_generic,
                                   lib$records$Precursor_type))
            potAdducts <- potAdducts[!is.na(potAdducts)]
        }
        
        potAdducts <- lapply(potAdducts, checkAndToAdduct, .var.name = "potAdducts")
        potAdductsPos <- potAdducts[sapply(potAdducts, slot, "charge") > 0]
        potAdductsNeg <- potAdducts[sapply(potAdducts, slot, "charge") < 0]
        lib$records[is.na(Precursor_type) & !is.na(ExactMass) & !is.na(PrecursorMZ) & !is.na(Ion_mode),
                    Precursor_type := withProg(.N, FALSE, mapply(ExactMass, PrecursorMZ, Ion_mode, FUN = function(em, pmz, im)
        {
            pa <- if (im == "POSITIVE") potAdductsPos else potAdductsNeg
            calcMZs <- calculateMasses(em, pa, "mz", err = FALSE) # set err to FALSE to ignore invalid adducts
            wh <- which(numLTE(abs(calcMZs - pmz), absMzDev))
            doProgress()
            # NOTE: multiple hits are ignored (=NA)
            return(if (length(wh) == 1) as.character(pa[[wh]]) else NA_character_)
        }))]
    }
    
    printf("Calculating missing precursor m/z values\n")
    lib$records[is.na(PrecursorMZ) & !is.na(ExactMass) & !is.na(Precursor_type),
                PrecursorMZ := withProg(.N, FALSE, mapply(ExactMass, Precursor_type, FUN = function(em, pt)
    {
        add <- tryCatch(as.adduct(pt), error = function(...) NULL)
        ret <- if (is.null(add)) NA_real_ else em + calculateMasses(em, add, type = "mz")
        doProgress()
        return(ret)
    }))]
    
    return(MSLibrary(records = lib$records[], spectra = lib$spectra, algorithm = "msp"))
}
