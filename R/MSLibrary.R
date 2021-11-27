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


loadMSPLibrary <- function(file, parseComments = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFileExists(file, "r", add = ac)
    aapply(checkmate::assertFlag, . ~ parseComments, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    printf("Parsing file... ")
    lib <- readMSP(normalizePath(file), parseComments)
    printf("Done!\n")
    
    printf("Converting to tables\n")
    lib$records <- as.data.table(lib$records)
    lib$spectra <- withProg(length(lib$spectra), FALSE, lapply(lib$spectra, function(s)
    {
        ret <- as.data.table(s)
        doProgress()
        return(ret)
    }))
    
    # C++ code sets "NA" as string, convert to NA
    for (j in seq_along(lib$records))
        set(lib$records, which(lib$records[[j]] == "NA"), j, NA_character_)
    
    # Ensure case of column names used by patRoon are consistent
    chCols <- c("Name", "SMILES", "InChI", "InChIKey", "Formula", "Precursor_Type", "Ion_mode")
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
    
    printf("Calculating missing formulas\n")
    lib$records[is.na(Formula) & !is.na(SMILES), Formula := withProg(.N, FALSE, sapply(SMILES, function(smi)
    {
        ret <- convertToFormulaBabel(smi, "smi", FALSE)
        doProgress()
        return(fifelse(nzchar(ret), ret, NA_character_))
    }))]

    printf("Calculating missing exact masses\n")
    lib$records[is.na(ExactMass) & !is.na(Formula), ExactMass := withProg(.N, FALSE, sapply(Formula, function(form)
    {
        ret <- tryCatch(getFormulaMass(form), error = function(...) NA_real_)
        doProgress()
        return(ret)
    }))]

    printf("Calculating missing precursor m/z values\n")
    lib$records[is.na(PrecursorMZ) & !is.na(ExactMass) & !is.na(Precursor_Type),
                PrecursorMZ := withProg(.N, FALSE, mapply(ExactMass, Precursor_Type, FUN = function(em, pt)
    {
        add <- tryCatch(as.adduct(pt), error = function(...) NULL)
        ret <- if (is.null(add)) NA_real_ else em + adductMZDelta(add)
        doProgress()
        return(ret)
    }))]
    
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
                       "[\\-\\+]{1}Hac" = "C2H4O2" # acetic acid
    )
    for (i in seq_along(adductMapping))
        lib$records[, Precursor_Type := sub(names(adductMapping)[i], adductMapping[i], Precursor_Type)]

    printf("Verify/Standardize adducts\n")
    lib$records[!is.na(Precursor_Type), Precursor_Type := withProg(.N, FALSE, sapply(Precursor_Type, function(pt)
    {
        ret <- tryCatch(as.character(as.adduct(pt)), error = function(...) NA_character_)
        doProgress()
        return(ret)
    }))]
    
    printf("Guessing missing adducts\n")
    potAdducts <- copy(GenFormAdducts()) # UNDONE: make optional?
    potAdducts <- unique(potAdducts, by = "adduct_generic")
    potAdducts <- potAdducts[molMult == 1] # UNDONE
    potAdducts[, delta := sapply(adduct_generic, function(a) adductMZDelta(as.adduct(a)))]
    lib$records[is.na(Precursor_Type) & !is.na(ExactMass) & !is.na(PrecursorMZ) & !is.na(Ion_mode),
                Precursor_Type := mapply(ExactMass, PrecursorMZ, Ion_mode, FUN = function(em, pmz, im)
    {
        d <- pmz - em
        pa <- potAdducts[((im == "POSITIVE" & charge > 0) |
                              (im == "NEGATIVE" & charge < 0)) & (abs(delta - d) <= 0.002)]$adduct_generic
        # NOTE: multiple hits are ignored (=NA)
        return(if (length(pa) == 1) pa else NA_character_)
    })]

    return(MSLibrary(records = lib$records[], spectra = lib$spectra, algorithm = "msp"))
}
