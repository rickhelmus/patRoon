#' @include main.R
#' @include workflow-step.R
NULL

makeDBIDsUnique <- function(records)
{
    if (anyDuplicated(records$DB_ID))
    {
        records[, DB_ID_orig := DB_ID]
        records[, DB_ID := make.unique(DB_ID, sep = "-")]
    }
    return(records)
}

sanitizeMSLibrary <- function(lib, potAdducts, absMzDev, calcSPLASH)
{
    printf("Converting to tables\n")
    lib$records <- as.data.table(lib$records)
    # UNDONE: not for now, takes a lot of time for large amounts
    # lib$spectra <- withProg(length(lib$spectra), FALSE, lapply(lib$spectra, function(s)
    # {
    #     ret <- as.data.table(s)
    #     doProgress()
    #     return(ret)
    # }))
    
    lib$spectra <- pruneList(lib$spectra, checkEmptyElements = TRUE)
    lib$annotations <- pruneList(lib$annotations, checkEmptyElements = TRUE)
    lib$annotations <- lib$annotations[names(lib$annotations) %chin% names(lib$spectra)]
    lib$records <- lib$records[DB_ID %chin% names(lib$spectra)]
    
    # C++ code sets "NA" as string, convert to NA. Similarly, library may have 'n/a' markers...
    for (j in seq_along(lib$records))
        set(lib$records, which(lib$records[[j]] %chin% c("NA", "n/a", "N/A")), j, NA_character_)
    
    # Ensure case of column names used by patRoon are consistent
    chCols <- c("Name", "SMILES", "InChI", "InChIKey", "Formula", "Precursor_type", "Ion_mode", "SPLASH")
    numCols <- c("ExactMass", "PrecursorMZ")
    allCols <- c(chCols, numCols)
    change <- match(tolower(allCols), tolower(names(lib$records)), nomatch = integer())
    setnames(lib$records, change, allCols)
    
    if (!all(c("Name", "DB_ID") %in% names(lib$records)))
        stop("MSP file misses mandatory Name and/or DB# data.")
    
    lib$records <- makeDBIDsUnique(lib$records)
    names(lib$spectra) <- lib$records$DB_ID
    
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
    suppressWarnings({
        lib$records[, (numCols) := lapply(.SD, as.numeric), .SDcols = numCols]
        if (!is.null(lib$records[["MW"]]))
            lib$records[, MW := as.numeric(MW)] # MW is not used by patRoon, but still make it numeric if present
    })
    
    # UNDONE: this is quite slow, skip for now?
    # printf("Sanitizing SMILES/InChI values... ")
    # lib$records[!is.na(SMILES), SMILES := babelConvert(SMILES, "smi", "smi", FALSE)]
    # lib$records[!is.na(InChI), InChI := babelConvert(InChI, "inchi", "inchi", FALSE)]
    # printf("Done!\n")
    
    printf("Clean up formulas...\n")
    lib$records[!is.na(Formula), Formula := gsub("\\[|\\]|\\+|\\-", "", Formula)] # remove ion species format ([formula]+/-)
    printf("Clearing invalid formulas...\n")
    lib$records[!verifyFormulas(Formula), Formula := NA_character_]
    
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
    
    printf("Verify/Standardize adducts\n")
    lib$records[!is.na(Precursor_type), Precursor_type := normalizeAdducts(Precursor_type, err = FALSE)]
    
    if (!isFALSE(potAdducts))
    {
        printf("Guessing missing adducts\n")
        
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
    
    if (calcSPLASH && any(is.na(lib$records$SPLASH)))
    {
        # UNDONE: this is rather slow... MassBank has SPLASH values, so not that important...
        printf("Calculating missing SPLASH values\n")
        checkPackage("splashR", "berlinguyinca/spectra-hash", "splashR")
        lib$records[is.na(SPLASH), SPLASH := withProg(.N, FALSE, sapply(lib$spectra, function(sp)
        {
            doProgress()
            splashR::getSplash(sp)
        }))]
    }
    
    return(lib)
}


#' @export
MSLibrary <- setClass("MSLibrary", slots = c(records = "data.table", spectra = "list", annotations = "list"),
                      contains = "workflowStep")

setMethod("initialize", "MSLibrary", function(.Object, ...)
{
    .Object <- callNextMethod(.Object, ...)
    .Object@annotations <- makeEmptyListNamed(.Object@annotations)
    return(.Object)
})

#' @export
setMethod("records", "MSLibrary", function(obj) obj@records)

#' @export
setMethod("spectra", "MSLibrary", function(obj) obj@spectra)

#' @export
setMethod("annotations", "MSLibrary", function(obj) obj@annotations)

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
        x <- delete(x, setdiff(names(x), i))
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
setMethod("delete", "MSLibrary", function(obj, i = NULL, j = NULL, ...)
{
    if (!is.function(i))
        i <- assertDeleteArgAndToChr(i, names(obj))
    checkmate::assert(
        checkmate::checkFunction(j, null.ok = TRUE),
        checkmate::checkIntegerish(j, any.missing = FALSE, null.ok = TRUE),
        .var.name = "j"
    )
    if (is.function(i) && is.function(j))
        stop("i and j cannot both be functions")

    if (length(i) == 0 || (!is.null(j) && length(j) == 0))
        return(obj) # nothing to remove...
    
    # i = vector remove specified records
    # i = function: remove returned records
    # j = vector/function: index of peaks to remove from spectra of i (all if i=NULL)
    
    if (is.null(j))
    {
        if (!is.function(i))
            obj@records <- obj@records[!DB_ID %chin% i]
        else
        {
            rm <- i(obj@records, ...)
            if (is.logical(rm))
                obj@records <- obj@records[!rm]
            else if (is.character(rm))
                obj@records <- obj@records[!rm %chin% DB_ID]
            else
                obj@records <- obj@records[setdiff(seq_len(nrow(obj@records)), rm)]
        }
        
        obj@spectra <- obj@spectra[obj@records$DB_ID]
        if (length(obj@annotations) > 0)
            obj@annotations <- obj@annotations[names(obj@annotations) %chin% obj@records$DB_ID]
    }
    else
    {
        for (rec in i)
        {
            if (is.function(j))
            {
                inds <- j(rec, obj@spectra[[rec]], obj@annotations[[rec]])
                if (is.logical(inds))
                    inds <- which(inds)
            }
            else # j = vector
                inds <- j[j <= nrow(obj@spectra)]
            if (length(inds) > 0)
            {
                obj@spectra[[rec]] <- obj@spectra[[rec]][-inds, , drop = FALSE]
                if (!is.null(obj@annotations[[rec]]))
                    obj@annotations[[rec]] <- obj@annotations[[rec]][-inds]
            }
        }
        obj@spectra <- pruneList(obj@spectra, checkEmptyElements = TRUE)
        obj@annotations <- pruneList(obj@annotations, checkEmptyElements = TRUE)
        obj@records <- obj@records[DB_ID %chin% names(obj@spectra)]
    }
    
    return(obj)
})

#' @export
setMethod("filter", "MSLibrary", function(obj, properties = NULL, massRange = NULL, mzRangeSpec = NULL,
                                          relMinIntensity = NULL, topMost = NULL, onlyAnnotated = FALSE,
                                          negate = FALSE)
{
    
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertList(properties, any.missing = FALSE, null.ok = TRUE, add = ac)
    if (!is.null(properties))
    {
        checkmate::assertNames(names(properties), type = "unique", subset.of = names(records(obj)), add = add)
        checkmate::qassertr(properties, "V")
    }
    aapply(assertRange, . ~ massRange + mzRangeSpec, null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertNumber(relMinIntensity, null.ok = TRUE, add = ac)
    checkmate::assertCount(topMost, positive = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertFlag(onlyAnnotated, add = ac)
    checkmate::assertFlag(negate, add = ac)
    checkmate::reportAssertions(ac)
    
    mark <- if (negate) function(x) !x else function(x) x

    if (onlyAnnotated)
    {
        if (negate)
            obj <- delete(obj, i = intersect(names(obj), names(annotations(obj))))
        else
            obj <- delete(obj, i = setdiff(names(obj), names(annotations(obj))))
    }
    
    if (!is.null(properties) || !is.null(massRange))
    {
        obj <- delete(obj, i = function(recs, ...)
        {
            recs <- copy(recs)
            recs[, keep := TRUE]
            
            if (!is.null(properties) && length(properties) > 0)
            {
                for (prop in names(properties))
                    recs[keep == TRUE, keep := mark(get(prop) %in% properties[[prop]])]
            }
            
            if (!is.null(massRange))
                recs[keep == TRUE, keep := !is.na(ExactMass) & mark(ExactMass %inrange% massRange)]
            
            return(!recs$keep)
        })
    }
    if (!is.null(mzRangeSpec) || !is.null(relMinIntensity))
    {
        obj <- delete(obj, j = function(rec, spec, ...)
        {
            keep <- rep(TRUE, nrow(spec))
            if (!is.null(mzRangeSpec))
                keep <- mark(spec[, "mz"] %inrange% mzRangeSpec)
            if (!is.null(relMinIntensity))
            {
                relInts <- spec[, "intensity"] / max(spec[, "intensity"])
                keep <- keep & mark(numGTE(relInts, relMinIntensity))
            }
            return(!keep)
        })
    }
    if (!is.null(topMost)) # NOTE: do after previous filters
    {
        if (!is.null(mzRangeSpec) || !is.null(relMinIntensity) || !is.null(topMost))
        obj <- delete(obj, j = function(rec, spec, ...)
        {
            err <- tryCatch(rep(TRUE, nrow(spec)), error=function(...) NULL)
            if (is.null(err)) browser()
            keep <- rep(TRUE, nrow(spec))
            if (nrow(spec) > topMost)
            {
                ord <- order(spec[, "intensity"], decreasing = !negate)
                keep <- seq_len(nrow(spec)) %in% ord[seq_len(topMost)] # NOTE: keep order
            }
            return(!keep)
        })
    }
    
    return(obj)
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

setMethod("merge", c("MSLibrary", "MSLibrary"), function(x, y, ...)
{
    if (length(x) == 0)
        return(y)
    else if (length(y) == 0)
        return(x)
    
    if (any(is.na(records(x)$SPLASH)) || any(is.na(records(y)$SPLASH)))
        stop("x/y has missing SPLASH values. Please load the library with calcSPLASH=TRUE")
    
    unY <- records(y)[!SPLASH %chin% records(x)$SPLASH]$DB_ID
    y <- y[unY]
    
    recordsAll <- rbind(records(x), records(y), fill = TRUE)
    specsAll <- c(spectra(x), spectra(y))
    
    recordsAll <- makeDBIDsUnique(recordsAll)
    names(specsAll) <- recordsAll$DB_ID
    
    return(MSLibrary(records = recordsAll[], spectra = specsAll, algorithm = "merged"))
})


loadMSPLibrary <- function(file, parseComments = TRUE, potAdducts = NULL, absMzDev = 0.002, calcSPLASH = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFileExists(file, "r", add = ac)
    aapply(checkmate::assertFlag, . ~ parseComments + calcSPLASH, fixed = list(add = ac))
    checkmate::assert(checkmate::checkNull(potAdducts),
                      checkmate::checkFALSE(potAdducts),
                      checkmate::checkCharacter(potAdducts, any.missing = FALSE, min.chars = 1),
                      checkmate::checkList(potAdducts, types = c("adduct", "character"), any.missing = FALSE),
                      .var.name = "potAdducts")
    checkmate::assertNumber(absMzDev, lower = 0, finite = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    lib <- readMSP(normalizePath(file), parseComments)
    lib <- sanitizeMSLibrary(lib, potAdducts, absMzDev, calcSPLASH)
    
    return(MSLibrary(records = lib$records[], spectra = lib$spectra, algorithm = "msp"))
}

loadMoNAJSONLibrary <- function(file, potAdducts = NULL, absMzDev = 0.002, calcSPLASH = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFileExists(file, "r", add = ac)
    checkmate::assertNumber(absMzDev, lower = 0, finite = TRUE, add = ac)
    checkmate::assertFlag(calcSPLASH, add = ac)
    checkmate::reportAssertions(ac)
    
    lib <- readMoNAJSON(normalizePath(file))
    lib <- sanitizeMSLibrary(lib, potAdducts, absMzDev, calcSPLASH)
    
    return(MSLibrary(records = lib$records[], spectra = lib$spectra, annotations = lib$annotations, algorithm = "json"))
}
