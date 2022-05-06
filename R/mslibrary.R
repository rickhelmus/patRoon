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

sanitizeMSLibrary <- function(lib, potAdducts, potAdductsLib, absMzDev, calcSPLASH)
{
    printf("Converting to tables... ")
    lib$records <- as.data.table(lib$records)
    lib$spectra <- Map(lib$spectraMZs, lib$spectraInts, lib$annotations, f = function(mz, int, ann)
    {
        # NOTE: use setDT(list(...)) here since it seems to be quite a bit faster than data.table(...)
        if (length(ann) == 0)
            setDT(list(mz = mz, intensity = int))
        else
            setDT(list(mz = mz, intensity = int, annotation = ann))
    })
    printf("Done!\n")
    
    lib$spectra <- pruneList(lib$spectra, checkZeroRows = TRUE)
    lib$records <- lib$records[DB_ID %chin% names(lib$spectra)]
    
    # C++ code sets "NA" as string, convert to NA. Similarly, library may have 'n/a' markers...
    for (j in seq_along(lib$records))
        set(lib$records, which(lib$records[[j]] %chin% c("NA", "n/a", "N/A")), j, NA_character_)
    
    # Ensure case of column names used by patRoon are consistent
    chCols <- c("Name", "SMILES", "InChI", "InChIKey", "formula", "Precursor_type", "Ion_mode", "SPLASH")
    numCols <- c("neutralMass", "PrecursorMZ")
    allCols <- c(chCols, numCols)
    colInds <- match(tolower(allCols), tolower(names(lib$records)))
    setnames(lib$records, colInds[!is.na(colInds)], allCols[!is.na(colInds)])
    
    # ExactMass --> neutralMass (for consistency with other patRoon code)
    emInd <- match("exactmass", tolower(names(lib$records)))
    if (!is.na(emInd) && !"neutralMass" %in% names(lib$records))
        setnames(lib$records, emInd, "neutralMass")
    
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
    
    printf("Clean up formulas... ")
    # remove ion species format ([formula]+/-)
    lib$records[!is.na(formula), formula := gsub("^\\[(.+)\\][[:digit:]]*[\\+\\-]+", "\\1", formula)]
    # remove trailing charges
    lib$records[!is.na(formula), formula := gsub("\\+|\\-$", "", formula)]
    lib$records[!is.na(formula), formula := gsub("\\[|\\]|\\+|\\-", "", formula)] # remove ion species format ([formula]+/-)
    printf("Done!\n")
    
    lib$records <- prepareChemTable(lib$records)
    
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
    
    printf("Verify/Standardize adducts... ")
    lib$records[!is.na(Precursor_type), Precursor_type := normalizeAdducts(Precursor_type, err = FALSE)]
    printf("Done!\n")
    
    if (!isFALSE(potAdducts))
    {
        printf("Guessing missing adducts...\n")
        
        if (is.null(potAdducts))
        {
            potAdducts <- union(GenFormAdducts()$adduct_generic, MetFragAdducts()$adduct_generic)
            if (potAdductsLib)
            {
                potAdducts <- union(potAdducts, lib$records$Precursor_type)
                potAdducts <- potAdducts[!is.na(potAdducts)]
            }
        }
        
        potAdducts <- lapply(potAdducts, checkAndToAdduct, .var.name = "potAdducts")
        potAdductsPos <- potAdducts[sapply(potAdducts, slot, "charge") > 0]
        potAdductsNeg <- potAdducts[sapply(potAdducts, slot, "charge") < 0]
        lib$records[is.na(Precursor_type) & !is.na(neutralMass) & !is.na(PrecursorMZ) & !is.na(Ion_mode),
                    Precursor_type := withProg(.N, FALSE, mapply(neutralMass, PrecursorMZ, Ion_mode, FUN = function(em, pmz, im)
                    {
                        pa <- if (im == "POSITIVE") potAdductsPos else potAdductsNeg
                        calcMZs <- calculateMasses(em, pa, "mz", err = FALSE) # set err to FALSE to ignore invalid adducts
                        wh <- which(numLTE(abs(calcMZs - pmz), absMzDev))
                        doProgress()
                        # NOTE: multiple hits are ignored (=NA)
                        return(if (length(wh) == 1) as.character(pa[[wh]]) else NA_character_)
                    }))]
    }
    
    printf("Calculating missing precursor m/z values...\n")
    lib$records[is.na(PrecursorMZ) & !is.na(neutralMass) & !is.na(Precursor_type),
                PrecursorMZ := withProg(.N, FALSE, mapply(neutralMass, Precursor_type, FUN = function(em, pt)
                {
                    add <- tryCatch(as.adduct(pt), error = function(...) NULL)
                    ret <- if (is.null(add)) NA_real_ else em + calculateMasses(em, add, type = "mz", err = FALSE)
                    doProgress()
                    return(ret)
                }))]
    
    if (calcSPLASH && any(is.na(lib$records$SPLASH)))
    {
        # UNDONE: this is rather slow... MassBank has SPLASH values, so not that important...
        printf("Calculating missing SPLASH values...\n")
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
MSLibrary <- setClass("MSLibrary", slots = c(records = "data.table", spectra = "list"),
                      contains = "workflowStep")

setMethod("initialize", "MSLibrary", function(.Object, ...)
{
    .Object <- callNextMethod(.Object, ...)
    return(.Object)
})

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
    allSpecs <- rbindlist(spectra(x), idcol = "DB_ID", fill = TRUE)
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
    }
    else
    {
        for (rec in i)
        {
            if (is.function(j))
            {
                inds <- j(rec, obj@spectra[[rec]])
                if (is.logical(inds))
                    inds <- which(inds)
            }
            else # j = vector
                inds <- j[j <= nrow(obj@spectra)]
            if (length(inds) > 0)
                obj@spectra[[rec]] <- obj@spectra[[rec]][-inds, , drop = FALSE]
        }
        obj@spectra <- pruneList(obj@spectra, checkZeroRows = TRUE)
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
        noAnnot <- sapply(spectra(obj), function(sp) is.null(sp[["annotation"]]))
        obj <- delete(obj, i = if (negate) !noAnnot else noAnnot)
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
                recs[keep == TRUE, keep := !is.na(neutralMass) & mark(neutralMass %inrange% massRange)]
            
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
setMethod("convertToSuspects", "MSLibrary", function(obj, adduct,
                                                     avgSpecParams = getDefAvgPListParams(minIntensityPre = 0,
                                                                                          minIntensityPost = 2,
                                                                                          topMost = 10),
                                                     collapse = TRUE, suspects = NULL)
{
    adduct <- checkAndToAdduct(adduct)
    
    ac <- checkmate::makeAssertCollection()
    assertAvgPListParams(avgSpecParams, add = ac)
    checkmate::assertFlag(collapse, add = ac)
    if (!is.null(suspects))
        assertSuspectList(suspects, FALSE, FALSE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (length(obj) == 0)
        stop("Cannot create suspect list: no data", call. = FALSE)

    libRecs <- records(obj)
    libRecs <- libRecs[Precursor_type == as.character(adduct)]
    if (nrow(libRecs) == 0)
        stop("No records found (for input adduct)", call. = FALSE)
    libSpecs <- spectra(obj)
    
    getAvgFrags <- function(specs, PrecursorMZ)
    {
        pls <- Map(specs, PrecursorMZ, f = function(sp, pmz)
        {
            sp <- copy(sp)[, c("mz", "intensity"), with = FALSE] # skip annotation, if present
            sp <- assignPrecursorToMSPeakList(sp, pmz)
            return(sp)
        })
        
        avgPL <- averageSpectra(pls, avgSpecParams$clusterMzWindow, avgSpecParams$topMost,
                                avgSpecParams$minIntensityPre, avgSpecParams$minIntensityPost,
                                avgSpecParams$avgFun, avgSpecParams$method, FALSE, avgSpecParams$retainPrecursorMSMS)
        doProgress()
        paste0(avgPL$mz, collapse = ";")
    }
        
    printf("Calculating MS/MS fragments...\n")
    
    if (!is.null(suspects))
    {
        ret <- if (is.data.table(suspects)) copy(suspects) else as.data.table(suspects)
        ret <- prepareSuspectList(ret, NULL, FALSE, calcMZs = FALSE)
        ret[, InChIKey1 := getIKBlock1(InChIKey)]
        
        if (any(is.na(ret$InChIKey1)))
            warning(paste("Ignored the following suspects because no InChIKey1 could be calculated:",
                           getStrListWithMax(ret$name, 10, ", ")))
        
        recs <- copy(libRecs)
        recs <- recs[!is.na(InChIKey) & !is.na(PrecursorMZ)]
        recs[, InChIKey1 := getIKBlock1(InChIKey)]
        
        withProg(uniqueN(ret$InChIKey1), FALSE, ret[, fragments_mz := sapply(InChIKey1, function(IK1)
        {
            recsSub <- recs[InChIKey1 == IK1]
            if (is.na(IK1) || nrow(recsSub) == 0)
                return("")
            return(getAvgFrags(libSpecs[recsSub$DB_ID], recsSub$PrecursorMZ))
        }), by = "InChIKey1"])
        
        ret[, InChIKey1 := NULL]
        
        frCount <- sum(nzchar(ret$fragments_mz))
        printf("Filled in fragments for %d/%d (%.2f%%) suspects\n", frCount, nrow(ret), frCount / nrow(ret) * 100)
    }
    else
    {
        ret <- copy(libRecs)
        
        if (collapse)
        {
            ret <- ret[!is.na(InChIKey) & !is.na(PrecursorMZ)]
            ret[, InChIKey1 := getIKBlock1(InChIKey)]

            frMZ <- withProg(uniqueN(ret$InChIKey1), FALSE, ret[, getAvgFrags(libSpecs[DB_ID], PrecursorMZ),
                                                                by = "InChIKey1"])
            setnames(frMZ, "V1", "fragments_mz")

            ret <- unique(ret, by = "InChIKey1")
            ret <- merge(ret, frMZ, by = "InChIKey1")
        }
        else
        {
            withProg(nrow(ret), FALSE, ret[, fragments_mz := sapply(libSpecs[ret$DB_ID], function(spec)
            {
                doProgress()
                paste0(spec[, "mz"], collapse = ";")
            })])
        }
        
        mapCols <- c(Name = "name",
                     SMILES = "SMILES",
                     InChI = "InChI",
                     InChIKey = "InChIKey",
                     formula = "formula",
                     Precursor_type = "adduct",
                     neutralMass = "neutralMass",
                     fragments_mz = "fragments_mz")
        mapCols <- mapCols[names(mapCols) %in% names(ret)]
        setnames(ret, names(mapCols), mapCols)
        ret <- ret[, mapCols, with = FALSE]
        ret <- prepareSuspectList(ret, NULL, FALSE, FALSE)
    }
    
    return(ret[])
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


loadMSPLibrary <- function(file, parseComments = TRUE, potAdducts = NULL, potAdductsLib = TRUE, absMzDev = 0.002,
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

loadMoNAJSONLibrary <- function(file, potAdducts = NULL, potAdductsLib = TRUE, absMzDev = 0.002, calcSPLASH = TRUE)
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
