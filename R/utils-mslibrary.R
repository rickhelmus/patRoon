# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

makeDBIDsUnique <- function(records)
{
    if (anyDuplicated(records$DB_ID))
    {
        records[, DB_ID_orig := DB_ID]
        records[, DB_ID := make.unique(DB_ID, sep = "-")]
    }
    return(records)
}

sanitizeMSLibrary <- function(lib, prefCalcChemProps, neutralChemProps, potAdducts, potAdductsLib, absMzDev, calcSPLASH)
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

    if (nrow(lib$records) == 0)
    {
        lib$spectra <- makeEmptyListNamed(lib$spectra)
        lib$records <- data.table(DB_ID = character())
        return(lib)
    }

    lib$records <- lib$records[DB_ID %chin% names(lib$spectra)]
    
    # C++ code sets "NA" as string, convert to NA. Similarly, library may have 'n/a' markers...
    for (j in seq_along(lib$records))
        set(lib$records, which(lib$records[[j]] %chin% c("NA", "n/a", "N/A")), j, NA_character_)
    
    # Ensure case of column names used by patRoon are consistent
    chCols <- c("Name", "SMILES", "InChI", "InChIKey", "formula", "Precursor_type", "Ion_mode", "Spectrum_type",
                "SPLASH")
    numCols <- c("neutralMass", "PrecursorMZ")
    allCols <- c(chCols, numCols)
    colInds <- match(tolower(allCols), tolower(names(lib$records)))
    setnames(lib$records, colInds[!is.na(colInds)], allCols[!is.na(colInds)])
    
    # ExactMass --> neutralMass (for consistency with other patRoon code)
    emInd <- match("exactmass", tolower(names(lib$records)))
    if (!is.na(emInd) && !"neutralMass" %in% names(lib$records))
        setnames(lib$records, emInd, "neutralMass")
    
    if (!all(c("Name", "DB_ID") %in% names(lib$records)))
        stop("Library file misses mandatory Name and/or DB# data.")
    
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
    printf("Done!\n")
    
    lib$records <- prepareChemTable(lib$records, prefCalcChemProps = prefCalcChemProps,
                                    neutralChemProps = neutralChemProps)
    
    # normalize polarity: ensure uppercase, sometimes shortened as P/N
    lib$records[, Ion_mode := toupper(Ion_mode)]
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
        
        if (isTRUE(potAdducts))
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
