# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

isValidMol <- function(mol) !is.null(mol) # && !is.na(mol)
emptyMol <- function() rcdk::parse.smiles("")[[1]]
isEmptyMol <- function(mol) rcdk::get.atom.count(mol) == 0

getMoleculesFromSMILES <- function(SMILES, doTyping = FALSE, doIsotopes = FALSE, emptyIfFails = FALSE,
                                   SMILESParser = NULL)
{
    # vectorization doesn't work if any of the SMILES are not OK
    # mols <- rcdk::parse.smiles(SMILES)
    mols <- lapply(SMILES, function(sm)
    {
        ret <- rcdk::parse.smiles(sm, smiles.parser = SMILESParser)[[1]]
        if (!isValidMol(ret))
            ret <- rcdk::parse.smiles(sm, kekulise = FALSE, smiles.parser = SMILESParser)[[1]] # might work w/out kekulization
        if (!isValidMol(ret))
        {
            warning(paste("Failed to parse SMILES:", sm))
            if (emptyIfFails)
                ret <- emptyMol()
        }
        else if (!isEmptyMol(ret))
        {
            if (doTyping)
            {
                rcdk::set.atom.types(ret)
                rcdk::do.aromaticity(ret)
            }
            if (doIsotopes)
                rcdk::do.isotopes(ret)
        }
        return(ret)
    })
    
    return(mols)
}

getNeutralMassFromSMILES <- function(SMILES, mustWork = TRUE)
{
    # Use OpenBabel to convert to (neutral) formulas and get mass from there,
    # which seems less buggy compared to parsing the SMILES with RCDK and then
    # using rcdk::get.mol2formula() to get the mass. The latter could give
    # troubles in tryCatch() calls on Linux with certain JVMs.
    forms <- babelConvert(SMILES, "smi", "formula", mustWork = mustWork)
    return(sapply(forms, function(f) if (length(f) > 0) getFormulaMass(f) else NA_character_,
                  USE.NAMES = FALSE))
}

distSMILES <- function(SMI1, SMI2, fpType, fpSimMethod)
{
    mols <- getMoleculesFromSMILES(list(SMI1, SMI2), doTyping = TRUE, emptyIfFails = TRUE)
    fps <- lapply(mols, rcdk::get.fingerprint, type = fpType)
    return(fingerprint::distance(fps[[1]], fps[[2]], fpSimMethod))
}

babelConvert <- function(input, inFormat, outFormat, appendFormula = FALSE, mustWork = TRUE, extraOpts = NULL)
{
    # NOTE: this functions supports formula, logP and MW outFormat as special cases

    if (outFormat == "smiles")
        outFormat <- "smi" # to make checks below easier
    
    stopOrWarn <- function(msg) do.call(if (mustWork) stop else warning, list(msg, call. = FALSE))

    input <- trimws(input)
    input[!nzchar(input)] <- NA_character_
    
    # NOTE: input with spaces should be removed as they will mess up parsing whitespace separated output of obabel    
    invalidInds <- which(!is.na(input) & grepl("[[:space:]]+", input))
    if (length(invalidInds))
    {
        for (i in invalidInds)
            stopOrWarn(sprintf("Failed to convert %d from %s to %s: input has whitespace", i, inFormat, outFormat))
        input[invalidInds] <- NA_character_
    }
        
    input[!nzchar(input)] <- NA_character_ # UNDONE: should we ignore empty strings or not?
    
    # only do non-NA and unique input
    indsToDo <- which(!is.na(input) & !duplicated(input))
    
    finishHandler <- function(cmd)
    {
        if (file.size(cmd$outFile) == 0)
            return(NULL)
        
        # read tab/space separated results
        
        if (outFormat == "smi" && appendFormula)
        {
            # NOTE: the output will be <SMILES><tab><index><space><formula> but fread can only handle one type of
            # separator. So we first read with a \t separation. This gets us in column 1 the SMILES and 2 the space
            # separated index+formula pairs. The latter is then split by another fread call with <space> as separator.
            
            ret <- fread(cmd$outFile, sep = "\t", header = FALSE)
            # BUG: fread thinks single string without newline for text arg is file name, just append newline
            if (nrow(ret) == 1)
                ret[[2]] <- paste0(ret[[2]], "\n")
            ret <- cbind(ret[, 1], fread(text = ret[[2]], sep = " ", header = FALSE))
        }
        else
            ret <- fread(cmd$outFile, sep = if (outFormat == "smi") "\t" else " ", header = FALSE)
        
        # NOTE: with formula/logP/MW output the the index is the first column
        if (outFormat %in% c("formula", "logP", "MW"))
            setnames(ret, 1, "index")
        else
            setnames(ret, seq_len(2), c("result", "index"))
        
        if (appendFormula || outFormat == "formula")
        {
            setnames(ret, ncol(ret), "formula")
            ret[, formula := sub("[\\+\\-]+$", "", formula)] # remove trailing positive/negative charge if present
            ret[, formula := gsub("D", "[2]H", formula, fixed = TRUE)] # handle deuteriums
        }
        else if (outFormat %in% c("logP", "MW"))
            setnames(ret, ncol(ret), outFormat)
        
        return(ret)
    }
    prepareHandler <- function(cmd)
    {
        inFile <- tempfile("obabel_in")
        fwrite(cmd$input, inFile, col.names = FALSE, sep = " ")
        outFile <- tempfile("obabel_out")
        args <- c(cmd$args, c(paste0("-i", inFormat), inFile, "-O", outFile))
        if (appendFormula || outFormat == "formula")
            args <- c(args, "--append", "formula")
        else if (outFormat %in% c("logP", "MW"))
            args <- c(args, "--append", outFormat)
        if (!is.null(extraOpts))
            args <- c(args, extraOpts)
        return(modifyList(cmd, list(command = getExtDepPath("openbabel"), args = args, outFile = outFile)))
    }
    errorHandler <- function(cmd, ...)
    {
        # UNDONE: actually check the exit code to see if the error was a segfault?
        indsToDo <<- c(indsToDo, cmd$input[[2]])
        FALSE
    }
    
    # NOTE: both the input and output is tagged with indices, which makes it much easier to see which conversions failed
    # see https://github.com/openbabel/openbabel/issues/2231
    
    mainArgs <- c(paste0("-o", if (outFormat %in% c("formula", "logP", "MW")) "txt" else outFormat), "-e")
    if (inFormat == "inchi")
        mainArgs <- c(mainArgs, "-an")
    if (outFormat == "inchi" || outFormat == "inchikey")
        mainArgs <- c(mainArgs, "-xw", "-xt")

    # NOTE: In rare cases obabel actually crashes on wrong input, instead of skipping it. Since we run in batches, we
    # cannot know which inputs are to blame. Thus, in the error handler we store all IDs of the batches that failed, and
    # keep try again with smaller batches until everything was processed or failed in a one by one run.
    
    mpm <- getOption("patRoon.MP.method", "classic")
    workersn <- if (mpm == "classic") getOption("patRoon.MP.maxProcs") else future::nbrOfWorkers()
    minBatchSize <- 1000 # put a minimum as overhead of creating processes is significant
    batchn <- max(1, min(workersn, round(length(indsToDo) / minBatchSize)))
    resultsList <- list()
    iterations <- 1
    repeat
    {
        batches <- splitInNBatches(indsToDo, batchn)
        indsToDo <- integer()
        cmdQueue <- lapply(seq_along(batches), function(bi)
        {
            b <- batches[[bi]]
            return(list(args = mainArgs, input = data.frame(input[b], b), logFile = sprintf("obabel-batch_%d-iter_%d.txt", bi, iterations)))
        })
        resultsList <- c(resultsList,  executeMultiProcess(cmdQueue, finishHandler = finishHandler,
                                                           prepareHandler = prepareHandler, errorHandler = errorHandler,
                                                           showProgress = FALSE, logSubDir = "obabel"))

        if (length(indsToDo) == 0 || all(lengths(batches) == 1))
            break # finished or tried all inputs one by one
                
        batchn <- batchn + 1 # try with more batches next iteration
        iterations <- iterations + 1
    }
    ret <- rbindlist(resultsList)

    if (nrow(ret) == 0)
    {
        r <- rep(NA_character_, length(input))
        ret <- if (outFormat == "formula")
            data.table(formula = r)
        else if (outFormat == "logP")
            data.table(logP = rep(NA_real_, length(input)))
        else if (outFormat == "MW")
            data.table(logP = rep(NA_character_, length(input)))
        else
            data.table(result = r)
        if (appendFormula)
            ret[, formula := r][]
    }
    else
    {
        # expand table for NA/non-unique input
        ret <- ret[match(input, input[index])]
        ret[, index := NULL][]
    }
    
    failed <- if (outFormat %in% c("formula", "logP", "MW"))
        which(is.na(ret[[outFormat]]) & !is.na(input))
    else
        which(is.na(ret$result) & !is.na(input))
    for (i in failed)
        stopOrWarn(sprintf("Failed to convert %d ('%s') from %s to %s", i, input[i], inFormat, outFormat))

    if (outFormat %in% c("formula", "logP", "MW"))
        return(ret[[outFormat]])
    if (appendFormula)
        return(ret) # return as table
    return(ret$result)
}

calculateLogP <- function(SMILES, method, mustWork)
{
    if (method == "rcdk")
    {
        doCalc <- function(mol)
        {
            rcdk::convert.implicit.to.explicit(mol)
            return(rcdk::get.xlogp(mol))
        }
        
        mols <- getMoleculesFromSMILES(SMILES, emptyIfFails = TRUE)
        ret <- mapply(mols, SMILES, FUN = function(mol, smi)
        {
            if (isEmptyMol(mol))
            {
                msg <- paste("Failed to parse SMILES to calculate XLogP for", smi)
                do.call(if (mustWork) stop else start, list(msg, call. = FALSE))
                return(NA_real_)
            }
            
            if (mustWork)
                return(doCalc(mol))
            return(tryCatch(doCalc(mol), error = function(...) NA_real_))
        })
        
    }
    else # obabel
        ret <- babelConvert(SMILES, "smi", "logP", mustWork = mustWork)
    
    return(ret)
}

prepareChemTable <- function(chemData, prefCalcChemProps, neutralChemProps, verbose = TRUE)
{
    if (verbose)
        printf("Calculating/Validating chemical data... ")
    
    chemData <- copy(chemData)
    smp <- rcdk::get.smiles.parser() # get re-usable instance, which per rcdk docs is faster
    
    # add missing input columns to simplify things a bit
    # trim white space to improve input data
    for (col in c("SMILES", "InChI", "InChIKey", "formula", "neutralMass"))
    {
        if (is.null(chemData[[col]]))
            chemData[, (col) := if (col == "neutralMass") NA_real_ else NA_character_]
        else if (is.character(chemData[[col]]))
            chemData[, (col) := trimws(get(col))]
    }
    
    convertedInChIs <- babelConvert(chemData$SMILES, "smi", "inchi", appendFormula = TRUE, mustWork = FALSE)
    convertedSMILES <- babelConvert(chemData$InChI, "inchi", "smi", appendFormula = TRUE, mustWork = FALSE)

    # clear input data that couldn't be converted: it's probably invalid
    chemData[is.na(convertedInChIs$result), SMILES := NA_character_]
    chemData[is.na(convertedSMILES$result), InChI := NA_character_]

    chemData[, SMILES := fifelse(!is.na(SMILES), SMILES, convertedSMILES$result)]
    chemData[, InChI := fifelse(!is.na(InChI), InChI, convertedInChIs$result)]
    
    convForms <- fifelse(!is.na(convertedSMILES$formula), convertedSMILES$formula, convertedInChIs$formula)
    
    if (neutralChemProps)
    {
        # UNDONE: InChI has its own neutralization algorithm, which seems more strict. Prefer this over OpenBabel, i.e.
        # by taking converted SMILES from InChIs?
        
        # NOTE: suppress warnings, since the "changed" option below only returns molecules that were changed, hence,
        # babelConvert() warns the rest failed to convert.
        suppressWarnings({
            SMILESNeut <- babelConvert(chemData$SMILES, "smi", "smi", appendFormula = TRUE, mustWork = FALSE,
                                       extraOpts = c("--neutralize", "changed"))
            InChINeut <- babelConvert(chemData$InChI, "inchi", "inchi", appendFormula = TRUE, mustWork = FALSE,
                                      extraOpts = c("--neutralize", "changed", "-xT", "nochg"))
        })
        
        # update data that was changed, clearout what needs to be updated
        takeSMIN <- !is.na(SMILESNeut$result)
        chemData[takeSMIN, c("SMILES", "InChIKey", "formula", "neutralMass") :=
                     .(SMILESNeut$result[takeSMIN], NA_character_, NA_character_, NA_real_)]
        takeInChIN <- !is.na(InChINeut$result)
        chemData[takeInChIN, c("InChI", "InChIKey", "formula", "neutralMass") :=
                     .(InChINeut$result[takeInChIN], NA_character_, NA_character_, NA_real_)]
        
        # mark those that were neutralized
        chemData[, molNeutralized := takeSMIN | takeInChIN]
        
        # update so neutralized formula data is used below
        convForms[takeSMIN] <- SMILESNeut$formula[takeSMIN]
        convForms[takeInChIN] <- InChINeut$formula[takeInChIN]
    }

    # clear invalid InChIKeys
    chemData[!is.na(InChIKey) & !grepl("^[[:upper:]]{14}\\-[[:upper:]]{10}\\-[[:upper:]]{1}$", InChIKey),
             InChIKey := NA_character_]
    
    # NOTE: for formula calculation, both RCDK and OpenBabel cannot output formulas with labeled isotopes. However,
    # OpenBabel does output deuterium as D.
    # --> For now don't consider any formulae from isotope labeled SMILES, unless the isotopes are just deuterium
    isSMIOKForFormula <- !is.na(chemData$SMILES) &
        !grepl("\\[[[:digit:]]+[[:upper:]]{1}[[:lower:]]*\\]", gsub("[2H]", "", chemData$SMILES, fixed = TRUE))
    
    takeConvertedIK <- !is.na(chemData$InChI)
    takeConvertedFormula <- !is.na(convForms) & isSMIOKForFormula
    if (!prefCalcChemProps)
    {
        takeConvertedIK <- takeConvertedIK & is.na(chemData$InChIKey)
        takeConvertedFormula <- takeConvertedFormula & is.na(chemData$formula)
    }
    
    chemData[takeConvertedIK, InChIKey := babelConvert(InChI, "inchi", "inchikey", mustWork = FALSE)]
    chemData[, formula := fifelse(takeConvertedFormula, convForms, formula)]
    
    convFormulaEq <- !is.na(chemData$formula) & !is.na(convForms) & chemData$formula == convForms
    
    # NOTE: use by to avoid duplicated calculations.
    # NOTE: For neutral mass calculation prefer SMILES over formulae as it's faster.
    # NOTE: However, for formulae that were not calculated and therefore need to be validated, the neutral mass is taken
    # from the formulae since rcdk::get.formula returns both.
    
    # Get neutral mass from SMILES: only consider if formula is (the same as) calculated (thus not verified below) or
    # unavailable.
    chemData[(prefCalcChemProps | is.na(neutralMass)) & isSMIOKForFormula &
                 (takeConvertedFormula | convFormulaEq | (is.na(formula) & !is.na(SMILES))), neutralMass := {
        ret <- tryCatch(rcdk::get.exact.mass(getMoleculesFromSMILES(SMILES[1], SMILESParser = smp)[[1]]),
                        error = function(...) NA_real_)
        if (is.na(ret) && !is.na(formula[1])) # failed, try from formula
            ret <- tryCatch(getFormulaMass(formula[1]), error = function(...) NA_real_)
        ret
    }, by = "SMILES"]
    
    
    # For non-calculated formulae: set to NA if invalid or normalize the format and get the neutralMass. Only
    # validate/normalize formulae if different than calculated or we lack a neutralMass.
    chemData[!is.na(formula) & !takeConvertedFormula & (!convFormulaEq | is.na(neutralMass)),
             c("formula", "neutralMassCalc") := {
        
        form <- tryCatch(rcdk::get.formula(formula[1]), error = function(...) NULL)
        if (is.null(form))
            list(NA_character_, NA_real_)
        else if (grepl("[", formula[1], fixed = TRUE))
            list(formula[1], form@mass) # don't normalize formulae with isotopes, as RCDK won't preserve them
        else            
            list(form@string, form@mass)
    }, by = "formula"]
    
    chemData[!is.na(neutralMassCalc) & (prefCalcChemProps | is.na(neutralMass)), neutralMass := neutralMassCalc]
    chemData[, neutralMassCalc := NULL]
    
    if (verbose)
        printf("Done!\n")
    
    return(chemData[])
}

predictRespFactorsSMILES <- function(fgSMILESTab, gInfo, calibrants, eluent, organicModifier, pHAq, concUnit)
{
    # UNDONE: OpenBabel references in ref docs
    
    baseHash <- makeHash(calibrants, eluent, organicModifier, pHAq)
    
    fgSMILESTab <- copy(fgSMILESTab)
    fgSMILESTab[, identifier := group]
    fgSMILESTab[, retention_time := gInfo[group, "rts"]]
    fgSMILESTab[, conc_M := NA_real_]
    fgSMILESTab[, area := 1] # NOTE: we set the area to one to easily get the response factor
    fgSMILESTab[, hash := makeHash(baseHash, retention_time, SMILES), by = seq_len(nrow(fgSMILESTab))] 
    
    cachedData <- loadCacheData("RF_SMILES", fgSMILESTab$hash, simplify = FALSE)
    fgSMILESTabTODO <- if (!is.null(cachedData)) fgSMILESTab[!hash %in% names(cachedData)] else fgSMILESTab
    
    MS2QRes <- NULL 
    if (nrow(fgSMILESTabTODO) > 0)
    {
        MS2QRes <- getMS2QuantRes(calibrants, fgSMILESTabTODO, eluent, organicModifier, pHAq, "")
        saveCacheData("MS2QMD", MS2QRes$MD, baseHash)
        setnames(MS2QRes$RFs, c("identifier", "RF_pred"), c("group", "RF_SMILES"))
        for (i in seq_len(nrow(MS2QRes$RFs)))
            saveCacheData("RF_SMILES", MS2QRes$RFs$RF_SMILES[i], fgSMILESTabTODO$hash[i])
    }
    
    if (!is.null(cachedData))
    {
        cachedRFs <- rbindlist(lapply(cachedData, function(cd) data.table(RF_SMILES = cd)), idcol = "hash")
        cachedRFs[, SMILES := fgSMILESTab$SMILES[match(hash, fgSMILESTab$hash)]]
        cachedRFs[, group := fgSMILESTab$group[match(hash, fgSMILESTab$hash)]]
        cachedRFs[, hash := NULL]
        
        if (is.null(MS2QRes))
        {
            MD <- loadCacheData("MS2QMD", baseHash)
            if (is.null(MD))
            {
                warning("Could not find cached calibration data! You may have an old cache file. ",
                        "Please clear any cached data, eg by running: clearCache(\"RF_SMILES\")", call. = FALSE)
                MD <- list()
            }
            MS2QRes <- list(RFs = cachedRFs, MD = MD)
        }
        else
        {
            MS2QRes$RFs <- rbind(MS2QRes$RFs, cachedRFs)
            MS2QRes$RFs <- MS2QRes$RFs[match(fgSMILESTab$SMILES, SMILES)] # sync order
        }
    }

    # NOTE: do unit conversion the last thing, so we can still use cached data if the user merely changed the unit
    MS2QRes$RFs[, MW := babelConvert(SMILES, "smi", "MW", mustWork = TRUE)]
    # NOTE: need to take the inverse before conversion
    MS2QRes$RFs[, RF_SMILES := 1/convertConc(1/RF_SMILES[1], "M", concUnit, MW[1]), by = "SMILES"]
    MS2QRes$RFs[, MW := NULL][]
    
    return(MS2QRes)
}

predictLC50SMILES <- function(SMILES, LC50Mode, concUnit)
{
    inp <- data.table(SMILES = SMILES)
    
    smp <- rcdk::get.smiles.parser()
    inp[, exactMass := tryCatch(rcdk::get.exact.mass(getMoleculesFromSMILES(SMILES[1], emptyIfFails = FALSE,
                                                                            SMILESParser = smp)[[1]]),
                                error = function(...) NA_real_), by = "SMILES"]
    if (any(is.na(inp$exactMass)))
    {
        warning("Ignoring the following SMILES because mass calculation failed: ",
                sprintf("%s (%d)", inp[is.na(exactMass)]$SMILES, which(is.na(inp$exactMass))))
        inp <- inp[!is.na(exactMass)]
    }
    
    if (nrow(inp) == 0)
        return(data.table(SMILES = character(), LC50_SMILES = character()))
    
    hashes <- sapply(inp$SMILES, makeHash, LC50Mode)
    cachedData <- loadCacheData("LC50_SMILES", hashes, simplify = FALSE)
    indsTODO <- if (!is.null(cachedData)) which(!hashes %in% names(cachedData)) else seq_along(hashes)
    hashesTODO <- hashes[indsTODO]
    
    LC50s <- NULL
    if (length(indsTODO) > 0)
    {
        suppressMessages(utils::capture.output(LC50s <- MS2Tox::LC50fromSMILES(inp[indsTODO], LC50Mode)))
        setDT(LC50s)
        setnames(LC50s, "LC50_predicted", "LC50_SMILES")
        for (i in seq_len(nrow(LC50s)))
            saveCacheData("LC50_SMILES", LC50s$LC50_SMILES[i], hashesTODO[i])
    }

    if (!is.null(cachedData))
    {
        cachedLC50s <- rbindlist(lapply(cachedData, function(cd) data.table(LC50_SMILES = cd)), idcol = "hash")
        cachedLC50s[, SMILES := inp$SMILES[match(hash, hashes)]]
        cachedLC50s[, hash := NULL]
        
        if (is.null(LC50s))
            LC50s <- cachedLC50s
        else
        {
            LC50s <- rbind(LC50s, cachedLC50s)
            LC50s <- LC50s[match(inp$SMILES, SMILES)] # sync order
        }
    }

    # NOTE: do unit conversion the last thing, so we can still use cached data if the user merely changed the unit
    LC50s[, MW := babelConvert(SMILES, "smi", "MW", mustWork = TRUE)]
    LC50s[, LC50_SMILES := convertConc(LC50_SMILES[1], "log mM", concUnit, MW[1]), by = "SMILES"]
    LC50s[, MW := NULL]
    
    return(LC50s[])
}
