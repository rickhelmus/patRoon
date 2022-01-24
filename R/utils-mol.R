isValidMol <- function(mol) !is.null(mol) # && !is.na(mol)
emptyMol <- function() rcdk::parse.smiles("")[[1]]
isEmptyMol <- function(mol) rcdk::get.atom.count(mol) == 0

getMoleculesFromSMILES <- function(SMILES, doTyping = FALSE, doIsotopes = FALSE, emptyIfFails = FALSE)
{
    # vectorization doesn't work if any of the SMILES are not OK
    # mols <- rcdk::parse.smiles(SMILES)
    mols <- lapply(SMILES, function(sm)
    {
        ret <- rcdk::parse.smiles(sm)[[1]]
        if (!isValidMol(ret))
            ret <- rcdk::parse.smiles(sm, kekulise = FALSE)[[1]] # might work w/out kekulization
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
    # NOTE: this functions supports formula outFormat as special case

    if (outFormat == "smiles")
        outFormat <- "smi" # to make checks below easier
        
    input[!nzchar(input)] <- NA_character_ # UNDONE: should we ignore empty strings or not?
    
    # only do non-NA and unique input
    indsToDo <- which(!is.na(input) & !duplicated(input))
    
    mpm <- getOption("patRoon.MP.method", "classic")
    batchn <- if (mpm == "classic") getOption("patRoon.MP.maxProcs") else future::nbrOfWorkers()
    minBatchSize <- 1000 # put a minimum as overhead of creating processes is significant
    batchn <- max(1, min(batchn, round(length(indsToDo) / minBatchSize)))
    
    batches <- splitInNBatches(indsToDo, batchn)
    
    # NOTE: both the input and output is tagged with indices, which makes it much easier to see which conversions failed
    # see https://github.com/openbabel/openbabel/issues/2231
    
    mainArgs <- c(paste0("-o", if (outFormat == "formula") "txt" else outFormat), "-e")
    if (inFormat == "inchi")
        mainArgs <- c(mainArgs, "-an")
    if (outFormat == "inchi")
        mainArgs <- c(mainArgs, "-xw", "-xt")
    cmdQueue <- lapply(seq_along(batches), function(bi)
    {
        b <- batches[[bi]]
        return(list(args = mainArgs, input = data.frame(input[b], b), logFile = paste0("obabel-batch_", bi, ".txt")))
    })
    
    resultsList <- executeMultiProcess(cmdQueue, finishHandler = function(cmd)
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
            ret <- cbind(ret[, 1], fread(text = ret[[2]], sep = " ", header = FALSE))
        }
        else
            ret <- fread(cmd$outFile, sep = if (outFormat == "smi") "\t" else " ", header = FALSE)
        
        # NOTE: with formula output the the index is the first column
        if (outFormat == "formula")
            setnames(ret, 1, "index")
        else
            setnames(ret, seq_len(2), c("result", "index"))
        
        if (appendFormula || outFormat == "formula")
        {
            setnames(ret, ncol(ret), "formula")
            ret[, formula := sub("[\\+\\-]+$", "", formula)] # remove trailing positive/negative charge if present
        }
        
        return(ret)
    }, prepareHandler = function(cmd)
    {
        inFile <- tempfile("obabel_in")
        fwrite(cmd$input, inFile, col.names = FALSE, sep = " ")
        outFile <- tempfile("obabel_out")
        args <- c(cmd$args, c(paste0("-i", inFormat), inFile, "-O", outFile))
        if (appendFormula || outFormat == "formula")
            args <- c(args, "--append", "formula")
        if (!is.null(extraOpts))
            args <- c(args, extraOpts)
        return(modifyList(cmd, list(command = getCommandWithOptPath("obabel", "obabel"),
                                    args = args, outFile = outFile)))
    }, showProgress = FALSE, logSubDir = "obabel")
    
    ret <- rbindlist(resultsList)
    
    if (nrow(ret) == 0)
    {
        r <- rep(NA_character_, length(input))
        ret <- if (outFormat == "formula") data.table(formula = r) else data.table(result = r)
        if (appendFormula)
            ret[, formula := r][]
    }
    else
    {
        # expand table for NA/non-unique input
        ret <- ret[match(input, input[index])]
        ret[, index := NULL][]
    }
    
    stopOrWarn <- function(msg) do.call(if (mustWork) stop else warning, list(msg, call. = FALSE))
    failed <- if (outFormat == "formula")
        which(is.na(ret$formula) & !is.na(input))
    else
        which(is.na(ret$result) & !is.na(input))
    for (i in failed)
        stopOrWarn(sprintf("Failed to convert %d ('%s') from %s to %s", i, input[i], inFormat, outFormat))

    if (outFormat == "formula")
        return(ret$formula)
    else if (appendFormula)
        return(ret) # return as table
    return(ret$result)
}

convertChemDataIfNeeded <- function(tab, destFormat, destCol, fromFormats, fromCols)
{
    hasData <- function(x) !is.na(x) & nzchar(x)
    missingInTab <- function(x) if (is.null(tab[[x]])) rep(TRUE, nrow(tab)) else !hasData(tab[[x]])
    
    countEntries <- function() if (is.null(tab[[destCol]])) 0 else sum(hasData(tab[[destCol]]))
    curEntryCount <- countEntries()
    if (curEntryCount < nrow(tab))
    {
        printf("Trying to calculate missing %s data... ", destCol)
        
        for (i in seq_along(fromFormats))
        {
            if (!is.null(tab[[fromCols[i]]]))
                tab[missingInTab(destCol) & !missingInTab(fromCols[i]),
                    (destCol) := babelConvert(get(fromCols[i]), fromFormats[i], destFormat, mustWork = FALSE)]
        }
        
        newEntryCount <- countEntries() - curEntryCount
        printf("Done! Filled in %d (%.1f%%) entries.\n", newEntryCount,
               if (newEntryCount > 0) newEntryCount * 100 / nrow(tab) else 0)
    }
    return(tab)
}

calculateXLogP <- function(SMILES, mustWork)
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
    
    return(ret)
}

prepareChemTable <- function(chemData)
{
    chemData <- copy(chemData)
    
    # UNDONE: skip missing input columns
    
    convertedInChIs <- babelConvert(chemData$SMILES, "smi", "inchi", appendFormula = TRUE, mustWork = FALSE)
    convertedSMILES <- babelConvert(chemData$InChI, "inchi", "smi", appendFormula = TRUE, mustWork = FALSE)
    
    # clear input data that couldn't be converted: it's probably invalid
    chemData[is.na(convertedInChIs$result), SMILES := NA_character_]
    chemData[is.na(convertedSMILES$result), InChI := NA_character_]
    
    # prefer calculated SMILES/InChIs/InChIKeys where possible
    chemData[, SMILES := fifelse(!is.na(convertedSMILES$result), convertedSMILES$result, SMILES)]
    chemData[, InChI := fifelse(!is.na(convertedInChIs$result), convertedInChIs$result, InChI)]
    chemData[!is.na(InChI), InChIKey := babelConvert(InChI, "inchi", "inchikey", mustWork = FALSE)]
    
    # clear invalid InChIKeys
    chemData[!is.na(InChIKey) & !grepl("^[[:upper:]]{14}\\-[[:upper:]]{9}\\-[[:upper:]]{1}$", InChIKey),
             InChIKey := NA_character_]
    
    # prefer calculated formulas
    chemData[, formula := fifelse(!is.na(convertedSMILES$formula), convertedSMILES$formula, formula)]
    chemData[, formula := fifelse(!is.na(convertedInChIs$formula), convertedInChIs$formula, formula)]
    
    # load non-calculated formulae to (1) NA if invalid and (2) normalize the format
    # NOTE: use by to avoid duplicated calculations
    chemData[is.na(SMILES) & !is.na(formula), formula := 
    {
        tryCatch(rcdk::get.formula(formula[1]), error = function(...) NA_character_)
    }, by = "formula"]
    
    # prefer calculated masses
    chemData[!is.na(formula), neutralMass := sapply(formula, getFormulaMass)]
    
    return(chemData)
}
