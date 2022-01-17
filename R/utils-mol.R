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
    forms <- convertToFormulaBabel(SMILES, "smi", mustWork = mustWork)
    return(sapply(forms, function(f) if (length(f) > 0) getFormulaMass(f) else NA_character_,
                  USE.NAMES = FALSE))
}

distSMILES <- function(SMI1, SMI2, fpType, fpSimMethod)
{
    mols <- getMoleculesFromSMILES(list(SMI1, SMI2), doTyping = TRUE, emptyIfFails = TRUE)
    fps <- lapply(mols, rcdk::get.fingerprint, type = fpType)
    return(fingerprint::distance(fps[[1]], fps[[2]], fpSimMethod))
}

babelConvert <- function(input, inFormat, outFormat, mustWork = TRUE, extraOpts = NULL)
{
    if (length(input) == 0)
        return(character())

    indsNoNA <- which(!is.na(input))
    
    mpm <- getOption("patRoon.MP.method", "classic")
    batchn <- if (mpm == "classic") getOption("patRoon.MP.maxProcs") else future::nbrOfWorkers()
    minBatchSize <- 1000 # put a minimum as overhead of creating processes is significant
    batchn <- max(1, min(batchn, round(length(indsNoNA) / minBatchSize)))
    
    batches <- splitInNBatches(indsNoNA, batchn)
    
    # NOTE: both the input and output is tagged with indices, which makes it much easier to see which conversions failed
    # see https://github.com/openbabel/openbabel/issues/2231
    
    mainArgs <- c("-an", paste0("-i", inFormat), paste0("-o", outFormat))
    cmdQueue <- lapply(seq_along(batches), function(bi)
    {
        b <- batches[[bi]]
        return(list(args = mainArgs, input = data.frame(input[b], b), logFile = paste0("obabel-batch_", bi, ".txt")))
    })

    results <- executeMultiProcess(cmdQueue, finishHandler = function(cmd)
    {
        # read space separated results; first column is result, second index
        ret <- fread(cmd$outFile, sep = " ")
        if (outFormat == "txt") # HACK HACK HACK: with txt (used by convertToFormulaBabel) the index is the first column
            setnames(ret, c("index", "result"))
        else
            setnames(ret, c("result", "index"))
        return(ret)
    }, prepareHandler = function(cmd)
    {
        inFile <- tempfile("obabel_in")
        fwrite(cmd$input, inFile, col.names = FALSE, sep = "\t")
        outFile <- tempfile("obabel_out")
        args <- c(cmd$args, c(inFile, "-O", outFile, "-xt", "-xw", "-e"))
        if (!is.null(extraOpts))
            args <- c(args, extraOpts)
        return(modifyList(cmd, list(command = getCommandWithOptPath("obabel", "obabel"),
                                    args = args, outFile = outFile)))
    }, showProgress = FALSE, logSubDir = "obabel")
    
    results <- rbindlist(results)
    ret <- rep(NA_character_, length(input))
    ret[results$index] <- results$result

    stopOrWarn <- function(msg) do.call(if (mustWork) stop else warning, list(msg, call. = FALSE))
    failed <- which(is.na(ret) & !is.na(input))
    if (length(failed) > 0)
    {
        for (i in failed)
            stopOrWarn(sprintf("Failed to convert %d ('%s') from %s to %s", i, input[i],
                               inFormat, outFormat))
    }
    
    return(ret)
}

convertToFormulaBabel <- function(input, inFormat, mustWork)
{
    ret <- babelConvert(input = input, inFormat = inFormat, outFormat = "txt", mustWork = mustWork,
                        extraOpts = c("--append", "formula"))
    ret <- sub("[\\+\\-]+$", "", ret) # remove trailing positive/negative charge if present
    return(ret)
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
        
        if (destFormat == "formula")
            doConv <- function(inp, f) convertToFormulaBabel(inp, f, mustWork = FALSE)
        else
            doConv <- function(inp, f) babelConvert(inp, f, destFormat, mustWork = FALSE)
        
        for (i in seq_along(fromFormats))
        {
            if (!is.null(tab[[fromCols[i]]]))
                tab[missingInTab(destCol) & !missingInTab(fromCols[i]), (destCol) := doConv(get(fromCols[i]), fromFormats[i])]
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
