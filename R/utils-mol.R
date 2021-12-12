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
    # Use batch conversion with a single input/output file. Note that obabel
    # will stop after an error. This can be overidden, however, then it is
    # unclear which entries failed. Hence, this option is not used, and when an
    # error occurs the batch conversion is simply restarted with the subsequent
    # entry.
    
    if (length(input) == 0)
        return(character())
    
    inputFile <- tempfile("obabel_inp", fileext = ".txt")
    outputFile <- tempfile("obabel_out", fileext = ".txt")
    doConversion <- function(inp)
    {
        cat(inp, file = inputFile, sep = "\n")
        
        args <- c(paste0("-i", inFormat), inputFile,
                  paste0("-o", outFormat), "-O", outputFile, "-xw")
        if (!is.null(extraOpts))
            args <- c(args, extraOpts)
        
        executeCommand(getCommandWithOptPath("obabel", "obabel"), args, stderr = FALSE)
        # each conversion is followed by a tab (why??) and newline. Read line
        # by line and remove tab afterwards.
        ret <- readLines(outputFile)
        return(trimws(ret, which = "right", whitespace = "\t"))
    }
    
    stopOrWarn <- function(msg)
    {
        do.call(if (mustWork) stop else warning, list(msg, call. = FALSE))
    }

    inpNA <- is.na(input)
    doInputs <- input[!inpNA]
    inputLen <- length(doInputs)
    conv <- character(inputLen)
    curIndex <- 1
    while(TRUE)
    {
        curRange <- seq(curIndex, inputLen)
        out <- doConversion(doInputs[curRange])
        outl <- length(out)
        
        if (outl > 0)
            conv[seq(curIndex, curIndex + outl - 1)] <- out
        else # silent fail
            conv[curIndex] <- NA_character_
        
        curIndex <- curIndex + outl + 1
        
        if (curIndex <= inputLen)
            stopOrWarn(sprintf("Failed to convert %d ('%s') from %s to %s", curIndex - 1, doInputs[curIndex - 1],
                               inFormat, outFormat))
        else
            break
    }
    
    # merge back in NA inputs
    ret <- character(length(input))
    ret[inpNA] <- NA_character_
    ret[!inpNA] <- conv
    
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
