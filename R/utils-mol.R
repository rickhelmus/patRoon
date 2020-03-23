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
        else
        {
            if (doTyping)
            {
                rcdk::do.typing(ret)
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
    # NOTE: molecules are converted to formula to get the mass instead of using
    # rcdk::get.exact.mass() so that the neutral mass can be obtained
    getMass <- function(m) rcdk::get.mol2formula(m)@mass
    
    sapply(getMoleculesFromSMILES(SMILES, doTyping = TRUE, doIsotopes = TRUE), function(m)
    {
        if (!mustWork)
        {
            if(!isValidMol(m))
                return(NA)
            
            # this could fail for some edge cases
            return(tryCatch(getMass(m), error = function(e) NA))
        }
        return(getMass(m))
    })
}

babelConvert <- function(input, inFormat, outFormat, mustWork = TRUE)
{
    # Use batch conversion with a single input/output file. Note that obabel
    # will stop after an error. This can be overidden, however, then it is
    # unclear which entries failed. Hence, this option is not used, and when an
    # error occurs the batch conversion is simply restarted with the subsequent
    # entry.
    
    inputFile <- tempfile("obabel_inp", fileext = ".txt")
    outputFile <- tempfile("obabel_out", fileext = ".txt")
    doConversion <- function(inp)
    {
        cat(inp, file = inputFile, sep = "\n")
        executeCommand(getCommandWithOptPath("obabel", "obabel"),
                       c(paste0("-i", inFormat), inputFile,
                         paste0("-o", outFormat), "-O", outputFile, "-xw"),
                       stderr = FALSE)
        # each conversion is followed by a tab (why??) and newline. Read line
        # by line and remove tab afterwards.
        ret <- readLines(outputFile)
        return(trimws(ret, which = "right", whitespace = "\t"))
    }
    
    inputLen <- length(input)
    ret <- character(inputLen)
    curIndex <- 1
    while(TRUE)
    {
        curRange <- seq(curIndex, inputLen)
        out <- doConversion(input[curRange])
        outl <- length(out)
        
        if (outl > 0)
            ret[seq(curIndex, curIndex + outl - 1)] <- out
        
        curIndex <- curIndex + outl + 1
        
        if (curIndex <= inputLen)
        {
            msg <- sprintf("Failed to convert %d ('%s')", curIndex - 1, input[curIndex - 1])
            if (mustWork)
                stop(msg)
            else
                warning(msg)
        }
        else
            break
    }
    
    return(ret)
}
