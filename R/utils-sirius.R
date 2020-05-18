#' @include main.R

getSiriusBin <- function()
{
    si <- Sys.info()

    if (si[["sysname"]] == "Linux")
        return("sirius")

    # windows
    if (si[["machine"]] == "x86-64")
        return("sirius-console-64")

    return("sirius-console-32")
}

isSIRIUSPre44 <- function()
{
    # SIRIUS 4.4 returns a version string when running with --version, older
    # versions don't actually report version info...
    return(!any(grepl("^SIRIUS 4\\.", executeCommand(getCommandWithOptPath(getSiriusBin(), "SIRIUS"),
                                                     "--version", stdout = TRUE, stderr = FALSE))))
}

getSiriusResultPath <- function(outPath, msFName, cmpName, isPre44)
{
    # format is resultno_specname_compoundname, older versions start with 1, newer with 0
    msFName <- basename(tools::file_path_sans_ext(msFName))
    return(list.files(outPath, pattern = sprintf("[0-9]+_%s_%s", msFName, cmpName), full.names = TRUE))
}

getSiriusFragFiles <- function(resultPath, isPre44)
{
    if (isPre44)
        pat <- "[:0-9:]+_([A-Za-z0-9]+).*\\.ms"
    else
        pat <- "([A-Za-z0-9]+).*\\.csv"
    return(list.files(file.path(resultPath, "spectra"), full.names = TRUE, pattern = pat))
}

getFormulaFromSiriusFragFile <- function(ffile, isPre44)
{
    if (isPre44)
        pat <- "[:0-9:]+_([A-Za-z0-9]+).*\\.ms"
    else
        pat <- "([A-Za-z0-9]+).*\\.csv"
    return(gsub(pat, "\\1", basename(ffile)))
}

makeSirMSFile <- function(plistMS, plistMSMS, parentMZ, compound, ionization, out)
{
    msFile <- file(out, "w")

    writeMeta <- function(var, data) cat(sprintf(">%s %s\n", var, data), file = msFile)

    writeMeta("compound", compound)
    writeMeta("parentmass", parentMZ)
    writeMeta("ionization", ionization)

    cat(">ms1peaks\n", file = msFile)
    write.table(plistMS[, c("mz", "intensity")], msFile, row.names = FALSE, col.names = FALSE)

    cat("\n>ms2peaks\n", file = msFile)
    write.table(plistMSMS[, c("mz", "intensity")], msFile, row.names = FALSE, col.names = FALSE)

    close(msFile)
}

unifySirNames <- function(sir)
{
    unNames <- c(# MonoisotopicMass = "neutralMass", UNDONE
                 smiles = "SMILES",
                 inchikey2D = "InChIKey1",
                 inchi = "InChI",
                 pubchemids = "identifier",
                 PubChemNumberPatents = "numberPatents",
                 score = "score",
                 molecularFormula = "formula",
                 xlogp = "XlogP",
                 name = "compoundName",
                 links = "libraryLinks",
                 
                 # some names were changed in 4.4 and new columns were added
                 # UNDONE: there is also a compound_dentifications.csv file with slightly different columns, use that?
                 formulaRank = "formulaRank",
                 InChI = "InChI",
                 InChIkey2D = "InChIKey1",
                 "CSI:FingerID_Score" = "score",
                 TreeIsotope_Score = "SIR_formulaScore" # UNDONE: better name?
                 )

    unNames <- unNames[names(unNames) %in% names(sir)] # filter out missing
    setnames(sir, names(unNames), unNames)

    return(sir[, unNames, with = FALSE]) # filter out any other columns
}

# get a command queue list that can be used with executeMultiProcess()
getSiriusCommand <- function(precursorMZ, MSPList, MSMSPList, profile, adduct, ppmMax, elements,
                             database, noise, withFingerID, fingerIDDatabase, topMost, extraOpts,
                             isPre44)
{
    outPath <- tempfile("sirius")
    # unlink(outPath, TRUE) # start with fresh output directory (otherwise previous results are combined)

    stopifnot(!file.exists(outPath))

    msFName <- tempfile("spec", fileext = ".ms")    
    ionization <- as.character(adduct, format = "sirius")
    mainArgs <- c("-p", profile,
                  "-e", elements,
                  "--ppm-max", ppmMax,
                  "-c", topMost)

    if (!is.null(database))
        mainArgs <- c(mainArgs, "-d", database)
    if (!is.null(noise))
        mainArgs <- c(mainArgs, "-n", noise)
    if (!is.null(extraOpts))
        mainArgs <- c(mainArgs, extraOpts)

    if (isPre44)
    {
        if (withFingerID)
            mainArgs <- c(mainArgs, "--fingerid", "--fingerid-db", fingerIDDatabase)
        args <- c(mainArgs, "-o", outPath, msFName)
    }
    else
    {
        args <- c("-o", outPath, "-i", msFName, "formula", mainArgs)
        if (withFingerID)
            args <- c(args, "structure", "--database", fingerIDDatabase)
    }

    cmpName <- "unknownCompound"
    makeSirMSFile(MSPList, MSMSPList, precursorMZ, cmpName, ionization, msFName)

    return(list(command = getCommandWithOptPath(getSiriusBin(), "SIRIUS"), args = args,
                outPath = outPath, msFName = msFName, cmpName = cmpName, isPre44 = isPre44))
}

# get a command queue list that can be used with executeMultiProcess()
getSiriusCommandBatch <- function(precursorMZs, MSPLists, MSMSPLists, profile, adduct, ppmMax, elements,
                                  database, noise, cores, withFingerID, fingerIDDatabase, topMost, extraOpts,
                                  isPre44)
{
    inPath <- tempfile("sirius_in")
    outPath <- tempfile("sirius_out")
    # unlink(outPath, TRUE) # start with fresh output directory (otherwise previous results are combined)
    stopifnot(!file.exists(outPath) && !file.exists(outPath))
    dir.create(inPath)
    
    ionization <- as.character(adduct, format = "sirius")
    cmpName <- "unknownCompound"
    
    msFNames <- mapply(precursorMZs, MSPLists, MSMSPLists, FUN = function(pmz, mspl, msmspl)
    {
        ret <- tempfile("spec", fileext = ".ms", tmpdir = inPath)
        makeSirMSFile(mspl, msmspl, pmz, cmpName, ionization, ret)
        return(ret)
    })

    mainArgs <- character()
    if (!is.null(cores))
        mainArgs <- c("--cores", cores)
    
    formArgs <- c("-p", profile,
                  "-e", elements,
                  "--ppm-max", ppmMax,
                  "-c", topMost)
    
    if (!is.null(database))
        formArgs <- c(formArgs, "-d", database)
    if (!is.null(noise))
        formArgs <- c(formArgs, "-n", noise)
    if (!is.null(extraOpts))
        formArgs <- c(formArgs, extraOpts)
    
    if (isPre44)
    {
        if (withFingerID)
            formArgs <- c(formArgs, "--fingerid", "--fingerid-db", fingerIDDatabase)
        args <- c(mainArgs, formArgs, "-o", outPath, inPath)
    }
    else
    {
        args <- c(mainArgs, "-o", outPath, "-i", inPath, "formula", formArgs)
        if (withFingerID)
            args <- c(args, "structure", "--database", fingerIDDatabase)
    }
    
    return(list(command = getCommandWithOptPath(getSiriusBin(), "SIRIUS"), args = args,
                outPath = outPath, msFNames = msFNames, cmpName = cmpName, isPre44 = isPre44))
}
