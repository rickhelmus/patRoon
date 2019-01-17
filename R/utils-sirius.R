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

getSiriusFragFiles <- function(resultPath)
{
    pat <- "[:0-9:]+_([A-Za-z0-9]+).*\\.ms"
    return(list.files(file.path(resultPath, "spectra"), full.names = TRUE, pattern = pat))
}

getFormulaFromSiriusFragFile <- function(ffile)
{
    pat <- "[:0-9:]+_([A-Za-z0-9]+).*\\.ms"
    gsub(pat, "\\1", basename(ffile))
}

makeSirMSFile <- function(plistMS, plistMSMS, parentMZ, compound, ionization, out)
{
    msFile <- file(out, "w")

    writeMeta <- function(var, data) cat(sprintf(">%s %s\n", var, data), file = msFile)

    writeMeta("compound", compound)
    writeMeta("parentmass", parentMZ)
    writeMeta("ionization", ionization)

    cat(">ms1peaks\n", file = msFile)
    write.table(plistMS, msFile, row.names = FALSE, col.names = FALSE)

    cat("\n>ms2peaks\n", file = msFile)
    write.table(plistMSMS, msFile, row.names = FALSE, col.names = FALSE)

    close(msFile)
}

unifySirNames <- function(sir)
{
    unNames <- c(# MonoisotopicMass = "neutralMass", UNDONE
                 smiles = "SMILES",
                 inchikey2D = "InChIKey1",
                 inchi = "InChi",
                 pubchemids = "identifier",
                 PubChemNumberPatents = "numberPatents",
                 score = "score",
                 molecularFormula = "formula",
                 xlogp = "XlogP",
                 name = "compoundName",
                 links = "libraryLinks")

    unNames <- unNames[names(unNames) %in% names(sir)] # filter out missing
    setnames(sir, names(unNames), unNames)

    return(sir[, unNames, with = FALSE]) # filter out any other columns
}

# get a command queue list that can be used with executeMultiProcess()
getSiriusCommand <- function(precursorMZ, MSPList, MSMSPList, profile, adduct, ppmMax, elements,
                             database, noise, withFingerID, fingerIDDatabase, topMost, extraOpts)
{
    outPath <- tempfile("sirius")
    # unlink(outPath, TRUE) # start with fresh output directory (otherwise previous results are combined)

    stopifnot(!file.exists(outPath))

    ionization <- as.character(adduct, format = "sirius")
    mainArgs <- c("-p", profile,
                  "-i", ionization,
                  "-e", elements,
                  "--ppm-max", ppmMax,
                  "-c", topMost,
                  "-o", outPath)

    if (!is.null(database))
        mainArgs <- c(mainArgs, "-d", database)
    if (!is.null(noise))
        mainArgs <- c(mainArgs, "-n", noise)
    if (withFingerID)
        mainArgs <- c(mainArgs, "--fingerid")
    if (!is.null(fingerIDDatabase))
        mainArgs <- c(mainArgs, "--fingerid-db", fingerIDDatabase)
    if (!is.null(extraOpts))
        mainArgs <- c(mainArgs, extraOpts)

    msFName <- tempfile("spec", fileext = ".ms")
    cmpName <- "unknownCompound"
    makeSirMSFile(MSPList, MSMSPList, precursorMZ, cmpName, ionization, msFName)

    return(list(command = getCommandWithOptPath(getSiriusBin(), "SIRIUS"), args = c(mainArgs, msFName),
                outPath = outPath, msFName = msFName, cmpName = cmpName))
}
