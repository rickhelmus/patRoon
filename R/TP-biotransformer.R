#' @include main.R
#' @include TP.R
NULL

# UNDONE: precursor --> parent?

#' @export
TPPredictionsBT <- setClass("TPPredictionsBT", contains = "TPPredictions")

setMethod("initialize", "TPPredictionsBT",
          function(.Object, ...) callNextMethod(.Object, algorithm = "biotransformer", ...))


getBTBin <- function()
{
    ret <- path.expand(getOption("patRoon.path.BioTransformer", ""))
    if (is.null(ret) || !nzchar(ret) || !file.exists(ret))
        stop("Please set the 'biotransformer' option with a (correct) path to the BioTransformer JAR file. Example: options(patRoon.path.BioTransformer = \"C:/biotransformerjar/biotransformer-1-0-8.jar\")")

    if (!nzchar(Sys.which("java")))
        stop("Please make sure that java is installed and its location is correctly set in PATH.")

    return(ret)
}

initBTCommand <- function(precursor, SMILES, type, steps, extraOpts, btBin, logFile)
{
    outFile <- tempfile("btresults", fileext = ".csv")

    args <- c("-b", type,
              "-k", "pred",
              "-ismi", SMILES,
              "-s", as.character(steps),
              "-ocsv", outFile,
              extraOpts)

    return(list(command = "java", args = c("-jar", btBin, args), logFile = logFile, outFile = outFile, precursor = precursor, SMILES = SMILES))
}

processBTResults <- function(cmd)
{
    if (!file.exists(cmd$outFile))
        return(data.table()) # no results

    ret <- fread(cmd$outFile, colClasses = c("Precursor ID" = "character", Synonyms = "character",
                                             "Molecular formula" = "character"))

    # UNDONE: transform column names, more?

    # Simplify/harmonize columns a bit
    setnames(ret,
             c("Molecular formula", "Major Isotope Mass"),
             c("formula", "mass"))
    
    # No need for these...
    # NOTE: cdk:Title seems the same as "Metabolite ID" column(?)
    ret[, c("Synonyms", "PUBCHEM_CID", "cdk:Title") := NULL]
    
    # BUG: BT somestimes doesn't fill in the formula. Calculate them manually
    ret[!nzchar(formula),
        formula := sapply(InChI, function(i) rcdk::get.mol2formula(rinchi::parse.inchi(i)[[1]])@string)]
    
    # Assign some unique identifier
    ret[, name := paste0(cmd$precursor, "-TP", seq_len(nrow(ret)))]
    
    ret[, RTDir := ifelse(ALogP < `Precursor ALogP`, -1, 1)]
    
    return(ret)
}

collapseBTResults <- function(pred)
{
    pred <- lapply(pred, function(p)
    {
        # merge duplicate compound rows, which can occur due to consecutive
        # reactions giving the same TP
        p <- copy(p)
        col <- "Precursor ID"
        p[!nzchar(get(col)), (col) := "parent"]
        p[, (col) := paste0(get(col), collapse = "/"), by = "InChIKey"]
        p <- unique(p, by = "InChIKey")
    })
    
    predAll <- rbindlist(pred, idcol = "precursor")
    
    # merge precursor and sub-precursor (ie from consecutive reactions)
    predAll[, precursor := paste0(precursor, " (", `Precursor ID`, ")")]

    # combine equal TPs from different precursors
    predAll[, precursor := paste0(precursor, collapse = ","), by = "InChIKey"]
    # ... and remove now duplicates
    predAll <- unique(predAll, by = "InChIKey")

    return(predAll)    
}

getPrecursorSuspList <- function(pred, adduct)
{
    addMZ <- adductMZDelta(adduct)
    ret <- copy(suspects(pred))
    ret[, mz := sapply(getMoleculesFromSMILES(SMILES, doTyping = TRUE, doIsotopes = TRUE), rcdk::get.exact.mass) + addMZ]
}

#' @export
predictTPsBioTransformer <- function(suspects, type = "env", steps = 2, extraOpts = NULL,
                                     logPath = file.path("log", "biotransformer"),
                                     maxProcAmount = getOption("patRoon.maxProcAmount"))
{
    # UNDONE: as long as BT may return empty formulas we need this
    checkPackage("rinchi", "CDK-R/rinchi")
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertDataFrame(suspects, any.missing = FALSE, min.rows = 1, add = ac)
    assertHasNames(suspects, c("name", "SMILES"), add = ac)
    checkmate::assertChoice(type, c("ecbased", "cyp450", "phaseII", "hgut", "superbio", "allHuman", "env"), add = ac)
    checkmate::assertCount(steps, positive = TRUE, add = ac)
    checkmate::assertCharacter(extraOpts, null.ok = TRUE, add = ac)
    assertMultiProcArgs(logPath, maxProcAmount, add = ac)
    checkmate::reportAssertions(ac)

    suspects <- as.data.table(suspects)

    if (any(!nzchar(suspects$name)))
        stop("'name' column contains one ore more empty values")
    if (any(!nzchar(suspects$SMILES)))
        stop("'SMILES' column contains one ore more empty values")
    
    cacheDB <- openCacheDBScope()
    baseHash <- makeHash(type, steps, extraOpts)
    setHash <- makeHash(suspects, baseHash)
    hashes <- sapply(suspects$SMILES, makeHash, baseHash)
    cachedSet <- loadCacheSet("predictTPsBT", setHash, cacheDB)

    btBin <- getBTBin()
    btPath <- dirname(btBin) # BT has to be executed from its own directory

    if (!is.null(logPath))
    {
        mkdirp(logPath)
        logPath <- normalizePath(logPath) # need full path since we will temporarily change the working directory.
    }

    cachedResults <- sapply(hashes, function(hash)
    {
        ret <- NULL
        if (!is.null(cachedSet))
            ret <- cachedSet[[hash]]
        if (is.null(ret))
            ret <- loadCacheData("predictTPsBT", hash, cacheDB)
        return(ret)
    }, simplify = FALSE)
    names(cachedResults) <- suspects$name
    cachedResults <- pruneList(cachedResults)

    doSuspects <- suspects[!name %in% names(cachedResults)]
    cmdQueue <- mapply(doSuspects$name, doSuspects$SMILES, SIMPLIFY = FALSE, FUN = function(n, sm)
    {
        logf <- if (!is.null(logPath)) file.path(logPath, paste0("biotr-", n, ".txt")) else NULL
        initBTCommand(n, sm, type = type, steps = steps, extraOpts = extraOpts, btBin = btBin, logFile = logf)
    })
    names(cmdQueue) <- doSuspects$name

    results <- list()

    if (length(cmdQueue) > 0)
    {
        results <- executeMultiProcess(cmdQueue, function(cmd)
        {
            res <- processBTResults(cmd)
            saveCacheData("predictTPsBT", res, hashes[[cmd$SMILES]], cacheDB)
            return(res)
        }, workDir = btPath, maxProcAmount = maxProcAmount)

        if (is.null(cachedSet))
            saveCacheSet("predictTPsBT", hashes, setHash, cacheDB)
    }

    if (length(cachedResults) > 0)
    {
        results <- c(results, cachedResults)
        results <- results[intersect(suspects$name, names(results))] # re-order
    }

    results <- pruneList(results, checkZeroRows = TRUE)
    suspects <- suspects[name %in% names(results)]
    
    return(TPPredictionsBT(suspects = suspects, predictions = results))
}

#' @export
setMethod("convertToMFDB", "TPPredictionsBT", function(pred, out, includePrec)
{
    # rinchi: remotes::install_github("CDK-R/rinchi", INSTALL_opts = "--no-multiarch")
    checkPackage("rinchi", "CDK-R/rinchi")
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertPathForOutput(out, overwrite = TRUE, add = ac) # NOTE: assert doesn't work on Windows...
    checkmate::assertFlag(includePrec, add = ac)
    checkmate::reportAssertions(ac)
    
    predAll <- collapseBTResults(pred@predictions)

    # set to MetFrag style names
    setnames(predAll,
             c("name", "formula", "mass", "Precursor Major Isotope Mass"),
             c("Identifier", "MolecularFormula", "MonoisotopicMass", "Precursor MonoisotopicMass"))

    mols <- suppressWarnings(rinchi::parse.inchi(predAll$InChI))
    nullMols <- which(sapply(mols, is.null))
    if (length(nullMols) > 0)
    {
        warning(paste("Failed to convert some InChI strings:\n",
                      paste0(sprintf("%s (%s)\n", predAll$InChI[nullMols], predAll$name[nullMols]),
                             collapse = "\n")))
        predAll <- predAll[-nullMols]
        mols <- mols[-nullMols]
    }
    
    predAll[, SMILES := sapply(mols, rcdk::get.smiles)]
        
    if (includePrec)
    {
        precs <- copy(suspects(pred))
        setnames(precs, "name", "Identifier")
        precs[, CompoundName := Identifier]
        
        mols <- getMoleculesFromSMILES(precs$SMILES, doTyping = TRUE, doIsotopes = TRUE)
        precs[, MolecularFormula := sapply(mols, function(m) rcdk::get.mol2formula(m)@string)]
        precs[, MonoisotopicMass := sapply(mols, rcdk::get.exact.mass)]
        
        precs[, InChI := sapply(SMILES, rinchi::get.inchi)]
        precs[, InChIKey := sapply(SMILES, rinchi::get.inchi.key)]
        
        predAll <- rbind(precs, predAll, fill = TRUE)
    }
    
    # Add required InChIKey1 column
    predAll[, InChIKey1 := sub("\\-.*", "", InChIKey)]

    # equalize identifiers and name's
    predAll[, CompoundName := Identifier]
    
    fwrite(predAll, out)
})

setMethod("linkPrecursorsToFGroups", "TPPredictionsBT", function(pred, fGroups, adduct, mzWindow)
{
    suspList <- getPrecursorSuspList(pred, adduct)
    scr <- screenSuspects(fGroups, suspList, mzWindow = mzWindow)
    scr <- scr[!is.na(group)] # remove non-hits
    return(scr[, c("name", "group")])
})
