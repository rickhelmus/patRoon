#' @include main.R
#' @include TP.R
NULL

# UNDONE: precursor --> parent?
# UNDONE: finalize new patRoon path options for BT and obabel (docs, verifyDependencies())

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
             c("formula", "neutralMass"))

    # No need for these...
    # NOTE: cdk:Title seems the same as "Metabolite ID" column(?)
    ret[, c("Synonyms", "PUBCHEM_CID", "cdk:Title") := NULL]

    # BUG: BT somestimes doesn't fill in the formula. Calculate them manually
    ret[!nzchar(formula), formula :=
    {
        SMI <- babelConvert(InChI, "inchi", "smi")
        mols <- getMoleculesFromSMILES(SMI)
        return(sapply(mols, function(m) rcdk::get.mol2formula(m)@string))
    }]

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
    predAll[, c("name", "precursor") := .(paste0(name, collapse = ","),
                                          paste0(precursor, collapse = ",")),
            by = "InChIKey"]

    # ... and remove now duplicates
    predAll <- unique(predAll, by = "InChIKey")

    return(predAll)
}

#' @export
predictTPsBioTransformer <- function(suspects = NULL, compounds = NULL, type = "env", steps = 2,
                                     extraOpts = NULL, skipInvalid = TRUE,
                                     logPath = file.path("log", "biotransformer"),
                                     maxProcAmount = getOption("patRoon.maxProcAmount"))
{
    # UNDONE: citations, also EnviPath
    
    if (is.null(suspects) && is.null(compounds))
        stop("Specify at least either the suspects or compounds argument.")
    
    checkmate::assertFlag(skipInvalid)

    ac <- checkmate::makeAssertCollection()
    assertSuspectList(suspects, adduct = NULL, skipInvalid = TRUE, checkAdduct = FALSE, add = ac)
    checkmate::assertClass(compounds, "compounds", null.ok = TRUE, add = ac)
    checkmate::assertChoice(type, c("ecbased", "cyp450", "phaseII", "hgut", "superbio", "allHuman", "env"), add = ac)
    checkmate::assertCount(steps, positive = TRUE, add = ac)
    checkmate::assertCharacter(extraOpts, null.ok = TRUE, add = ac)
    assertMultiProcArgs(logPath, maxProcAmount, add = ac)
    checkmate::reportAssertions(ac)

    if (!is.null(suspects))
    {
        suspects <- prepareSuspectList(suspects, adduct, skipInvalid, calcMZs = FALSE)
        if (is.null(suspects[["SMILES"]]))
            stop("No SMILES information available for suspects. Please include either SMILES or InChI columns.")
    }
    else
        suspects <- data.table()
    
    if (!is.null(compounds))
    {
        compTab <- as.data.table(compounds)
        if (!is.null(compTab[["compoundName"]]))
            compTab[, name := ifelse(nzchar(compoundName), compoundName, identifier)]
        else
            setnames(compTab, "Identifier", "name")
        suspects <- rbind(suspects, compTab[, c("name", "SMILES"), with = FALSE])
    }

    # UNDONE: move to suspect lists checks
    if (any(!nzchar(suspects$name)))
        stop("'name' column contains one or more empty values")

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
    # UNDONE: cache?

    ac <- checkmate::makeAssertCollection()
    checkmate::assertPathForOutput(out, overwrite = TRUE, add = ac) # NOTE: assert doesn't work on Windows...
    checkmate::assertFlag(includePrec, add = ac)
    checkmate::reportAssertions(ac)

    cat("Collapsing results... ")
    predAll <- collapseBTResults(pred@predictions)
    cat("Done!\n")

    # set to MetFrag style names
    setnames(predAll,
             c("name", "formula", "mass", "Precursor Major Isotope Mass"),
             c("Identifier", "MolecularFormula", "MonoisotopicMass", "Precursor MonoisotopicMass"))


    cat("Calculating SMILES... ")
    predAll[, SMILES := babelConvert(InChI, "inchi", "smi")]
    cat("Done!\n")

    if (includePrec)
    {
        cat("Adding and calculating precursor information... ")

        precs <- copy(suspects(pred))
        setnames(precs, "name", "Identifier")
        precs[, CompoundName := Identifier]

        mols <- getMoleculesFromSMILES(precs$SMILES, doTyping = TRUE, doIsotopes = TRUE)
        precs[, MolecularFormula := sapply(mols, function(m) rcdk::get.mol2formula(m)@string)]
        precs[, MonoisotopicMass := sapply(mols, rcdk::get.exact.mass)]

        precs[, InChI := babelConvert(SMILES, "smi", "inchi")]
        precs[, InChIKey := babelConvert(SMILES, "smi", "inchikey", )]

        predAll <- rbind(precs, predAll, fill = TRUE)

        cat("Done!\n")
    }

    # Add required InChIKey1 column
    predAll[, InChIKey1 := sub("\\-.*", "", InChIKey)]

    # equalize identifiers and names
    predAll[, CompoundName := Identifier]

    keepCols <- c("Identifier", "MolecularFormula", "MonoisotopicMass", "Precursor MonoisotopicMass",
                  "SMILES", "InChI", "InChIKey", "InChIKey1", "ALogP") # UNDONE: more?

    fwrite(predAll[, keepCols, with = FALSE], out)
})

setMethod("linkPrecursorsToFGroups", "TPPredictionsBT", function(pred, fGroupsScr)
{
    return(screenInfo(fGroupsScr)[name %in% names(pred), c("name", "group"), with = FALSE])
})

#' @export
setMethod("filter", "TPPredictionsBT", function(obj, removeEqualFormulas = FALSE, negate = FALSE)
{
    # UNDONE: move to base class?

    ac <- checkmate::makeAssertCollection()
    checkmate::assertFlag(removeEqualFormulas, add = ac)
    checkmate::assertFlag(negate, add = ac)
    checkmate::reportAssertions(ac)

    if (length(obj) == 0)
        return(obj)

    oldn <- length(obj)

    hash <- makeHash(obj, removeEqualFormulas, negate)
    cache <- loadCacheData("filterTPs", hash)
    if (!is.null(cache))
        obj <- cache
    else
    {
        if (removeEqualFormulas)
        {
            mols <- getMoleculesFromSMILES(obj@suspects$SMILES, doTyping = TRUE, doIsotopes = TRUE)
            pforms <- sapply(mols, function(m) rcdk::get.mol2formula(m)@string)
            obj@predictions <- mapply(pforms, obj@predictions, SIMPLIFY = FALSE, FUN = function(pform, pred)
            {
                if (negate)
                    return(pred[formula == pform])
                return(pred[formula != pform])
            })
        }

        obj@predictions <- pruneList(obj@predictions, checkZeroRows = TRUE)

        saveCacheData("filterTPs", obj, hash)
    }

    newn <- length(obj)
    printf("Done! Filtered %d (%.2f%%) TPs. Remaining: %d\n", oldn - newn, if (oldn == 0) 0 else (1-(newn/oldn))*100, newn)

    return(obj)
})
