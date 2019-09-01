#' @include main.R
NULL

getBTBin <- function()
{
    ret <- path.expand(getOption("patRoon.path.BioTransformer", ""))
    if (is.null(ret) || !nzchar(ret) || !file.exists(ret))
        stop("Please set the 'biotransformer' option with a (correct) path to the BioTransformer JAR file. Example: options(patRoon.path.BioTransformer = \"C:/biotransformerjar/biotransformer-1-0-8.jar\")")

    if (!nzchar(Sys.which("java")))
        stop("Please make sure that java is installed and its location is correctly set in PATH.")

    return(ret)
}

initBTCommand <- function(SMILES, type, steps, extraOpts, btBin, logFile)
{
    outFile <- tempfile("btresults", fileext = ".csv")

    args <- c("-b", type,
              "-k", "pred",
              "-ismi", SMILES,
              "-s", as.character(steps),
              "-ocsv", outFile,
              extraOpts)

    return(list(command = "java", args = c("-jar", btBin, args), logFile = logFile, outFile = outFile, SMILES = SMILES))
}

processBTResults <- function(cmd)
{
    if (!file.exists(cmd$outFile))
        return(data.table()) # no results

    ret <- fread(cmd$outFile, colClasses = c("Precursor ID" = "character"))

    # UNDONE: transform column names, more?

    return(ret)
}

predictTPsBioTransformer <- function(suspects, type = "env", steps = 2, extraOpts = NULL,
                                     logPath = file.path("log", "biotransformer"),
                                     maxProcAmount = getOption("patRoon.maxProcAmount"))
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertDataFrame(suspects, any.missing = FALSE, min.rows = 1, add = ac)
    assertHasNames(suspects, c("name", "SMILES"), add = ac)
    checkmate::assertChoice(type, c("ecbased", "cyp450", "phaseII", "hgut", "superbio", "allHuman", "env"), add = ac)
    checkmate::assertCount(steps, positive = TRUE, add = ac)
    checkmate::assertCharacter(extraOpts, null.ok = TRUE, add = ac)
    assertMultiProcArgs(logPath, maxProcAmount, add = ac)
    checkmate::reportAssertions(ac)

    suspects <- as.data.table(suspects)

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
        initBTCommand(sm, type = type, steps = steps, extraOpts = extraOpts, btBin = btBin, logFile = logf)
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

    # UNDONE: make nice S4 class?

    return(results)
}

convertTPPredictionToMFDB <- function(pred, out)
{
    pred <- pruneList(pred, checkZeroRows = TRUE)

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

    # set to MetFrag style names, also ensure that minimally required columns are present
    setnames(predAll,
             c("Synonyms", "Molecular formula", "Major Isotope Mass", "Precursor Major Isotope Mass"),
             c("CompoundName", "MolecularFormula", "MonoisotopicMass", "Precursor MonoisotopicMass"))

    # Add required InChIKey1 column
    predAll[, InChIKey1 := sub("\\-.*", "", InChIKey)]

    # seems the same as "Metabolite ID" column(?)
    predAll[, `cdk:Title` := NULL]

    # merge precursor and sub-precursor (ie from consecutive reactions)
    predAll[, precursor := paste0(precursor, " (", `Precursor ID`, ")")]
    predAll[, `Precursor ID` := NULL]

    # combine equal TPs from different precursors
    predAll[, precursor := paste0(precursor, collapse = ","), by = "InChIKey"]
    # ... and remove duplicates to not bother MetFrag with them
    predAll <- unique(predAll, by = "InChIKey")

    # Assign some unique identifier
    predAll[, Identifier := paste0("TP", seq_len(nrow(predAll)))]

    fwrite(predAll, out)
}
