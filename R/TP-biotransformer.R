#' @include main.R
#' @include TP.R
NULL

#' @export
TPPredictionsBT <- setClass("TPPredictionsBT", contains = "TPPredictions")

setMethod("initialize", "TPPredictionsBT",
          function(.Object, ...) callNextMethod(.Object, algorithm = "biotransformer", ...))


getBaseBTCmd <- function(precursor, SMILES, type, steps, fpType, fpSimMethod, extraOpts, baseHash)
{
    mainArgs <- c("-b", type,
                  "-k", "pred",
                  "-ismi", SMILES,
                  "-s", as.character(steps),
                  extraOpts)
    
    return(list(command = "java", args = mainArgs, logFile = paste0("biotr-", precursor, ".txt"), precursor = precursor,
                SMILES = SMILES, fpType = fpType, fpSimMethod = fpSimMethod,
                hash = makeHash(precursor, SMILES, baseHash)))
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

BTMPFinishHandler <- function(cmd)
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
    
    # BUG: BT sometimes doesn't fill in the formula. Calculate them manually
    ret[!nzchar(formula), formula :=
            {
                mols <- patRoon:::getMoleculesFromSMILES(SMILES)
                return(sapply(mols, function(m) rcdk::get.mol2formula(m)@string))
            }]
    
    # Assign some unique identifier
    ret[, name := paste0(cmd$precursor, "-TP", seq_len(nrow(ret)))]
    
    ret[, RTDir := fifelse(ALogP < `Precursor ALogP`, -1, 1)]
    
    ret[, similarity := mapply(`Precursor SMILES`, SMILES, FUN = patRoon:::distSMILES,
                               MoreArgs = list(fpType = cmd$fpType, fpSimMethod = cmd$fpSimMethod))][]
    
    return(ret)
}

BTMPPrepareHandler <- function(cmd)
{
    btBin <- path.expand(getOption("patRoon.path.BioTransformer", ""))
    if (is.null(btBin) || !nzchar(btBin) || !file.exists(btBin))
        stop("Please set the 'biotransformer' option with a (correct) path to the BioTransformer JAR file. Example: options(patRoon.path.BioTransformer = \"C:/biotransformerjar/biotransformer-2.0.3.jar\")")
    
    if (!nzchar(Sys.which("java")))
        stop("Please make sure that java is installed and its location is correctly set in PATH.")
    
    outFile <- tempfile("btresults", fileext = ".csv")
    
    cmd$args <- c("-jar", btBin, cmd$args, "-ocsv", outFile)
    cmd$outFile <- outFile
    cmd$workDir <- dirname(btBin)
    return(cmd)
}

#' @export
predictTPsBioTransformer <- function(suspects, type = "env", steps = 2, extraOpts = NULL, adduct = NULL,
                                     skipInvalid = TRUE, fpType = "extended", fpSimMethod = "tanimoto")
{
    checkmate::assert(
        checkmate::checkClass(suspects, "data.frame"),
        checkmate::checkClass(suspects, "compounds"),
        checkmate::checkClass(suspects, "featureGroupsScreening"),
        checkmate::checkClass(suspects, "featureGroupsScreeningSet"),
        .var.name = "suspects"
    )

    ac <- checkmate::makeAssertCollection()
    if (is.data.frame(suspects))
        assertSuspectList(suspects, needsAdduct = FALSE, skipInvalid = TRUE, add = ac)
    checkmate::assertChoice(type, c("ecbased", "cyp450", "phaseII", "hgut", "superbio", "allHuman", "env"), add = ac)
    checkmate::assertCount(steps, positive = TRUE, add = ac)
    checkmate::assertCharacter(extraOpts, null.ok = TRUE, add = ac)
    checkmate::assertFlag(skipInvalid, add = ac)
    aapply(checkmate::assertString, . ~ fpType + fpSimMethod, min.chars = 1, fixed = list(add = ac))
    checkmate::reportAssertions(ac)

    suspects <- getTPSuspects(suspects, adduct, skipInvalid)

    baseHash <- makeHash(type, steps, extraOpts, adduct, skipInvalid, fpType, fpSimMethod)
    setHash <- makeHash(suspects, baseHash)
    
    cmdQueue <- Map(suspects$name, suspects$SMILES, f = getBaseBTCmd,
                    MoreArgs = list(type = type, steps = steps, extraOpts = extraOpts, fpType = fpType,
                                    fpSimMethod = fpSimMethod, baseHash = baseHash))

    results <- list()

    if (length(cmdQueue) > 0)
    {
        results <- executeMultiProcess(cmdQueue, finishHandler = patRoon:::BTMPFinishHandler,
                                       prepareHandler = patRoon:::BTMPPrepareHandler,
                                       cacheName = "predictTPsBT", setHash = setHash, logSubDir = "biotransformer")
    }

    results <- pruneList(results, checkZeroRows = TRUE)
    suspects <- suspects[name %in% names(results)]

    return(TPPredictionsBT(suspects = suspects, predictions = results))
}

#' @export
setMethod("convertToMFDB", "TPPredictionsBT", function(pred, out, includePrec)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertPathForOutput(out, overwrite = TRUE, add = ac) # NOTE: assert doesn't work on Windows...
    checkmate::assertFlag(includePrec, add = ac)
    checkmate::reportAssertions(ac)

    cat("Collapsing results... ")
    predAll <- collapseBTResults(pred@predictions)
    cat("Done!\n")

    doConvertToMFDB(predAll, suspects(pred), out, includePrec)
})

setMethod("linkPrecursorsToFGroups", "TPPredictionsBT", function(pred, fGroups)
{
    return(screenInfo(fGroups)[name %in% names(pred), c("name", "group"), with = FALSE])
})

#' @export
setMethod("filter", "TPPredictionsBT", function(obj, removeEqualFormulas = FALSE, minSimilarity = NULL,
                                                negate = FALSE)
{
    # UNDONE: move to base class?

    ac <- checkmate::makeAssertCollection()
    checkmate::assertFlag(removeEqualFormulas, add = ac)
    checkmate::assertNumber(minSimilarity, lower = 0, finite = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertFlag(negate, add = ac)
    checkmate::reportAssertions(ac)

    if (length(obj) == 0)
        return(obj)

    oldn <- length(obj)

    hash <- makeHash(obj, removeEqualFormulas, minSimilarity, negate)
    cache <- loadCacheData("filterTPs", hash)
    if (!is.null(cache))
        obj <- cache
    else
    {
        if (removeEqualFormulas)
        {
            obj@predictions <- Map(suspects(obj)$formula, obj@predictions, f = function(pform, pred)
            {
                if (negate)
                    return(pred[formula == pform])
                return(pred[formula != pform])
            })
        }
        
        if (!is.null(minSimilarity))
        {
            pred <- if (negate) function(x) x < minSimilarity else function(x) numGTE(x, minSimilarity)
            obj@predictions <- lapply(obj@predictions, function(p) p[pred(similarity)])
        }

        obj@predictions <- pruneList(obj@predictions, checkZeroRows = TRUE)

        saveCacheData("filterTPs", obj, hash)
    }

    newn <- length(obj)
    printf("Done! Filtered %d (%.2f%%) TPs. Remaining: %d\n", oldn - newn, if (oldn == 0) 0 else (1-(newn/oldn))*100, newn)

    return(obj)
})
