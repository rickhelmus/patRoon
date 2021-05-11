#' @include main.R
#' @include TP.R
NULL

#' @export
transformationProductsBT <- setClass("transformationProductsBT", contains = "transformationProducts")

setMethod("initialize", "transformationProductsBT",
          function(.Object, ...) callNextMethod(.Object, algorithm = "biotransformer", ...))


getBaseBTCmd <- function(parent, SMILES, type, steps, fpType, fpSimMethod, extraOpts, baseHash)
{
    mainArgs <- c("-b", type,
                  "-k", "pred",
                  "-ismi", SMILES,
                  "-s", as.character(steps),
                  extraOpts)
    
    return(list(command = "java", args = mainArgs, logFile = paste0("biotr-", parent, ".txt"), parent = parent,
                SMILES = SMILES, fpType = fpType, fpSimMethod = fpSimMethod,
                hash = makeHash(parent, SMILES, baseHash)))
}

collapseBTResults <- function(prod)
{
    prod <- lapply(prod, function(p)
    {
        # merge duplicate compound rows, which can occur due to consecutive
        # transformation giving the same TP
        p <- copy(p)
        col <- "parent_ID"
        p[!nzchar(get(col)), (col) := "parent"]
        p[, (col) := paste0(get(col), collapse = "/"), by = "InChIKey"]
        p <- unique(p, by = "InChIKey")
    })

    prodAll <- rbindlist(prod, idcol = "parent")

    # merge parent and sub-parent (ie from consecutive transformation)
    prodAll[, parent := paste0(parent, " (", parent_ID, ")")]

    # combine equal TPs from different parents
    prodAll[, c("name", "parent") := .(paste0(name, collapse = ","), paste0(parent, collapse = ",")), by = "InChIKey"]

    # ... and remove now duplicates
    prodAll <- unique(prodAll, by = "InChIKey")

    return(prodAll)
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
    setnames(ret, sub("^Precursor ", "parent_", names(ret)))
    setnames(ret, c("Reaction", "Reaction ID"), c("transformation", "transformation_ID"))
    
    # No need for these...
    # NOTE: cdk:Title seems the same as "Metabolite ID" column(?)
    ret[, c("Synonyms", "PUBCHEM_CID", "cdk:Title") := NULL]
    
    # BUG: BT sometimes doesn't fill in the formula. Calculate them manually
    ret[!nzchar(formula), formula := {
        mols <- patRoon:::getMoleculesFromSMILES(SMILES)
        return(sapply(mols, function(m) rcdk::get.mol2formula(m)@string))
    }]
    
    # Assign some unique identifier
    ret[, name := paste0(cmd$parent, "-TP", seq_len(nrow(ret)))]
    
    ret[, retDir := fifelse(ALogP < parent_ALogP, -1, 1)]
    
    ret[, similarity := mapply(parent_SMILES, SMILES, FUN = patRoon:::distSMILES,
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
generateTPsBioTransformer <- function(parents, type = "env", steps = 2, extraOpts = NULL, adduct = NULL,
                                      skipInvalid = TRUE, fpType = "extended", fpSimMethod = "tanimoto")
{
    checkmate::assert(
        checkmate::checkClass(parents, "data.frame"),
        checkmate::checkClass(parents, "compounds"),
        checkmate::checkClass(parents, "featureGroupsScreening"),
        checkmate::checkClass(parents, "featureGroupsScreeningSet"),
        .var.name = "parents"
    )

    ac <- checkmate::makeAssertCollection()
    if (is.data.frame(parents))
        assertSuspectList(parents, needsAdduct = FALSE, skipInvalid = TRUE, add = ac)
    checkmate::assertChoice(type, c("ecbased", "cyp450", "phaseII", "hgut", "superbio", "allHuman", "env"), add = ac)
    checkmate::assertCount(steps, positive = TRUE, add = ac)
    checkmate::assertCharacter(extraOpts, null.ok = TRUE, add = ac)
    checkmate::assertFlag(skipInvalid, add = ac)
    aapply(checkmate::assertString, . ~ fpType + fpSimMethod, min.chars = 1, fixed = list(add = ac))
    checkmate::reportAssertions(ac)

    parents <- getTPParents(parents, adduct, skipInvalid)

    baseHash <- makeHash(type, steps, extraOpts, adduct, skipInvalid, fpType, fpSimMethod)
    setHash <- makeHash(parents, baseHash)
    
    cmdQueue <- Map(parents$name, parents$SMILES, f = getBaseBTCmd,
                    MoreArgs = list(type = type, steps = steps, extraOpts = extraOpts, fpType = fpType,
                                    fpSimMethod = fpSimMethod, baseHash = baseHash))

    results <- list()

    if (length(cmdQueue) > 0)
    {
        results <- executeMultiProcess(cmdQueue, finishHandler = patRoon:::BTMPFinishHandler,
                                       prepareHandler = patRoon:::BTMPPrepareHandler,
                                       cacheName = "generateTPsBT", setHash = setHash, logSubDir = "biotransformer")
    }

    results <- pruneList(results, checkZeroRows = TRUE)
    parents <- parents[name %in% names(results)]

    return(transformationProductsBT(parents = parents, products = results))
}

#' @export
setMethod("convertToMFDB", "transformationProductsBT", function(TPs, out, includeParents)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertPathForOutput(out, overwrite = TRUE, add = ac) # NOTE: assert doesn't work on Windows...
    checkmate::assertFlag(includeParents, add = ac)
    checkmate::reportAssertions(ac)

    cat("Collapsing results... ")
    prodAll <- collapseBTResults(TPs@products)
    cat("Done!\n")

    doConvertToMFDB(prodAll, parents(TPs), out, includeParents)
})

setMethod("linkParentsToFGroups", "transformationProductsBT", function(TPs, fGroups)
{
    return(screenInfo(fGroups)[name %in% names(TPs), c("name", "group"), with = FALSE])
})

#' @export
setMethod("filter", "transformationProductsBT", function(obj, removeEqualFormulas = FALSE, minSimilarity = NULL,
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
            # NOTE: obj@products should be first arg to Map to keep names...
            obj@products <- Map(obj@products, parents(obj)$formula, f = function(prod, pform)
            {
                if (negate)
                    return(prod[formula == pform])
                return(prod[formula != pform])
            })
        }
        
        if (!is.null(minSimilarity))
        {
            pred <- if (negate) function(x) x < minSimilarity else function(x) numGTE(x, minSimilarity)
            obj@products <- lapply(obj@products, function(p) p[pred(similarity)])
        }

        obj@products <- pruneList(obj@products, checkZeroRows = TRUE)
        obj@parents <- obj@parents[name %in% names(obj@products)]

        saveCacheData("filterTPs", obj, hash)
    }

    newn <- length(obj)
    printf("Done! Filtered %d (%.2f%%) TPs. Remaining: %d\n", oldn - newn, if (oldn == 0) 0 else (1-(newn/oldn))*100, newn)

    return(obj)
})
