#' @include main.R
#' @include TP.R
NULL

#' Class to store transformation products (TPs) predicted by BioTransformer
#'
#' This class is used to store prediction results that are generated with
#' \href{http://biotransformer.ca/}{BioTransformer}.
#'
#' Objects from this class are generate with \code{\link{generateTPsBioTransformer}}. This class is derived from the
#' \code{\link{transformationProducts}} base class, please see its documentation for more details.
#'
#' @param obj,TPs \code{transformationProductsBTs} object to be accessed
#'
#' @seealso The base class \code{\link{transformationProducts}} for more relevant methods and
#'   \code{\link{generateTPs}}
#'
#' @references \insertRef{DjoumbouFeunang2019}{patRoon} \cr\cr \insertRef{Wicker2015}{patRoon}
#'
#' @templateVar class transformationProductsBT
#' @template class-hierarchy
#'
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

#' Obtain transformation products (TPs) with BioTransformer
#'
#' Uses \href{http://biotransformer.ca/}{BioTransformer} to predict TPs
#'
#' @templateVar algo BioTransformer
#' @templateVar do obtain transformation products
#' @templateVar generic generateTPs
#' @templateVar algoParam biotransformer
#' @template algo_generator
#'
#' @details Structural similarities between the parent and its TPs are calculated, which can be used to
#'   \link[=filter,transformationProductsBT-method]{filter} the results.
#'
#'   In order to use this function the \file{.jar} command line utility should be installed and specified in the
#'   \code{\link[=patRoon-package]{patRoon.path.BioTransformer}} option. The \file{.jar} file can be obtained via
#'   \url{https://bitbucket.org/djoumbou/biotransformer/src/master}.
#'
#' @param type The type of prediction. Valid values are: \code{"env"}, \code{"ecbased"}, \code{"cyp450"},
#'   \code{"phaseII"}, \code{"hgut"}, \code{"superbio"}, \code{"allHuman"}. Sets the \command{-b} command line option.
#' @param steps The number of steps for the predictions. Sets the \command{-s} command line option.
#' @param extraOpts A \code{character} with extra command line options passed to the \command{biotransformer.jar} tool.
#' @param MP If \code{TRUE} then multiprocessing is enabled. Since \command{BioTransformer} supports native
#'   parallelization, additional multiprocessing generally doesn't lead to significant reduction in computational times.
#'   Furthermore, enabling multiprocessing can lead to very high CPU/RAM usage.
#'
#' @return The TPs are stored in an object from the \code{\link{transformationProductsBT}} class.
#' 
#' @template tp_gen-scr
#' @template fp-args
#'
#' @templateVar what \code{generateTPsBioTransformer}
#' @template uses-multiProc
#'
#' @references \insertRef{DjoumbouFeunang2019}{patRoon} \cr\cr \insertRef{Wicker2015}{patRoon} \cr\cr
#'   \addCitations{rcdk}{1}
#'
#' @export
generateTPsBioTransformer <- function(parents, type = "env", steps = 2, extraOpts = NULL,
                                      skipInvalid = TRUE, fpType = "extended", fpSimMethod = "tanimoto", MP = FALSE)
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
    checkmate::assertFlag(MP, add = ac)
    checkmate::reportAssertions(ac)

    parents <- getTPParents(parents, skipInvalid)

    baseHash <- makeHash(type, steps, extraOpts, skipInvalid, fpType, fpSimMethod)
    setHash <- makeHash(parents, baseHash)
    
    cmdQueue <- Map(parents$name, parents$SMILES, f = getBaseBTCmd,
                    MoreArgs = list(type = type, steps = steps, extraOpts = extraOpts, fpType = fpType,
                                    fpSimMethod = fpSimMethod, baseHash = baseHash))

    results <- list()

    if (length(cmdQueue) > 0)
    {
        if (!MP)
            withr::local_options(list(patRoon.MP.maxProcs = 1))
        results <- executeMultiProcess(cmdQueue, finishHandler = patRoon:::BTMPFinishHandler,
                                       prepareHandler = patRoon:::BTMPPrepareHandler,
                                       cacheName = "generateTPsBT", setHash = setHash, logSubDir = "biotransformer")
    }

    results <- pruneList(results, checkZeroRows = TRUE)
    parents <- parents[name %in% names(results)]

    return(transformationProductsBT(parents = parents, products = results))
}

#' @templateVar class transformationProductsBT
#' @template convertToMFDB
#' @export
setMethod("convertToMFDB", "transformationProductsBT", function(TPs, out, includeParents = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertPathForOutput(out, overwrite = TRUE, add = ac) # NOTE: assert doesn't work on Windows...
    checkmate::assertFlag(includeParents, add = ac)
    checkmate::reportAssertions(ac)

    cat("Collapsing results... ")
    prodAll <- if (length(TPs) > 0) collapseBTResults(TPs@products) else data.table()
    cat("Done!\n")

    doConvertToMFDB(prodAll, parents(TPs), out, includeParents)
})

setMethod("linkParentsToFGroups", "transformationProductsBT", function(TPs, fGroups)
{
    return(screenInfo(fGroups)[name %in% names(TPs), c("name", "group"), with = FALSE])
})

#' @describeIn transformationProductsBT Performs rule-based filtering of the \command{BioTransformer} predictions.
#'   Useful to simplify and clean-up the data.
#'
#' @param removeDuplicates If \code{TRUE} then the TPs of a parent with duplicate structures (\acronym{SMILES}) are
#'   removed. Such duplicates may occur when different transformation pathways yield the same TPs. The first TP
#'   candidate with duplicate structure will be kept.
#' @param removeParentIsomers If \code{TRUE} then TPs with an equal formula as their parent (isomers) are removed.
#' @param removeTPIsomers If \code{TRUE} then all TPs with equal formula as any sibling TPs (isomers) are removed.
#'   Unlike \code{removeDuplicates}, \emph{all} TP candidates are removed (including the first match). This filter
#'   automatically sets \code{removeDuplicates=TRUE} to avoid complete removal of TPs with equal structure.
#' @param minSimilarity Minimum structure similarity (\samp{0-1}) that a TP should have relative to its parent. For
#'   details on how these similarities are calculated, see the \code{\link{generateTPsBioTransformer}} function. May be
#'   useful under the assumption that parents and TPs who have a high structural similarity, also likely have a high
#'   MS/MS spectral similarity (which can be evaluated after componentization with \code{\link{generateComponentsTPs}}.
#' @param negate If \code{TRUE} then filters are performed in opposite manner.
#'
#' @return \code{filter} returns a filtered \code{transformationProductsBT} object.
#'
#' @export
setMethod("filter", "transformationProductsBT", function(obj, removeParentIsomers = FALSE, removeTPIsomers = FALSE,
                                                         removeDuplicates = FALSE, minSimilarity = NULL, negate = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertFlag, . ~ removeParentIsomers + removeTPIsomers + removeDuplicates + negate,
           fixed = list(add = ac))
    checkmate::assertNumber(minSimilarity, lower = 0, finite = TRUE, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    if (length(obj) == 0)
        return(obj)

    if (removeTPIsomers)
        removeDuplicates <- TRUE
    
    oldn <- length(obj)

    hash <- makeHash(obj, removeParentIsomers, removeTPIsomers, removeDuplicates, minSimilarity, negate)
    cache <- loadCacheData("filterTPs", hash)
    if (!is.null(cache))
        obj <- cache
    else
    {
        if (removeParentIsomers || removeTPIsomers || removeDuplicates)
        {
            # NOTE: obj@products should be first arg to Map to keep names...
            obj@products <- Map(obj@products, parents(obj)$formula, f = function(prod, pform)
            {
                if (removeDuplicates)
                    prod <- if (negate) prod[duplicated(SMILES)] else prod[!duplicated(SMILES)]
                if (removeParentIsomers)
                    prod <- if (negate) prod[formula == pform] else prod[formula != pform]
                if (removeTPIsomers)
                {
                    df <- getDuplicatedStrings(prod$formula)
                    prod <- if (negate) prod[formula %chin% df] else prod[!formula %chin% df]
                }
                return(prod)
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
