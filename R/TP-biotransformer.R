#' @include main.R
#' @include TP-structure.R
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
transformationProductsBT <- setClass("transformationProductsBT", contains = "transformationProductsStructure")

setMethod("initialize", "transformationProductsBT",
          function(.Object, ...) callNextMethod(.Object, algorithm = "biotransformer", ...))


getBaseBTCmd <- function(parent, SMILES, type, steps, extraOpts, baseHash)
{
    mainArgs <- c("-b", type,
                  "-k", "pred",
                  "-ismi", SMILES,
                  "-s", as.character(steps),
                  extraOpts)
    
    return(list(command = "java", args = mainArgs, logFile = paste0("biotr-", parent, ".txt"), parent = parent,
                SMILES = SMILES, hash = makeHash(parent, SMILES, baseHash)))
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
    # Simplify/harmonize columns a bit
    setnames(ret,
             c("Molecular formula", "Major Isotope Mass", "Reaction", "Reaction ID", "Metabolite ID", "Precursor ID",
               "Precursor ALogP", "Enzyme(s)", "Biosystem"),
             c("formula", "neutralMass", "transformation", "transformation_ID", "ID", "parent_ID", "parent_ALogP",
               "enzyme", "biosystem"))
    ret[!nzchar(parent_ID), parent_ID := NA]
    for (col in c("ID", "parent_ID"))
        set(ret, i = NULL, j = col, value = as.integer(sub("^BTM", "", ret[[col]])))
    
    # No need for these...
    # NOTE: cdk:Title seems the same as "Metabolite ID" column(?)
    ret[, c("Synonyms", "PUBCHEM_CID", "cdk:Title") := NULL]
    ret[, (grep("^Precursor ", names(ret), value = TRUE)) := NULL]
    
    # BUG: BT sometimes doesn't fill in the formula. Calculate them manually
    ret[!nzchar(formula), formula := {
        mols <- patRoon:::getMoleculesFromSMILES(SMILES)
        return(sapply(mols, function(m) rcdk::get.mol2formula(m)@string))
    }]
    
    # Assign some unique identifier
    ret[, name := paste0(cmd$parent, "-TP", ID)]

    # NOTE: take the _original_ parent ALogP as reference
    parALogP <- ret[is.na(parent_ID)]$parent_ALogP[1]
    ret[, retDir := fifelse(ALogP < parALogP, -1, 1)]

    setcolorder(ret, c("name", "ID", "parent_ID", "SMILES", "InChI", "InChIKey", "formula", "neutralMass"))
    
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
                                      skipInvalid = TRUE, calcSims = FALSE, fpType = "extended",
                                      fpSimMethod = "tanimoto", MP = FALSE)
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
    aapply(checkmate::assertFlag, . ~ skipInvalid + calcSims + MP, fixed = list(add = ac))
    aapply(checkmate::assertString, . ~ fpType + fpSimMethod, min.chars = 1, fixed = list(add = ac))
    checkmate::reportAssertions(ac)

    parents <- getTPParents(parents, skipInvalid)

    baseHash <- makeHash(type, steps, extraOpts, skipInvalid, fpType, fpSimMethod)
    setHash <- makeHash(parents, baseHash)
    
    cmdQueue <- Map(parents$name, parents$SMILES, f = getBaseBTCmd,
                    MoreArgs = list(type = type, steps = steps, extraOpts = extraOpts, baseHash = baseHash))

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

    return(transformationProductsBT(calcSims = calcSims, fpType = fpType, fpSimMethod = fpSimMethod, parents = parents,
                                    products = results))
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

    # NOTE: need to override method since BT results need to be collapsed
    cat("Collapsing results... ")
    prodAll <- if (length(TPs) > 0) collapseBTResults(TPs@products) else data.table()
    cat("Done!\n")

    doConvertToMFDB(prodAll, parents(TPs), out, includeParents)
})
