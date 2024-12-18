# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
#' @include TP-structure.R
NULL

#' @rdname transformationProductsStructure-class
transformationProductsBT <- setClass("transformationProductsBT", contains = "transformationProductsStructure")

setMethod("initialize", "transformationProductsBT",
          function(.Object, ...) callNextMethod(.Object, algorithm = "biotransformer", ...))


getBaseBTCmd <- function(parent, SMILES, type, generations, extraOpts, baseHash, maxExpGenerations, neutralizeTPs)
{
    mainArgs <- c("-b", type,
                  "-k", "pred",
                  "-ismi", SMILES,
                  "-s", as.character(generations),
                  extraOpts)
    
    return(list(command = "java", args = mainArgs, logFile = paste0("biotr-", parent, ".txt"), parent = parent,
                SMILES = SMILES, maxExpGenerations = maxExpGenerations, neutralizeTPs = neutralizeTPs,
                hash = makeHash(parent, SMILES, baseHash)))
}

BTMPFinishHandler <- function(cmd)
{
    if (!file.exists(cmd$outFile))
        return(data.table()) # no results
    
    results <- fread(cmd$outFile, colClasses = c("Precursor ID" = "character", Synonyms = "character",
                                                 "Molecular formula" = "character"))
    
    # Simplify/harmonize columns a bit
    setnames(results,
             c("Molecular formula", "Major Isotope Mass", "Reaction", "Reaction ID", "Metabolite ID", "Precursor ID",
               "Precursor ALogP", "Enzyme(s)", "Biosystem"),
             c("formula", "neutralMass", "transformation", "transformation_ID", "chem_ID", "parent_chem_ID", "parent_ALogP",
               "enzyme", "biosystem"))
    results[!nzchar(parent_chem_ID), parent_chem_ID := 0L]
    for (col in c("chem_ID", "parent_chem_ID"))
        set(results, i = NULL, j = col, value = as.integer(sub("^BTM", "", results[[col]])))
    # results[is.na(parent_chem_ID), parent_chem_ID := 0L]
    
    # BT only specifies to which chemical structure a TP is, not if it the parent came from a specific route. For now
    # assume the TP could be from any of the parents.
    curTPID <- 0
    processChilds <- function(parID, parChemID, generation)
    {
        if (generation > cmd$maxExpGenerations)
            return(data.table())
        resSub <- copy(results[parent_chem_ID == parChemID])
        resSub[, c("ID", "parent_ID", "generation") := .(curTPID + seq_len(nrow(resSub)), parID, generation)]
        curTPID <<- curTPID + nrow(resSub)
        return(rbindlist(c(list(resSub), Map(resSub$ID, resSub$chem_ID, f = processChilds,
                                             MoreArgs = list(generation = generation + 1)))))
    }
    
    ret <- processChilds(0, 0, 1)
    ret[parent_ID == 0, parent_ID := NA_integer_]

    # No need for these...
    # NOTE: cdk:Title seems the same as "Metabolite ID" column(?)
    ret[, c("Synonyms", "PUBCHEM_CID", "cdk:Title", "parent_chem_ID") := NULL]
    ret[, (grep("^Precursor ", names(ret), value = TRUE)) := NULL]
    
    # BUG: BT sometimes doesn't fill in the formula. Calculate them manually
    ret[!nzchar(formula), formula := {
        mols <- patRoon:::getMoleculesFromSMILES(SMILES)
        return(sapply(mols, function(m) rcdk::get.mol2formula(m)@string))
    }]
    
    # Assign some unique identifier
    ret[, name := paste0(cmd$parent, "-TP", chem_ID)]

    # NOTE: take the _original_ parent ALogP as reference
    parALogP <- ret[is.na(parent_ID)]$parent_ALogP[1]
    ret[, retDir := fifelse(ALogP < parALogP, -1, 1)]

    setcolorder(ret, c("name", "ID", "parent_ID", "chem_ID", "SMILES", "InChI", "InChIKey", "formula", "neutralMass"))
    
    ret <- prepareChemTable(ret, prefCalcChemProps = FALSE, neutralChemProps = cmd$neutralizeTPs, verbose = FALSE)
    
    return(ret)
}

BTMPPrepareHandler <- function(cmd)
{
    btBin <- getExtDepPath("biotransformer")
    
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
#' @details In order to use this function the \file{.jar} command line utility should be installed and specified in the
#'   \code{\link[=patRoon-package]{patRoon.path.BioTransformer}} option. The \file{.jar} file can be obtained via
#'   \url{https://bitbucket.org/djoumbou/biotransformer/src/master}. Alternatively, the \pkg{patRoonExt} package can be
#'   installed to automatically install/configure the necessary files.
#'
#' @param type The type of prediction. Valid values are: \code{"env"}, \code{"ecbased"}, \code{"cyp450"},
#'   \code{"phaseII"}, \code{"hgut"}, \code{"superbio"}, \code{"allHuman"}. Sets the \command{-b} command line option.
#' @param generations The number of generations (steps) for the predictions. Sets the \command{-s} command line option.
#'   More generations may be reported, see the \verb{Hierarchy expansion} section below.
#' @param maxExpGenerations The maximum number of generations during hierarchy expansion, see below.
#' @param extraOpts A \code{character} with extra command line options passed to the \command{biotransformer.jar} tool.
#' @param MP If \code{TRUE} then multiprocessing is enabled. Since \command{BioTransformer} supports native
#'   parallelization, additional multiprocessing generally doesn't lead to significant reduction in computational times.
#'   Furthermore, enabling multiprocessing can lead to very high CPU/RAM usage.
#'
#' @return The TPs are stored in an object derived from the \code{\link{transformationProductsStructure}} class.
#'
#' @section Hierarchy expansion: \command{BioTransformer} only reports the direct parent for a TP, not
#'   the complete pathway. For instance, consider the following results: \itemize{
#'
#'   \item parent --> TP1
#'
#'   \item parent --> TP2
#'
#'   \item TP1 --> TP2
#'
#'   \item TP2 --> TP3
#'
#'   }
#'   
#'   In this case, TP3 may be formed either as: \itemize{
#'   
#'   \item parent --> TP1 --> TP2 --> TP3
#'   
#'   \item parent --> TP2 --> TP3
#'   
#'   }
#'   
#'   For this reason, \pkg{patRoon} simply expands the hierarchy and assumes that all routes are possible. For instance,
#'   \verb{
#'      Parent    
#'      /-  -\    
#'    /-      -\  
#'   -          - 
#'  TP1        TP2
#'   |          | 
#'   |          | 
#'  TP2        TP3
#'   |            
#'   |            
#'  TP3
#'   }
#'   
#'   Note that this may result in pathways with more generations than defined by the \code{generations} argument. Thus,
#'   the \code{maxExpGenerations} argument is used to avoid excessive expansions.
#'
#' @template tp_gen-scr
#' @template tp_gen-sim
#' @template fp-args
#'
#' @templateVar whatCP parent suspect list
#' @template chemPropCalc
#'
#' @templateVar what \code{generateTPsBioTransformer}
#' @template uses-multiProc
#'
#' @references \insertRef{DjoumbouFeunang2019}{patRoon} \cr\cr \insertRef{Wicker2015}{patRoon}
#'
#' @export
generateTPsBioTransformer <- function(parents, type = "env", generations = 2, maxExpGenerations = generations + 2,
                                      extraOpts = NULL, skipInvalid = TRUE, prefCalcChemProps = TRUE,
                                      neutralChemProps = FALSE, neutralizeTPs = TRUE, calcSims = FALSE,
                                      fpType = "extended", fpSimMethod = "tanimoto", MP = FALSE)
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
    checkmate::assertCount(generations, positive = TRUE, add = ac)
    checkmate::assertCount(maxExpGenerations, positive = TRUE, add = ac)
    checkmate::assertCharacter(extraOpts, null.ok = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ skipInvalid + prefCalcChemProps + neutralChemProps + neutralizeTPs + calcSims + MP,
           fixed = list(add = ac))
    aapply(checkmate::assertString, . ~ fpType + fpSimMethod, min.chars = 1, fixed = list(add = ac))
    checkmate::reportAssertions(ac)

    parents <- getTPParents(parents, skipInvalid, prefCalcChemProps, neutralChemProps)

    baseHash <- makeHash(type, generations, maxExpGenerations, extraOpts, prefCalcChemProps, neutralChemProps,
                         neutralizeTPs, skipInvalid, fpType, fpSimMethod)
    setHash <- makeHash(parents, baseHash)
    
    cmdQueue <- Map(parents$name, parents$SMILES, f = getBaseBTCmd,
                    MoreArgs = list(type = type, generations = generations, maxExpGenerations = maxExpGenerations,
                                    neutralizeTPs = neutralizeTPs, extraOpts = extraOpts, baseHash = baseHash))

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
