#' @include main.R
#' @include compounds.R
#' @include utils-sirius.R
#' @include feature_groups-set.R
NULL

processSIRIUSCompounds <- function(msFName, outPath, MSMS, database, adduct, topMost)
{
    resultPath <- getSiriusResultPath(outPath, msFName)
    summary <- file.path(resultPath, "structure_candidates.tsv")
    results <- scRanges <- data.table()
    
    if (length(summary) != 0 && file.exists(summary))
    {
        results <- fread(summary)
        results <- unifySirNames(results)
        
        # NOTE: so far SIRIUS only has one score
        if (nrow(results) > 0)
            scRanges <- list(score = range(results$score))
        
        # manually sort results: SIRIUS seems to sometimes randomly order results with equal scores
        # NOTE: do this before topMost filter step below to ensure proper filtering
        setorderv(results, c("score", "identifier"), order = c(-1, 1))
        
        if (!is.null(topMost))
        {
            if (nrow(results) > topMost)
                results <- results[seq_len(topMost)] # results should already be sorted on score
        }
        
        # identifiers to characters, as they sometimes are and sometimes aren't depending on if multiple are present
        results[, identifier := as.character(identifier)]

        # NOTE: fragment info is based on SIRIUS results, ie from formula
        # prediction and not by compounds! Hence, results are the same for
        # candidates with the same formula.
        fragFiles <- getSiriusFragFiles(resultPath)
        for (ff in fragFiles)
        {
            fragInfo <- fread(ff)
            fragInfo[, c("rel.intensity", "exactmass", "intensity") := NULL]
            fragInfo[, ionization := gsub(" ", "", ionization)]
            fragInfo[, PLID := sapply(mz, function(omz) MSMS[which.min(abs(omz - mz))]$ID)]

            # each frag file always contains the precursor (even in input doesn't) --> use this to figure out which
            # candidate(s) it belongs to
            wh <- which(results$neutral_formula %in% fragInfo$formula) # UNDONE: check if it's really the precursor?
            
            if (length(wh) > 0)
            {
                # sirius neutralizes fragments, make them ion again
                if (!is.null(fragInfo[["implicitAdduct"]]))
                    ionImpAdducts <- ifelse(nzchar(fragInfo$implicitAdduct),
                                            mapply(paste0("+", fragInfo$implicitAdduct, "]"), fragInfo$ionization,
                                                   FUN = sub, MoreArgs = list(pattern = "\\]")),
                                            fragInfo$ionization)
                else
                    ionImpAdducts <- fragInfo$ionization
                setnames(fragInfo, "formula", "formula_SIR")
                fragInfo[, ion_formula := mapply(formula_SIR, ionImpAdducts, FUN = calculateIonFormula)]
                
                ionform <- calculateIonFormula(results$neutral_formula[wh[1]], adduct)
                fragInfo[, neutral_loss := sapply(ion_formula, subtractFormula, formula1 = ionform)]
                
                set(results, wh, "fragInfo", list(list(fragInfo)))
                set(results, wh, "explainedPeaks", nrow(fragInfo))
            }
        }
        
        if (is.null(results[["fragInfo"]]))
        {
            # warning(sprintf("no fragment info for %s", cmd$gName))
            results[, fragInfo := list(rep(list(data.table(mz = numeric(0), ion_formula = character(0),
                                                           neutral_loss = character(0),  score = numeric(0),
                                                           PLID = numeric(0))),
                    nrow(results)))]
            results[, explainedPeaks := 0]
        }
        
        results[, database := database][]
    }
    
    ret <- list(comptab = results, scRanges = scRanges)
    return(ret)
}

#' Compound annotation with SIRIUS
#'
#' Uses \href{https://bio.informatik.uni-jena.de/software/sirius/}{SIRIUS} in combination with
#' \href{https://www.csi-fingerid.uni-jena.de/}{CSI:FingerID} for compound annotation.
#'
#' @templateVar algo SIRIUS
#' @templateVar do generate compound candidates
#' @templateVar generic generateCompounds
#' @templateVar algoParam sirius
#' @template algo_generator
#'
#' @details Similar to \code{\link{generateFormulasSIRIUS}}, candidate formulae are generated with SIRIUS. These results
#'   are then feed to CSI:FingerID to acquire candidate structures. This method requires the availability of MS/MS data,
#'   and feature groups without it will be ignored.
#'
#' @param fingerIDDatabase Database specifically used for \command{CSI:FingerID}. If \code{NULL}, the value of the
#'   \code{formulaDatabase} parameter will be used or \code{"pubchem"} when that is also \code{NULL}. Sets the
#'   \option{--fingerid-db} option.
#' @param topMostFormulas Do not return more than this number of candidate formulae. Note that only compounds for these
#'   formulae will be searched. Sets the \option{--candidates} commandline option.
#' @param verbose If \code{TRUE} then more output is shown in the terminal.
#'
#' @templateVar ident TRUE
#' @template sirius-args
#' @template sirius_form-args
#' @template adduct-arg
#' @template comp_algo-args
#' 
#' @inheritParams generateCompounds
#'
#' @inherit generateCompounds return
#' 
#' @templateVar what \code{generateCompoundsSIRIUS}
#' @template uses-multiProc
#'
#' @templateVar what generateCompoundsSIRIUS
#' @template main-rd-method
#' @export
setMethod("generateCompoundsSIRIUS", "featureGroups", function(fGroups, MSPeakLists, relMzDev = 5, adduct = NULL,
                                                               elements = "CHNOP",
                                                               profile = "qtof", formulaDatabase = NULL,
                                                               fingerIDDatabase = "pubchem", noise = NULL, errorRetries = 2,
                                                               cores = NULL, topMost = 100, topMostFormulas = 5,
                                                               extraOptsGeneral = NULL, extraOptsFormula = NULL, verbose = TRUE,
                                                               splitBatches = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertNumber(relMzDev, lower = 0, finite = TRUE, add = ac)
    aapply(checkmate::assertString, . ~ elements + profile + fingerIDDatabase, fixed = list(add = ac))
    aapply(checkmate::assertString, . ~ formulaDatabase + fingerIDDatabase, null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertNumber(noise, lower = 0, finite = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertCount(errorRetries, add = ac)
    checkmate::assertCount(cores, positive = TRUE, null.ok = TRUE, add = ac)
    aapply(checkmate::assertCount, . ~ topMost + topMostFormulas, positive = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertCharacter, . ~ extraOptsGeneral + extraOptsFormula, null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertFlag(verbose, add = ac)
    checkmate::assertFlag(splitBatches, add = ac)
    checkmate::reportAssertions(ac)

    if (length(fGroups) == 0)
        return(compounds(algorithm = "sirius"))
    
    adduct <- checkAndToAdduct(adduct, fGroups)
    if (is.null(fingerIDDatabase))
        fingerIDDatabase <- if (!is.null(formulaDatabase)) formulaDatabase else "pubchem"

    gCount <- length(fGroups)
    printf("Processing %d feature groups with SIRIUS-CSI:FingerID...\n", gCount)
    
    results <- doSIRIUS(fGroups, MSPeakLists, FALSE, profile, adduct, relMzDev, elements,
                        formulaDatabase, noise, cores, TRUE, fingerIDDatabase, topMost, extraOptsGeneral, extraOptsFormula,
                        verbose, "compoundsSIRIUS", patRoon:::processSIRIUSCompounds,
                        list(database = fingerIDDatabase, topMost = topMost),
                        splitBatches)
    
    # prune empty/NULL results
    if (length(results) > 0)
        results <- results[sapply(results, function(r) !is.null(r$comptab) && nrow(r$comptab) > 0, USE.NAMES = FALSE)]
    
    if (verbose)
    {
        printf("-------\n")
        ngrp <- length(results)
        printf("Assigned %d compounds to %d feature groups (%.2f%%).\n", sum(unlist(lapply(results, function(r) nrow(r$comptab)))),
               ngrp, if (gCount == 0) 0 else ngrp * 100 / gCount)
    }

    return(compounds(groupAnnotations = lapply(results, "[[", "comptab"), scoreTypes = "score",
                     scoreRanges = lapply(results, "[[", "scRanges"),
                     algorithm = "sirius"))
})

#' @template featAnnSets-gen_args
#' @rdname generateCompoundsSIRIUS
#' @export
setMethod("generateCompoundsSIRIUS", "featureGroupsSet", function(fGroups, MSPeakLists, relMzDev = 5, adduct = NULL,
                                                                  ..., setThreshold = 0, setThresholdAnn = 0)
{
    generateCompoundsSet(fGroups, MSPeakLists, adduct, generateCompoundsSIRIUS, relMzDev = relMzDev, ...,
                         setThreshold = setThreshold, setThresholdAnn = setThresholdAnn)
})
