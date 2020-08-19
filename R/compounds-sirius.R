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
    
    # NOTE: SIRIUS frag results are identified by 'neutral adducts', which is
    # the (de)protonated form of a formula with adduct, eg [M+NH4]+ yields M+NH3
    addChr <- as.character(adduct)
    neutralFormIsAdductForm <- addChr == "[M+H]+" || addChr == "[M-H]-" || addChr == "[M]"
    fragAdduct <- as.adduct(if (adduct@charge > 0) "[M+H]+" else "[M-H]-")
    
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
        
        if (neutralFormIsAdductForm)
            neutralAdductForms <- results$formula
        else
        {
            # get neutral formula and (de)protonate to get neutral adduct formula
            neutralAdductForms <- calculateNeutralFormula(calculateIonFormula(results$formula, adduct),
                                                          fragAdduct)
        }
        
        # NOTE: fragment info is based on SIRIUS results, ie from formula
        # prediction and not by compounds! Hence, results are the same for
        # candidates with the same formula.
        fragFiles <- getSiriusFragFiles(resultPath)
        for (ff in fragFiles)
        {
            neutralAdductFormFrag <- getFormulaFromSiriusFragFile(ff)
            
            if (neutralAdductFormFrag %in% neutralAdductForms) # may not be there if filtered out or no compound was found
            {
                fragInfo <- fread(ff)
                fragInfo[, c("rel.intensity", "exactmass", "ionization") := NULL]
                fragInfo[, PLIndex := sapply(mz, function(omz) which.min(abs(omz - MSMS$mz)))]
                
                wh <- which(neutralAdductForms == neutralAdductFormFrag)
                if (length(wh) > 0)
                {
                    fragInfo[, neutral_loss := sapply(formula, subtractFormula,
                                                      formula1 = results$formula[wh[1]])]
                    
                    # sirius neutralizes fragments, make them ion again
                    fragInfo[, formula := calculateIonFormula(formula, ..adduct)]
                    
                    set(results, wh, "fragInfo", list(list(fragInfo)))
                    set(results, wh, "explainedPeaks", nrow(fragInfo))
                }
            }
        }
        
        if (is.null(results[["fragInfo"]]))
        {
            # warning(sprintf("no fragment info for %s", cmd$gName))
            results[, fragInfo := list(rep(list(data.table()), nrow(results)))]
            results[, explainedPeaks := 0]
        }
        
        results[, database := database][]
    }
    
    ret <- list(comptab = results, scRanges = scRanges)
    return(ret)
}

#' @details \code{generateCompoundsSIRIUS} uses
#'   \href{https://bio.informatik.uni-jena.de/software/sirius/}{SIRIUS} in
#'   combination with \href{https://www.csi-fingerid.uni-jena.de/}{CSI:FingerID}
#'   for compound identification. Similar to
#'   \code{\link{generateFormulasSIRIUS}}, candidate formulae are generated with
#'   SIRIUS. These results are then feed to CSI:FingerID to acquire candidate
#'   structures. This method requires the availability
#'   of MS/MS data, and feature groups without it will be ignored.
#'
#' @templateVar ident TRUE
#' @template sirius-args
#'
#' @templateVar genForm FALSE
#' @template form-args
#'
#' @return \code{generateCompoundsSIRIUS} returns a \code{\link{compounds}}
#'   object.
#' @param fingerIDDatabase Database specifically used for
#'   \command{CSI:FingerID}. If \code{NULL}, the value of the
#'   \code{formulaDatabase} parameter will be used or \code{"pubchem"} when that
#'   is also \code{NULL}. Sets the \option{--fingerid-db} option.
#' @param topMostFormulas Do not return more than this number of candidate
#'   formulae. Note that only compounds for these formulae will be searched.
#'   Sets the \option{--candidates} commandline option.
#' @param verbose If \code{TRUE} then more output is shown in the terminal.
#'
#' @rdname compound-generation
#' @export
setMethod("generateCompoundsSIRIUS", "featureGroups", function(fGroups, MSPeakLists, relMzDev = 5, adduct = "[M+H]+",
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

    adduct <- checkAndToAdduct(adduct)
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

    return(compounds(compounds = lapply(results, "[[", "comptab"), scoreTypes = "score",
                     scoreRanges = lapply(results, "[[", "scRanges"),
                     algorithm = "sirius"))
})

setMethod("generateCompoundsSIRIUS", "featureGroupsSet", function(fGroups, MSPeakLists, ..., setThreshold = 0.75)
{
    generateCompoundsSet(fGroups, MSPeakLists, generateCompoundsSIRIUS, ..., setThreshold = setThreshold)
})
