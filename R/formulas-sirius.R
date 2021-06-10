#' @include main.R
#' @include mspeaklists.R
#' @include formulas.R
#' @include utils-sirius.R
NULL

# callback for executeMultiProcess()
processSIRIUSFormulas <- function(msFName, outPath, adduct, ...)
{
    noResult <- forms <- data.table::data.table(neutral_formula = character(), ion_formula = character(),
                                                neutralMass = numeric(), ion_formula_mz = numeric(),
                                                adduct = character(), score = numeric(), MSMSScore = numeric(),
                                                isoScore = numeric(), explainedPeaks = integer(),
                                                explainedIntensity = numeric(), fragInfo = list())
    
    
    resultPath <- patRoon:::getSiriusResultPath(outPath, msFName)
    summary <- file.path(resultPath, "formula_candidates.tsv")
    if (length(summary) == 0 || length(summary) == 0 || !file.exists(summary))
        forms <- noResult
    else
    {
        forms <- fread(summary)
        fragFiles <- patRoon:::getSiriusFragFiles(resultPath)
        
        if (nrow(forms) == 0 || length(fragFiles) == 0)
            forms <- noResult
        else
        {
            data.table::setnames(forms,
                     c("molecularFormula", "precursorFormula", "SiriusScore",
                       "TreeScore", "IsotopeScore", "numExplainedPeaks", "medianMassErrorFragmentPeaks(ppm)",
                       "medianAbsoluteMassErrorFragmentPeaks(ppm)", "massErrorPrecursor(ppm)"),
                     c("neutral_formula", "neutral_adduct_formula", "score", "MSMSScore", "isoScore",
                       "explainedPeaks", "error_frag_median", "error_frag_median_abs", "error"))
            
            ionImpAdductsCached <- patRoon:::makeEmptyListNamed(list())
            fragInfoList <- lapply(fragFiles, function(ff)
            {
                fragInfo <- data.table::fread(ff)
                data.table::setnames(fragInfo,
                         c("exactmass", "formula"),
                         c("ion_formula_mz", "ion_formula_SIR"))
                fragInfo[, rel.intensity := NULL]
                fragInfo[, ionization := gsub(" ", "", ionization, fixed = TRUE)]

                # sirius neutralizes fragments, make them ion again
                if (!is.null(fragInfo[["implicitAdduct"]]))
                    ionImpAdducts <- ifelse(nzchar(fragInfo$implicitAdduct),
                                            mapply(paste0("+", fragInfo$implicitAdduct, "]"), fragInfo$ionization,
                                                   FUN = sub, MoreArgs = list(pattern = "\\]")),
                                            fragInfo$ionization)
                else
                    ionImpAdducts <- fragInfo$ionization
                
                notCached <- setdiff(ionImpAdducts, names(ionImpAdductsCached))
                if (length(notCached) > 0)
                    ionImpAdductsCached <<- c(ionImpAdductsCached, sapply(notCached, patRoon:::as.adduct,
                                                                          format = "sirius", simplify = FALSE))
                ionImpAdducts <- ionImpAdductsCached[ionImpAdducts]

                fragInfo[, ion_formula := mapply(ion_formula_SIR, ionImpAdducts, FUN = patRoon:::calculateIonFormula)]
                if (!is.null(fragInfo[["implicitAdduct"]]))
                {
                    # 'correct' formula masses: SIRIUS subtract implicit adduct from it
                    fragInfo[nzchar(implicitAdduct), ion_formula_mz := ion_formula_mz +
                                 sapply(implicitAdduct, patRoon:::getFormulaMass)]
                }
                return(fragInfo)
            })
            names(fragInfoList) <- sapply(fragFiles, patRoon:::getFormulaFromSiriusFragFile)
            
            # initialize all with empty fragInfos
            if (length(fragInfoList) > 0)
                forms[match(names(fragInfoList), neutral_adduct_formula), fragInfo := fragInfoList]
            else
                forms[, fragInfo := list()]

            forms <- patRoon:::addMiscFormulaInfo(forms, adduct)
            
            forms[, rank := NULL]
            
            # Precursor is always present in MS/MS spectrum: it's added by SIRIUS if necessarily (with zero intensity).
            # Remove it and use its mz to get ion_formula_mz
            forms[, ion_formula_mz := mapply(ion_formula, fragInfo,
                                             FUN = function(form, fi) fi$ion_formula_mz[form == fi$ion_formula])]
            forms[, fragInfo := Map(ion_formula, fragInfo, f = function(form, fi)
            {
                fi <- fi[intensity != 0 | ion_formula != form]
                fi[, intensity := NULL]
                return(fi)
            })]
            forms[, explainedPeaks := sapply(fragInfo, nrow)] # update
            
            # set nice column order
            data.table::setcolorder(forms, c("neutral_formula", "ion_formula", "neutral_adduct_formula", "neutralMass",
                                             "ion_formula_mz", "error", "error_frag_median", "error_frag_median_abs",
                                             "adduct", "score", "MSMSScore", "isoScore", "explainedPeaks",
                                             "explainedIntensity"))
            
            forms <- patRoon:::rankFormulaTable(forms)
        }
    }
    return(forms)
}

#' @details \code{generateFormulasSIRIUS} uses
#'   \href{https://bio.informatik.uni-jena.de/software/sirius/}{SIRIUS} to
#'   generate chemical formulae. Similarity of measured and theoretical isotopic
#'   patterns will be used for scoring candidates. Note that \command{SIRIUS}
#'   requires availability of MS/MS data.
#'
#' @templateVar ident FALSE
#' @template sirius-args
#' 
#' @param verbose If \code{TRUE} then more output is shown in the terminal.
#'
#' @rdname formula-generation
#' @export
setMethod("generateFormulasSIRIUS", "featureGroups", function(fGroups, MSPeakLists, relMzDev = 5,
                                                              adduct = NULL, elements = "CHNOP",
                                                              profile = "qtof", database = NULL, noise = NULL,
                                                              cores = NULL, topMost = 100, extraOptsGeneral = NULL,
                                                              extraOptsFormula = NULL, calculateFeatures = TRUE,
                                                              featThreshold = 0, featThresholdAnn = 0.75,
                                                              absAlignMzDev = 0.002, verbose = TRUE,
                                                              splitBatches = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertNumber(relMzDev, lower = 0, finite = TRUE, add = ac)
    aapply(checkmate::assertString, . ~ elements + profile, fixed = list(add = ac))
    checkmate::assertString(database, null.ok = TRUE, add = ac)
    checkmate::assertNumber(noise, lower = 0, finite = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertCount(cores, positive = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertCount(topMost, positive = TRUE, add = ac)
    aapply(checkmate::assertCharacter, . ~ extraOptsGeneral + extraOptsFormula, null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertFlag(calculateFeatures, add = ac)
    aapply(checkmate::assertNumber, . ~ featThreshold + featThresholdAnn + absAlignMzDev, lower = 0, upper = 1,
           fixed = list(add = ac))
    checkmate::assertFlag(verbose, add = ac)
    checkmate::assertFlag(splitBatches, add = ac)
    checkmate::reportAssertions(ac)
    
    if (length(fGroups) == 0)
        return(formulas(algorithm = "sirius"))
    
    adduct <- checkAndToAdduct(adduct, fGroups)
    gNames <- names(fGroups)
    gCount <- length(fGroups)
    
    printf("Processing %d feature groups with SIRIUS...\n---\n", gCount)
    formTable <- doSIRIUS(fGroups, MSPeakLists, calculateFeatures, profile, adduct, relMzDev, elements,
                          database, noise, cores, FALSE, NULL, topMost, extraOptsGeneral, extraOptsFormula,
                          verbose, "formulasSIRIUS", patRoon:::processSIRIUSFormulas, NULL,
                          splitBatches)
        
    if (calculateFeatures)
    {
        if (length(formTable) > 0)
        {
            formTable <- lapply(formTable, pruneList, checkZeroRows = TRUE)
            groupFormulas <- generateGroupFormulasByConsensus(formTable,
                                                              lapply(groupFeatIndex(fGroups), function(x) sum(x > 0)),
                                                              featThreshold, featThresholdAnn, gNames)
        }
        else
            groupFormulas <- list()
    }
    else
    {
        groupFormulas <- pruneList(formTable, checkZeroRows = TRUE)
        formTable <- list()
    }
    
    groupFormulas <- setFormulaPLID(groupFormulas, MSPeakLists, absAlignMzDev)
    
    if (verbose)
    {
        printf("-------\n")
        if (calculateFeatures && length(formTable) > 0)
        {
            fCounts <- sapply(formTable, countUniqueFormulas)
            fTotCount <- sum(fCounts)
            printf("Formula statistics:\n")
            printf("%s: %d (%.1f%%)\n", names(formTable), fCounts, if (fTotCount == 0) 0 else fCounts * 100 / fTotCount)
            printf("Total: %d\n", fTotCount)
        }
        ngrp <- length(groupFormulas)
        printf("Assigned %d unique formulas to %d feature groups (%.2f%% coverage).\n", countUniqueFormulas(groupFormulas),
               ngrp, if (gCount == 0) 0 else ngrp * 100 / gCount)
    }
    
    return(formulas(groupAnnotations = groupFormulas, featureFormulas = formTable, algorithm = "sirius"))
})

#' @rdname formula-generation
#' @export
setMethod("generateFormulasSIRIUS", "featureGroupsSet", function(fGroups, MSPeakLists, ..., setThreshold = 0,
                                                                 setThresholdAnn = 0.75)
{
    generateFormulasSet(fGroups, MSPeakLists, generateFormulasSIRIUS, ..., setThreshold = setThreshold,
                        setThresholdAnn = setThresholdAnn)
})
