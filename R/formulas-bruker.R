#' @include main.R
#' @include utils-bruker.R
NULL

simplifyDAFormula <- function(formula)
{
    # Remove spaces
    formula <- gsub(" ", "", formula)

    # Removes any element counts of 1: match all "1"'s after a character but not before another number
    return(gsub("(?<=[[:alpha:]])1{1}(?![[:digit:]])", "", formula, perl = TRUE))
}

#' @details \code{generateFormulasDA} uses Bruker DataAnalysis to generate
#'   chemical formulae. This method supports scoring based on overlap between
#'   measured and theoretical isotopic patterns (both MS and MS/MS data) and the
#'   presence of 'fitting' MS/MS fragments. The method will iterate through all
#'   features (or "Compounds" in DataAnalysis terms) and call
#'   \command{SmartFormula} (and \command{SmartFormula3D} if MS/MS data is
#'   available) to generate all formulae. Parameters affecting formula
#'   calculation have to be set in advance within the DataAnalysis method for
#'   each analysis (\emph{e.g.} by \code{\link{setDAMethod}}). This method
#'   requires that features were obtained with
#'   \code{\link{findFeaturesBruker}}. Unlike other algorithms, there is no
#'   prior need to generate \code{\link[=MSPeakLists]{MS peak lists}}.
#'
#' @param precursorMzSearchWindow Search window for \emph{m/z} values (+/- the
#'   feature \emph{m/z}) used to find back feature data of precursor/parent ions
#'   from MS/MS spectra (this data is not readily available from
#'   \command{SmartFormula3D} results).
#'
#' @template dasaveclose-args
#'
#' @template DA-restart-note
#'
#' @rdname formula-generation
#' @export
setMethod("generateFormulasDA", "featureGroups", function(fGroups, precursorMzSearchWindow = 0.002,
                                                          MSMode = "both", adduct, featThreshold = 0,
                                                          featThresholdAnn = 0.75, save = TRUE, close = save)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    checkmate::assertNumber(precursorMzSearchWindow, lower = 0, finite = TRUE, add = ac)
    checkmate::assertChoice(MSMode, c("ms", "msms", "both"), add = ac)
    aapply(checkmate::assertNumber, . ~ featThreshold + featThresholdAnn, lower = 0, upper = 1, fixed = list(add = ac))
    assertDACloseSaveArgs(close, save, add = ac)
    checkmate::reportAssertions(ac)

    adduct <- checkAndToAdduct(adduct)

    DA <- getDAApplication()
    anaInfo <- analysisInfo(fGroups)
    ftind <- groupFeatIndex(fGroups)
    gTable <- groupTable(fGroups)
    gCount <- length(fGroups)
    gNames <- names(fGroups)

    hideDAInScope()

    fts <- featureTable(fGroups)
    fTable <- list()

    cacheDB <- openCacheDBScope() # open manually so caching code doesn't need to on each R/W access
    setHash <- makeHash(fGroups, precursorMzSearchWindow, MSMode)
    cachedSet <- loadCacheSet("formulasBruker", setHash, cacheDB)
    formHashes <- vector("character", nrow(anaInfo) * gCount)
    formHashCount <- 0

    for (anai in seq_len(nrow(anaInfo)))
    {
        ana <- anaInfo$analysis[anai]
        DAAnaInd <- getDAFileIndex(DA, ana, anaInfo$path[anai])
        checkDAFMFCompounds(DA, fts[[ana]], DAAnaInd, TRUE)

        printf("Loading all formulas from analysis '%s'...\n", ana)

        baseHash <- makeHash(fGroups, ana, precursorMzSearchWindow, MSMode)

        cmpds <- DA[["Analyses"]][[DAAnaInd]][["Compounds"]]

        prog <- openProgBar(0, gCount)

        for (grpi in seq_len(gCount))
        {
            grp <- gNames[grpi]
            fti <- ftind[[grp]][anai]

            if (fti == 0)
                next

            hash <- makeHash(baseHash, grp)
            formHashCount <- formHashCount + 1
            formHashes[formHashCount] <- hash

            flist <- NULL
            if (!is.null(cachedSet))
                flist <- cachedSet[[hash]]
            if (is.null(flist))
                flist <- loadCacheData("formulasBruker", hash, cacheDB)

            if (is.null(flist))
            {
                cmpd <- cmpds[[fts[[ana]]$ID[fti]]]

                ftable <- vector("list", gCount)
                ftableCount <- 0

                if (MSMode != "msms")
                {
                    cmpd$SmartFormula()

                    SMFResults <- cmpd[[1]][["SmartFormulaResults"]]
                    SMFResultsCount <- SMFResults[["Count"]]

                    if (SMFResultsCount > 0)
                    {
                        for (k in 1:SMFResultsCount)
                        {
                            SMFResult <- SMFResults[[k]]
                            SMFResultCount <- SMFResult[["Count"]]

                            if (SMFResultCount < 1 || all.equal(SMFResult[["m_over_z"]], fts[[ana]]$mz[fti]) != TRUE)
                                next

                            dt <- data.table(neutral_formula = character(SMFResultCount),
                                             formula = character(SMFResultCount),
                                             formula_mz = numeric(SMFResultCount),
                                             error = numeric(SMFResultCount), mSigma = numeric(SMFResultCount),
                                             score = numeric(SMFResultCount), byMSMS = FALSE)

                            for (l in 1:SMFResultCount)
                            {
                                SMFResultItem <- SMFResult[[l]]
                                dt[l, c("neutral_formula", "formula", "formula_mz", "error", "mSigma", "dbe", "score") :=
                                       .(simplifyDAFormula(SMFResultItem[["NeutralFormula"]]),
                                         simplifyDAFormula(SMFResultItem[["SumFormula"]]), SMFResultItem[["m_over_z"]],
                                         SMFResultItem[["Error"]], SMFResultItem[["Sigma"]] * 1000,
                                         SMFResultItem[["RingsAndDoubleBonds"]], SMFResultItem[["Score"]])]
                            }

                            ftableCount <- ftableCount + 1
                            ftable[[ftableCount]] <- dt
                            break # Shouldn't be any more results
                        }
                    }
                }

                ftable3D <- vector("list", gCount)
                ftable3DCount <- 0

                if (MSMode != "ms" && cmpd[["Count"]] > 1)
                {
                    cmpd$SmartFormula3D()

                    SMF3DResults <- cmpd[[2]][["SmartFormula3DResults"]]
                    SMF3DResultsCount <- SMF3DResults[["Count"]]

                    if (SMF3DResultsCount > 0)
                    {
                        ftable3DFrag <- vector("list", SMF3DResultsCount)
                        ftable3DFragCount <- 0

                        for (k in 1:SMF3DResultsCount)
                        {
                            SMF3DResult <- SMF3DResults[[k]]
                            SMF3DResultCount <- SMF3DResult[["Count"]]
                            parent <- SMF3DResult[["ParentResults"]]

                            if (SMF3DResultCount < 1 || abs(parent[["m_over_z"]] - fts[[ana]]$mz[fti]) > precursorMzSearchWindow)
                                next

                            form <- simplifyDAFormula(parent[["SumFormula"]])
                            nform <- calculateNeutralFormula(form, adduct)

                            dt <- data.table(neutral_formula = nform, formula = form, formula_mz = parent[["m_over_z"]],
                                             error = parent[["Error"]], mSigma = parent[["Sigma"]] * 1000,
                                             dbe = parent[["RingsAndDoubleBonds"]], score = parent[["Score"]],
                                             byMSMS = TRUE, frag_formula = character(SMF3DResultCount),
                                             frag_mz = numeric(SMF3DResultCount),
                                             frag_error = numeric(SMF3DResultCount), frag_mSigma = numeric(SMF3DResultCount),
                                             neutral_loss = character(SMF3DResultCount), frag_dbe = numeric(SMF3DResultCount),
                                             frag_score = numeric(SMF3DResultCount))

                            for (l in 1:SMF3DResultCount)
                            {
                                frag <- SMF3DResult[[l]]
                                fform <- simplifyDAFormula(frag[["SumFormula"]])
                                dt[l, c("frag_formula", "frag_mz", "frag_error", "frag_mSigma", "neutral_loss",
                                        "frag_dbe", "frag_score") :=
                                       .(fform, frag[["m_over_z"]], frag[["Error"]], frag[["Sigma"]] * 1000,
                                         subtractFormula(form, fform), frag[["RingsAndDoubleBonds"]], frag[["Score"]])]
                            }


                            ftable3DFragCount <- ftable3DFragCount + 1
                            ftable3DFrag[[ftable3DFragCount]] <- dt
                        }

                        if (ftable3DFragCount > 0)
                        {
                            ftable3DCount <- ftable3DCount + 1
                            ftable3D[[ftable3DCount]] <- rbindlist(ftable3DFrag[1:ftable3DFragCount])
                        }
                    }
                }

                flist <- list()

                if (ftableCount > 0)
                    flist <- rbindlist(ftable[1:ftableCount])

                if (ftable3DCount > 0)
                {
                    MSMSFlist <- rbindlist(ftable3D[1:ftable3DCount])
                    if (ftableCount > 0)
                        flist <- flist[!formula %in% MSMSFlist[["formula"]]]
                    flist <- rbind(flist, MSMSFlist, fill = TRUE)
                }

                saveCacheData("formulasBruker", flist, hash, cacheDB)
            }

            if (!is.null(flist) && length(flist) > 0)
            {
                if (is.null(fTable[[ana]]))
                    fTable[[ana]] <- list()

                fTable[[ana]][[grp]] <- flist
            }

            if ((grpi %% 5) == 0)
                setTxtProgressBar(prog, grpi)
        }

        setTxtProgressBar(prog, gCount)
        close(prog)

        closeSaveDAFile(DA, DAAnaInd, close, save)

        ngrp <- length(fTable[[ana]])
        printf("Loaded %d formulas for %d features (%.2f%%).\n", countUniqueFormulas(fTable[[ana]]),
               ngrp, if (gCount == 0) 0 else ngrp * 100 / gCount)
    }

    if (is.null(cachedSet))
        saveCacheSet("formulasBruker", formHashes[1:formHashCount], setHash, cacheDB)

    fTable <- pruneList(sapply(fTable, function(ft) ft[sapply(ft, nrow) > 0], simplify = FALSE), TRUE)

    if (length(fTable) > 0)
        groupFormulas <- generateGroupFormulasByConsensus(fTable, lapply(ftind, function(x) sum(x > 0)),
                                                          featThreshold, featThresholdAnn,
                                                          gNames, "analysis_from", "analyses", "featCoverage",
                                                          "featCoverageAnn")
    else
        groupFormulas <- list()

    return(formulas(formulas = groupFormulas, featureFormulas = fTable, algorithm = "bruker"))
})

setMethod("generateFormulasDA", "featureGroupsSet", function(fGroups, ..., setThreshold = 0, setThresholdAnn = 0.75)
{
    generateFormulasSet(fGroups, generateFormulasDA, ..., setThreshold = setThreshold,
                        setThresholdAnn = setThresholdAnn)
})
