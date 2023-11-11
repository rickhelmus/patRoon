#' @include main.R
#' @include mspeaklists.R
#' @include formulas.R
NULL

getGenFormBin <- function()
{
    gbin <- if (Sys.info()[["sysname"]] == "Windows") "GenForm.exe" else "GenForm"
    pOpt <- getOption("patRoon.path.GenForm")
    if (!is.null(pOpt) && nzchar(pOpt))
        ret <- file.path(pOpt, gbin)
    else
    {
        ret <- system.file("bin", Sys.getenv("R_ARCH"), gbin, package = "patRoon")
        
        # HACK: try GenForm binary from local src folder if devtools::load_all() was used
        if (!nzchar(ret) && requireNamespace("pkgload", quietly = TRUE) && !is.null(pkgload::dev_meta("patRoon")))
            ret <- normalizePath(file.path(system.file(".", package = "patRoon"), "..", "src", gbin))    }

    if (!file.exists(ret))
        stop(sprintf("GenForm binary does not exist: %s", ret))

    return(ret)
}

makeGenFormFilesNames <- function()
{
    MSFile <- tempfile("MSPList", fileext = ".txt")
    MSMSFile <- tempfile("MSMSPList", fileext = ".txt")
    outFile <- tempfile("formulas", fileext = ".txt")
    return(list(MSFile = MSFile, MSMSFile = MSMSFile, outFile = outFile))
}

writeGenFormFiles <- function(MSPList, MSMSPList, MSFile, MSMSFile)
{
    writePList <- function(f, pl) fwrite(pl[, c("mz", "intensity")], f, quote = FALSE,
                                         sep = "\t", row.names = FALSE, col.names = FALSE)

    writePList(MSFile, MSPList)

    if (!is.null(MSMSPList))
        writePList(MSMSFile, MSMSPList)
}

makeGenFormCmdQueue <- function(groupPeakLists, annTable, baseHash, MSMode, isolatePrec, adduct, topMost,
                                mainArgs, ana)
{
    gNames <- names(groupPeakLists)
    
    fgAdd <- getFGroupAdducts(gNames, annTable, adduct, "genform")
    
    pruneList(mapply(gNames, fgAdd$grpAdducts, fgAdd$grpAdductsChr, FUN = function(grp, add, addChr)
    {
        if (is.null(groupPeakLists[[grp]]) || is.null(groupPeakLists[[grp]][["MS"]]))
            return(NULL) # MS or MSMS spectrum probably filtered away
        
        hasMSMS <- !is.null(groupPeakLists[[grp]][["MSMS"]]) && MSMode != "ms"
        if (MSMode == "msms" && !hasMSMS)
            return(NULL)
        
        plms <- groupPeakLists[[grp]][["MS"]]
        plmz <- plms[precursor == TRUE, mz]
        if (length(plmz) == 0)
            return(NULL) # no precursor
        
        if (!isFALSE(isolatePrec))
        {
            if (is.logical(isolatePrec) && isolatePrec)
                isolatePrec <- getDefIsolatePrecParams(z = abs(add@charge))
            else
                isolatePrec$z = abs(add@charge)
            
            plms <- isolatePrecInMSPeakList(plms, isolatePrec, negate = FALSE)
        }

        args <- c(mainArgs, paste0("ion=", addChr))
        hash <- makeHash(groupPeakLists[[grp]], baseHash, isolatePrec, addChr)
               
        return(list(args = args, PLMZ = plmz, MSPL = plms, MSMSPL = groupPeakLists[[grp]][["MSMS"]],
                    hash = hash, group = grp, isMSMS = hasMSMS, adduct = add,
                    topMost = topMost, MSMode = MSMode, ana = ana))
    }, SIMPLIFY = FALSE))
}

GenFormMPFinishHandler <- function(cmd)
{
    f <- patRoon:::processGenFormResultFile(cmd$outFile, cmd$isMSMS, cmd$adduct, cmd$topMost)
    if (is.null(f))
        f <- data.table::data.table()
    else if (cmd$MSMode == "msms")
    {
        # note that even if MSMS data is available we may get MS only
        # formula in case no peaks could be explained
        f <- f[explainedPeaks > 0]
    }
    return(f)
}

GenFormMPTimeoutHandler <- function(cmd, retries)
{
    warning(paste("Formula calculation timed out for", cmd$group,
                  if (!is.null(cmd[["ana"]])) sprintf("(analysis '%s')", cmd$ana) else ""),
            call. = FALSE)
    return(FALSE)
}

GenFormMPPrepareHandler <- function(cmd)
{
    gfBin <- patRoon:::getGenFormBin()
    
    MSFile <- tempfile("MSPList", fileext = ".txt")
    outFile <- tempfile("formulas", fileext = ".txt")
    
    writePList <- function(f, pl) data.table::fwrite(pl[, c("mz", "intensity")], f, quote = FALSE,
                                                     sep = "\t", row.names = FALSE, col.names = FALSE)
    writePList(MSFile, cmd$MSPL)
    
    MSMSFile <- NULL
    if (!is.null(cmd[["MSMSPL"]]))
    {
        MSMSFile <- tempfile("MSMSPList", fileext = ".txt")
        writePList(MSMSFile, cmd$MSMSPL)
    }
    
    args <- c(sprintf("ms=%s", MSFile), sprintf("m=%f", cmd$PLMZ),
              sprintf("out=%s", outFile))
    if (cmd$MSMode != "ms" && cmd$isMSMS)
        args <- c(args, sprintf("msms=%s", MSMSFile), "analyze")
    
    cmd$args <- c(cmd$args, args)
    return(c(cmd, list(command = gfBin, MSFile = MSFile,
                       MSMSFile = MSMSFile, outFile = outFile)))
}

runGenForm <- function(mainArgs, groupPeakLists, annTable, MSMode, isolatePrec,
                       baseHash, setHash, gNames, adduct, topMost,
                       batchSize, timeout, ana)
{
    cmdQueue <- makeGenFormCmdQueue(groupPeakLists, annTable, baseHash, MSMode, isolatePrec, adduct, topMost,
                                    mainArgs, ana)
    
    ret <- list()
    if (length(cmdQueue) > 0)
    {
        ret <- executeMultiProcess(cmdQueue, finishHandler = patRoon:::GenFormMPFinishHandler,
                                   timeoutHandler = patRoon:::GenFormMPTimeoutHandler,
                                   prepareHandler = patRoon:::GenFormMPPrepareHandler,
                                   waitTimeout = 10, batchSize = batchSize, procTimeout = timeout,
                                   cacheName = "formulasGenForm", setHash = setHash)
    }
    
    return(pruneList(ret, checkZeroRows = TRUE))
}

getEmptyGFFragInfo <- function() data.table(mz = numeric(), ion_formula = character(), dbe = numeric(),
                                            ion_formula_mz = numeric(), error = numeric())

processGenFormMSMSResultFile <- function(file)
{
    formsMSMS <- data.table::fread(file, sep = "\t", fill = TRUE, header = FALSE) # fill should be set with strange formatting of file
    
    # first column: either formula precursor or fragment mass
    # second column: either dbe or formula fragment/neutral loss
    # third column: either precursor formula mz or dbe fragment
    # fourth column: either mass deviation precursor or fragment formula mz
    # fifth column: either MS score precursor or mass deviation fragment
    # sixth column: either MSMS score for precursor or NA
    # seventh column: either combined score for precursor or NA
    
    # first split MS and MSMS results
    formsMSMS[, isPrecursor := !is.na(formsMSMS[[6]])]
    formsMSMS[, precursorGroup := cumsum(isPrecursor)] # assign unique ID for all precursors+their fragments
    
    # generate fragInfos
    fragInfos <- formsMSMS[isPrecursor == FALSE]
    setnames(fragInfos, seq_len(5), c("mz", "ion_formula", "dbe", "ion_formula_mz", "error"))
    
    # apply Hill sorting
    fragInfos[, ion_formula := sapply(ion_formula, sortFormula)]
    
    # If multiple results exist for a mass then subsequent results after the first have empty mz column.
    # Fill in those empty masses
    for (fi in seq_len(nrow(fragInfos)))
    {
        if (!nzchar(fragInfos$mz[fi]))
            set(fragInfos, fi, "mz", fragInfos$mz[fi - 1])
    }
    
    fragInfos[, c("V6", "V7", "isPrecursor") := NULL]
    
    # set correct column types
    for (col in c("mz", "dbe", "ion_formula_mz", "error"))
        data.table::set(fragInfos, j = col, value = as.numeric(fragInfos[[col]]))
    data.table::set(fragInfos, j = "ion_formula", value = as.character(fragInfos$ion_formula))
    
    fragInfoList <- split(fragInfos, by = "precursorGroup", keep.by = FALSE) # NOTE: will be named by precursorGroup
    
    ret <- formsMSMS[isPrecursor == TRUE]
    data.table::setnames(ret, seq_len(7),
                         c("neutral_formula", "dbe", "ion_formula_mz", "error", "isoScore", "MSMSScore", "combMatch"))
    # initialize all with empty fragInfos
    ret[, fragInfo := list(list(getEmptyGFFragInfo()))]
    
    if (length(fragInfoList) > 0)
        ret[match(as.integer(names(fragInfoList)), precursorGroup), fragInfo := fragInfoList]
    
    ret[, c("isPrecursor", "precursorGroup") := NULL]
    
    return(ret)
}

processGenFormResultFile <- function(file, isMSMS, adduct, topMost)
{
    if (file.size(file) == 0)
        return(NULL)
    
    if (!isMSMS)
    {
        forms <- data.table::fread(file, header = FALSE)
        data.table::setnames(forms, c("neutral_formula", "dbe", "ion_formula_mz", "error", "isoScore"))
        forms[, fragInfo := list(list(patRoon:::getEmptyGFFragInfo()))]
    }
    else
        forms <- patRoon:::processGenFormMSMSResultFile(file)
    
    if (is.null(forms) || nrow(forms) == 0)
        return(NULL)
    
    forms <- patRoon:::rankFormulaTable(forms)
    
    # select topMost after ranking
    if (!is.null(topMost) && nrow(forms) > topMost)
        forms <- forms[seq_len(topMost)]

    forms[, neutral_formula := sapply(neutral_formula, sortFormula)] # GenForm doesn't seem to use Hill sorting
    
    forms <- addMiscFormulaInfo(forms, adduct)
    
    # set correct column types
    numCols <- intersect(c("error", "dbe", "isoScore", "MSMSScore", "combMatch"), names(forms))
    for (col in numCols)
        data.table::set(forms, j = col, value = as.numeric(forms[[col]]))
    data.table::set(forms, j = "neutral_formula", value = as.character(forms$neutral_formula))

    # set nice column order
    data.table::setcolorder(forms, c("neutral_formula", "ion_formula", "neutralMass", "ion_formula_mz", "error", "dbe",
                                     "isoScore", "explainedPeaks"))
    
    return(forms)
}

# UNDONE: fuzzy formulas

#' Generate formula with GenForm
#'
#' Uses \href{https://sourceforge.net/projects/genform/}{GenForm} to generate chemical formula candidates.
#'
#' @templateVar algo GenForm
#' @templateVar do generate formula candidates
#' @templateVar generic generateFormulas
#' @templateVar algoParam genform
#' @template algo_generator
#'
#' @details When MS/MS data is available it will be used to score candidate formulae by presence of 'fitting' fragments.
#'
#' @param relMzDev Maximum relative deviation between the measured and candidate formula \emph{m/z} values (in ppm).
#'   Sets the \option{ppm} command line option.
#' @param elements Elements to be considered for formulae calculation. This will heavily affects the number of
#'   candidates! Always try to work with a minimal set by excluding elements you don't expect. Sets the \option{el}
#'   command line option.
#' @param hetero Only consider formulae with at least one hetero atom. Sets the \option{het} commandline option.
#' @param oc Only consider organic formulae (\emph{i.e.} with at least one carbon atom). Sets the \option{oc}
#'   commandline option.
#' @param thrMS,thrMSMS,thrComb Sets the thresholds for the \command{GenForm} MS score (\code{isoScore}), MS/MS score
#'   (\code{MSMSScore}) and combined score (\code{combMatch}). Sets the \option{thms}/\option{thmsms}/\option{thcomb}
#'   command line options, respectively. Set to \code{NULL} for no threshold.
#' @param maxCandidates If this number of candidates are found then \command{GenForm} aborts any further formula
#'   calculations. The number of candidates is determined \emph{after} any formula filters, hence, the properties and
#'   'quality' of the candidates is influenced by options such as \code{oc} and \code{thrMS} arguments. Note that this
#'   is different than \code{topMost}, which selects the candidates after \command{GenForm} finished. Sets the
#'   \option{max} command line option. Set to \samp{0} or \code{Inf} for no maximum.
#' @param extraOpts An optional character vector with any other command line options that will be passed to
#'   \command{GenForm}. See the \verb{GenForm options} section for all available command line options.
#' @param MSMode Whether formulae should be generated only from MS data (\code{"ms"}), MS/MS data (\code{"msms"}) or
#'   both (\code{"both"}). Selecting \code{"both"} will fall back to formula calculation with only MS data in case no
#'   MS/MS data is available.
#' @param isolatePrec Settings used for isolation of precursor mass peaks and their isotopes. This isolation is highly
#'   important for accurate isotope scoring of candidates, as non-relevant mass peaks will dramatically decrease the
#'   score. The value of \code{isolatePrec} should either be a \code{list} with parameters (see the
#'   \code{\link[=filter,MSPeakLists-method]{filter method}} for \code{MSPeakLists} for more details), \code{TRUE} for
#'   default parameters or \code{FALSE} for no isolation (\emph{e.g.} when you already performed isolation with the
#'   \code{filter} method). The \code{z} parameter (charge) is automatically deduced from the adduct used for annotation
#'   (unless \code{isolatePrec=FALSE}), hence any custom \code{z} setting is ignored.
#' @param timeout Maximum time (in seconds) that a \command{GenForm} command is allowed to execute. If this time is
#'   exceeded a warning is emitted and the command is terminated. See the notes section for more information on the need
#'   of timeouts.
#' @param batchSize Maximum number of \command{GenForm} commands that should be run sequentially in each parallel
#'   process. Combining commands with short runtimes (such as \command{GenForm}) can significantly increase parallel
#'   performance. For more information see \code{\link{executeMultiProcess}}. Note that this is ignored if
#'   \option{patRoon.MP.method="future"}.
#'
#' @template adduct-arg
#' @templateVar algo genform
#' @template form_algo-args
#'
#' @inheritParams generateFormulas
#'
#' @inherit generateFormulas return
#'
#' @section GenForm options: Below is a list of options (generated by running \command{GenForm} without commandline
#'   options) which can be set by the \code{extraOpts} parameter.
#'
#' @eval paste0("@@section GenForm options: \\preformatted{", patRoon:::readAllFile(system.file("misc", "genform.txt",
#'   package = "patRoon")), "}")
#'
#' @templateVar what \code{generateFormulasGenForm}
#' @template uses-multiProc
#'
#' @section Parallelization: When \code{futures} are used for parallel processing (\code{patRoon.MP.method="future"}),
#'   calculations with \command{GenForm} are done with batch mode disabled (see \code{batchSize} argument), which
#'   generally limit overall performance.
#'
#' @note This function always sets the \option{exist} and \option{oei} \command{GenForm} command line options.
#'
#'   Formula calculation with \command{GenForm} may produce an excessive number of candidates for high \emph{m/z} values
#'   (\emph{e.g.} above 600) and/or many elemental combinations (set by \code{elements}). In this scenario formula
#'   calculation may need a very long time. Timeouts are used to avoid excessive computational times by terminating long
#'   running commands (set by the \code{timeout} argument).
#'
#' @references \insertRef{Meringer2011}{patRoon}
#'
#' @templateVar what generateFormulasGenForm
#' @template main-rd-method
#' @export
setMethod("generateFormulasGenForm", "featureGroups", function(fGroups, MSPeakLists,
                                                               specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
                                                               relMzDev = 5, adduct = NULL,
                                                               elements = "CHNOP", hetero = TRUE, oc = FALSE,
                                                               thrMS = NULL, thrMSMS = NULL, thrComb = NULL,
                                                               maxCandidates = Inf, extraOpts = NULL,
                                                               calculateFeatures = TRUE, featThreshold = 0,
                                                               featThresholdAnn = 0.75, absAlignMzDev = 0.002,
                                                               MSMode = "both", isolatePrec = TRUE, timeout = 120,
                                                               topMost = 50, batchSize = 8)
{
    if (is.infinite(maxCandidates))
        maxCandidates <- 0
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    assertSpecSimParams(specSimParams, add = ac)
    aapply(checkmate::assertNumber, . ~ relMzDev + timeout, lower = 0, finite = TRUE,
           fixed = list(add = ac))
    aapply(checkmate::assertString, . ~ elements, fixed = list(add = ac))
    aapply(checkmate::assertFlag, . ~ hetero + oc + calculateFeatures, fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ thrMS + thrMSMS + thrComb, null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertCount(maxCandidates, add = ac)
    aapply(checkmate::assertNumber, . ~ featThreshold + featThresholdAnn + absAlignMzDev, lower = 0, upper = 1,
           fixed = list(add = ac))
    checkmate::assertChoice(MSMode, c("ms", "msms", "both"), add = ac)
    checkmate::assertCharacter(extraOpts, null.ok = TRUE, add = ac)
    checkmate::assertCount(topMost, positive = TRUE, add = ac)
    checkmate::assertCount(batchSize, positive = TRUE, add = ac)
    
    if (!is.logical(isolatePrec))
        assertPListIsolatePrecParams(isolatePrec, add = ac)
    
    checkmate::reportAssertions(ac)
    
    if (length(fGroups) == 0)
        return(formulas(algorithm = "genform"))
    
    adduct <- checkAndToAdduct(adduct, fGroups)
    
    gInfo <- groupInfo(fGroups)
    anaInfo <- analysisInfo(fGroups)
    featIndex <- groupFeatIndex(fGroups)
    fTable <- featureTable(fGroups)
    gCount <- length(fGroups)
    gNames <- names(fGroups)
    annTable <- annotations(fGroups)
    
    MSPeakLists <- MSPeakLists[, gNames] # only do relevant
    
    mainArgs <- c("exist",
                  "oei",
                  "noref",
                  "dbe",
                  "cm",
                  sprintf("ppm=%f", relMzDev),
                  sprintf("el=%s", elements),
                  sprintf("max=%d", maxCandidates),
                  extraOpts)
    if (!is.null(thrMS))
        mainArgs <- c(mainArgs, sprintf("thms=%f", thrMS))
    if (!is.null(thrMSMS))
        mainArgs <- c(mainArgs, sprintf("thmsms=%f", thrMSMS))
    if (!is.null(thrComb))
        mainArgs <- c(mainArgs, sprintf("thcomb=%f", thrComb))
    if (hetero)
        mainArgs <- c(mainArgs, "het")
    if (oc)
        mainArgs <- c(mainArgs, "oc")
    
    formTable <- list()
    baseHash <- makeHash(mainArgs, MSMode, isolatePrec)
    setHash <- makeHash(fGroups, MSPeakLists, baseHash)
    
    # ana is optional and not used when only calculating group average formulas
    doGenForm <- function(groupPeakLists, ana)
    {
        if (!is.null(ana))
            printf("Loading all formulas for analysis '%s'...\n", ana)
        else
            printf("Loading all formulas...\n")
        
        ann <- if (nrow(annTable) > 0) annTable[match(names(groupPeakLists), group)] else annTable
        forms <- runGenForm(mainArgs, groupPeakLists, ann, MSMode, isolatePrec, baseHash,
                            setHash, gNames, adduct, topMost, batchSize, timeout, ana)
        
        printf("Loaded %d formulas for %d %s (%.2f%%).\n", countUniqueFormulas(forms), length(forms),
               if (!is.null(ana)) "features" else "feature groups",
               if (gCount == 0) 0 else length(forms) * 100 / gCount)
        
        return(forms)
    }
    
    if (calculateFeatures)
    {
        pLists <- peakLists(MSPeakLists)
        
        formTable <- sapply(seq_along(anaInfo$analysis), function(anai)
        {
            ana <- anaInfo$analysis[anai]
            if (is.null(pLists[[ana]]))
                return(NULL)
            
            pl <- pLists[[ana]]
            ftinds <- sapply(gNames, function(grp) featIndex[[grp]][anai])
            doFGroups <- gNames[ftinds != 0] # prune missing
            doFGroups <- intersect(doFGroups, names(pl))
            
            return(doGenForm(pl[doFGroups], ana))
        }, simplify = FALSE)
        
        names(formTable) <- anaInfo$analysis
        formTable <- pruneList(formTable, TRUE)

        if (length(formTable) > 0)
            groupFormulas <- generateGroupFormulasByConsensus(formTable, lapply(featIndex, function(x) sum(x > 0)),
                                                              featThreshold, featThresholdAnn, gNames)
        else
            groupFormulas <- list()
    }
    else
    {
        groupFormulas <- doGenForm(averagedPeakLists(MSPeakLists), NULL)
        formTable <- list()
    }
    
    groupFormulas <- setFormulaPLID(groupFormulas, MSPeakLists, absAlignMzDev)
    
    return(formulas(groupAnnotations = groupFormulas, featureFormulas = formTable, algorithm = "genform",
                    MSPeakLists = MSPeakLists, specSimParams = specSimParams))
})

#' @template featAnnSets-gen_args
#' @rdname generateFormulasGenForm
#' @export
setMethod("generateFormulasGenForm", "featureGroupsSet", function(fGroups, MSPeakLists,
                                                                  specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
                                                                  relMzDev = 5, adduct = NULL,
                                                                  ..., setThreshold = 0, setThresholdAnn = 0,
                                                                  setAvgSpecificScores = FALSE)
{
    generateFormulasSet(fGroups, MSPeakLists, specSimParams, adduct, generateFormulasGenForm, relMzDev = relMzDev, ...,
                        setThreshold = setThreshold, setThresholdAnn = setThresholdAnn,
                        setAvgSpecificScores = setAvgSpecificScores)
})
