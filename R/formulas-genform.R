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

makeGenFormCmdQueue <- function(groupPeakLists, baseHash, MSMode, isolatePrec, adduct, topMost,
                                mainArgs, ana)
{
    pruneList(sapply(names(groupPeakLists), function(grp)
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
        
        if (!is.logical(isolatePrec))
            plms <- isolatePrecInMSPeakList(plms, isolatePrec, negate = FALSE)
 
        hash <- makeHash(groupPeakLists[[grp]], baseHash)
               
        return(list(args = mainArgs, PLMZ = plmz, MSPL = plms, MSMSPL = groupPeakLists[[grp]][["MSMS"]],
                    hash = hash, group = grp, isMSMS = hasMSMS, adduct = adduct,
                    topMost = topMost, MSMode = MSMode, ana = ana))
    }, simplify = FALSE))
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
        f <- f[byMSMS == TRUE]
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

runGenForm <- function(mainArgs, groupPeakLists, MSMode, isolatePrec,
                       baseHash, setHash, gNames, adduct, topMost,
                       batchSize, timeout, ana)
{
    cmdQueue <- makeGenFormCmdQueue(groupPeakLists, baseHash, MSMode, isolatePrec, adduct, topMost,
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

    fMS <- formsMSMS[isPrecursor == TRUE, ]
    data.table::setnames(fMS, 1:7, c("neutral_formula", "dbe", "formula_mz", "error", "isoScore", "MSMSScore", "combMatch"))
    fMS[, isPrecursor := NULL]
    data.table::setkey(fMS, "precursorGroup")

    fMSMS <- formsMSMS[isPrecursor == FALSE, ]

    setnames(fMSMS, 1:5, c("frag_mz", "frag_formula", "frag_dbe", "frag_formula_mz", "frag_error"))

    # If multiple results exist for a mass then subsequent results after the first have empty mz column.
    # Fill in those empty masses
    for (fi in seq_len(nrow(fMSMS)))
    {
        if (!nzchar(fMSMS$frag_mz[fi]))
            set(fMSMS, fi, "frag_mz", fMSMS$frag_mz[fi - 1])
    }

    fMSMS[, c("V6", "V7", "isPrecursor") := NULL]
    data.table::setkey(fMSMS, "precursorGroup")

    formsMSMS <- merge(fMS, fMSMS, all.x = TRUE) # re-join (but ensure that MS only formulae in fMS are kept)
    formsMSMS[, byMSMS := !is.na(frag_formula)]
    formsMSMS[, precursorGroup := NULL]

    # clear out MSMS scores for formulae w/out MSMS explanations
    # UNDONE: do we need to discern candidates w/out MSMS data and w/ MSMS data but no explanations?
    formsMSMS[byMSMS == FALSE, c("MSMSScore", "combMatch") := NA_real_]

    return(formsMSMS)
}

processGenFormResultFile <- function(file, isMSMS, adduct, topMost)
{
    if (file.size(file) == 0)
        return(NULL)

    if (!isMSMS)
    {
        forms <- data.table::fread(file, header = FALSE)
        data.table::setnames(forms, c("neutral_formula", "dbe", "formula_mz", "error", "isoScore"))
        forms[, byMSMS := FALSE]
    }
    else
        forms <- patRoon:::processGenFormMSMSResultFile(file)

    if (is.null(forms) || nrow(forms) == 0)
        return(NULL)

    forms <- patRoon:::rankFormulaTable(forms)
    
    # select topMost after ranking
    if (!is.null(topMost) && data.table::uniqueN(forms, by = "neutral_formula") > topMost)
    {
        forms[, unFormID := .GRP, by = "neutral_formula"]
        forms <- forms[unFormID <= topMost][, unFormID := NULL]
    }
    
    forms[, neutral_formula := sapply(neutral_formula, sortFormula)] # GenForm doesn't seem to use Hill sorting
    forms[, formula := calculateIonFormula(neutral_formula, adduct)]

    # set correct column types
    numCols <- intersect(c("error", "dbe", "isoScore", "frag_mz", "frag_error",
                           "frag_dbe", "MSMSScore", "combMatch"), names(forms))
    for (col in numCols)
        data.table::set(forms, j = col, value = as.numeric(forms[[col]]))
    
    chrCols <- intersect(c("formula", "frag_formula", "neutral_formula"), names(forms))
    for (col in chrCols)
        data.table::set(forms, j = col, value = as.character(forms[[col]]))
    
    if (!is.null(forms[["frag_formula"]]))
    {
        forms[byMSMS == TRUE, neutral_loss := as.character(Vectorize(subtractFormula)(formula, frag_formula))]
        forms[byMSMS == TRUE, frag_formula := Vectorize(sortFormula)(frag_formula)]
    }

    # set nice column order
    data.table::setcolorder(forms, c("neutral_formula", "formula", "formula_mz", "error", "dbe", "isoScore", "byMSMS"))

    return(forms)
}

# UNDONE: fuzzy formulas

#' @details \code{generateFormulasGenForm} uses
#'   \href{https://sourceforge.net/projects/genform/}{GenForm} to generate
#'   chemical formulae. When MS/MS data is available it will be used to score
#'   candidate formulae by presence of 'fitting' fragments.
#'
#' @section GenForm options: Below is a list of options (generated by running
#'   \command{GenForm} without commandline options) which can be set by the
#'   \code{extraOpts} parameter.
#'
#' @eval paste0("@@section GenForm options: \\preformatted{",
#'   patRoon:::readAllFile(system.file("misc", "genform.txt", package =
#'   "patRoon")), "}")
#'
#' @param hetero Only consider formulae with at least one hetero atom. Sets the
#'   \option{het} commandline option.
#' @param oc Only consider organic formulae (\emph{i.e.} with at least one
#'   carbon atom). Sets the \option{oc} commandline option.
#' @param isolatePrec Settings used for isolation of precursor mass peaks and
#'   their isotopes. This isolation is highly important for accurate isotope
#'   scoring of candidates, as non-relevant mass peaks will dramatically
#'   decrease the score. The value of \code{isolatePrec} should either be a
#'   \code{list} with parameters (see the
#'   \code{\link[=filter,MSPeakLists-method]{filter method}} for
#'   \code{MSPeakLists} for more details), \code{TRUE} for default parameters
#'   (the \code{z} parameter is automatically deduced from the \code{adduct}
#'   argument) or \code{FALSE} for no isolation (\emph{e.g.} when you already
#'   performed isolation with the \code{filter} method).
#' @param timeout Maximum time (in seconds) that a \command{GenForm} command is
#'   allowed to execute. If this time is exceeded a warning is emitted and the
#'   command is terminated. See the notes section for more information on the
#'   need of timeouts.
#' @param batchSize Maximum number of \command{GenForm} commands that should be
#'   run sequentially in each parallel process. Combining commands with short
#'   runtimes (such as \command{GenForm}) can significantly increase parallel
#'   performance. For more information see \code{\link{executeMultiProcess}}.
#'   Note that this is ignored if \option{patRoon.MP.method="future"}.
#' 
#' @section Parallelization: When \code{futures} are used for parallel processing
#'   (\code{patRoon.MP.method="future"}), calculations with \command{GenForm}
#'   are done with batch mode is disabled (see \code{batchSize} argument), which
#'   generally limit overall performance.
#'   
#' @note \code{generateFormulasGenForm} always sets the \option{exist} and
#'   \option{oei} \command{GenForm} commandline options.
#'
#'   Formula calculation with \command{GenForm} may produce an excessive number
#'   of candidates for high \emph{m/z} values (\emph{e.g.} above 600) and/or
#'   many elemental combinations (set by \code{elements}). In this scenario
#'   formula calculation may need a very long time. Timeouts are used to avoid
#'   excessive computational times by terminating long running commands (set by
#'   the \code{timeout} argument).
#'
#' @references \insertRef{Meringer2011}{patRoon}
#'
#' @rdname formula-generation
#' @export
setMethod("generateFormulasGenForm", "featureGroups", function(fGroups, MSPeakLists, relMzDev = 5, adduct = "[M+H]+",
                                                               elements = "CHNOP", hetero = TRUE, oc = FALSE, extraOpts = NULL,
                                                               calculateFeatures = TRUE, featThreshold = 0.75, MSMode = "both",
                                                               isolatePrec = TRUE, timeout = 120, topMost = 50,
                                                               batchSize = 8)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    aapply(checkmate::assertNumber, . ~ relMzDev + timeout, lower = 0, finite = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertString, . ~ elements, fixed = list(add = ac))
    aapply(checkmate::assertFlag, . ~ hetero + oc + calculateFeatures, fixed = list(add = ac))
    checkmate::assertNumber(featThreshold, lower = 0, finite = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertChoice(MSMode, c("ms", "msms", "both"), add = ac)
    checkmate::assertCharacter(extraOpts, null.ok = TRUE, add = ac)
    checkmate::assertCount(topMost, positive = TRUE, add = ac)
    checkmate::assertCount(batchSize, positive = TRUE, add = ac)

    if (!is.logical(isolatePrec))
         assertPListIsolatePrecParams(isolatePrec, add = ac)

    checkmate::reportAssertions(ac)

    adduct <- checkAndToAdduct(adduct)

    if (is.logical(isolatePrec) && isolatePrec)
        isolatePrec <- getDefIsolatePrecParams(z = abs(adduct@charge))
    
    gInfo <- groupInfo(fGroups)
    anaInfo <- analysisInfo(fGroups)
    featIndex <- groupFeatIndex(fGroups)
    fTable <- featureTable(fGroups)
    gCount <- length(fGroups)
    gNames <- names(fGroups)

    MSPeakLists <- MSPeakLists[, gNames] # only do relevant
    
    mainArgs <- c("exist",
                  "oei",
                  "noref",
                  "dbe",
                  "cm",
                  sprintf("ion=%s", as.character(adduct, format = "genform")),
                  sprintf("ppm=%f", relMzDev),
                  sprintf("el=%s", elements),
                  extraOpts)

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

        forms <- runGenForm(mainArgs, groupPeakLists, MSMode, isolatePrec,
                            baseHash, setHash, gNames, adduct, topMost, batchSize, timeout, ana)

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
            
            ftinds <- sapply(gNames, function(grp) featIndex[[grp]][anai])
            doFGroups <- gNames[ftinds != 0] # prune missing
            
            return(doGenForm(pLists[[ana]][doFGroups], ana))
        }, simplify = FALSE)
        
        names(formTable) <- anaInfo$analysis
        formTable <- pruneList(formTable, TRUE)

        if (length(formTable) > 0)
            groupFormulas <- generateGroupFormulasByConsensus(formTable, featThreshold, gNames, "analysis", "anaCoverage")
        else
            groupFormulas <- list()
    }
    else
    {
        groupFormulas <- doGenForm(averagedPeakLists(MSPeakLists), NULL)
        formTable <- list()
    }

    return(formulas(formulas = groupFormulas, featureFormulas = formTable, algorithm = "genform"))
})

setMethod("generateFormulasGenForm", "featureGroupsSet", function(fGroups, MSPeakLists, ..., setThreshold = 0.75)
{
    setArgs <- assertAndGetMSPLSetsArgs(fGroups, MSPeakLists)
    generateFormulasSet(fGroups, generateFormulasGenForm, ..., setArgs = setArgs, setThreshold = setThreshold)
})
