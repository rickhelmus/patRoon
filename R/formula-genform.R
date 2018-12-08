#' @include main.R
#' @include mspeaklists.R
#' @include formulas.R
NULL

getGenFormBin <- function()
{
    gbin <- if (Sys.info()[["sysname"]] == "Linux") "GenForm" else "GenForm.exe"
    pOpt <- getOption("patRoon.path.GenForm")
    if (!is.null(pOpt) && nzchar(pOpt))
        ret <- file.path(pOpt, gbin)
    else
    {
        ret <- system.file("bin", Sys.getenv("R_ARCH"), gbin, package = "patRoon")

        if (!nzchar(ret) &&
            ((requireNamespace("pkgload", quietly = TRUE) && !is.null(pkgload::dev_meta("patRoon"))) ||
             (requireNamespace("devtools", quietly = TRUE) && !is.null(devtools::dev_meta("patRoon")))))
            ret <- normalizePath(file.path(system.file(".", package = "patRoon"), "..", "src", gbin))
    }

    if (!file.exists(ret))
        stop(sprintf("GenForm binary does not exist: %s", ret))

    return(ret)
}

initGenFormFiles <- function(MSPList, MSMSPList)
{
    writePList <- function(f, pl) fwrite(pl[, c("mz", "intensity")], f, quote = FALSE,
                                         sep = "\t", row.names = FALSE, col.names = FALSE)
    
    if (is.null(MSPList))
        return(NULL)
    
    MSFile <- tempfile("MSPList", fileext = ".txt")
    writePList(MSFile, MSPList)
    
    MSMSFile <- NULL
    if (!is.null(MSMSPList))
    {
        MSMSFile <- tempfile("MSMSPList", fileext = ".txt")
        writePList(MSMSFile, MSMSPList)
    }
    
    outFile <- tempfile("formulas", fileext = ".txt")
    
    return(list(MSFile = MSFile, MSMSFile = MSMSFile, outFile = outFile))
}

makeGenFormFilesNames <- function()
{
    MSFile <- tempfile("MSPList", fileext = ".txt")
    MSMSFile <- tempfile("MSMSPList", fileext = ".txt")
    outFile <- tempfile("formulas", fileext = ".txt")
    return(list(MSFile = MSFile, MSMSFile = MSMSFile, outFile = outFile))
}

writeGenFormFiles <- function(MSPList, MSMSPList, workFiles)
{
    writePList <- function(f, pl) fwrite(pl[, c("mz", "intensity")], f, quote = FALSE,
                                         sep = "\t", row.names = FALSE, col.names = FALSE)
    
    # file may have been already created when calculating MS only formulas    
    if (!file.exists(workFiles$MSFile))
        writePList(workFiles$MSFile, MSPList)
    
    if (!is.null(MSMSPList))
        writePList(workFiles$MSMSFile, MSMSPList)
}

makeGenFormCmdQueue <- function(gfBin, mainArgs, featMZs, groupPeakLists, workFiles, hashes, doMSMS)
{
    pruneList(sapply(names(groupPeakLists), function(grp)
    {
        if (is.null(groupPeakLists[[grp]]) || is.null(groupPeakLists[[grp]][["MS"]]))
            return(NULL) # MS or MSMS spectrum probably filtered away
        
        if (doMSMS && is.null(groupPeakLists[[grp]][["MSMS"]]))
            return(NULL)
        
        plmz <- getMZFromMSPeakList(featMZs[grp], groupPeakLists[[grp]][["MS"]])
        
        args <- c(mainArgs, sprintf("ms=%s", workFiles[[grp]]$MSFile), sprintf("m=%f", plmz),
                  sprintf("out=%s", workFiles[[grp]]$outFile))
        
        if (doMSMS)
        {
            if (is.null(workFiles[[grp]]$MSMSFile))
                return(NULL)
            args <- c(args, sprintf("msms=%s", workFiles[[grp]]$MSMSFile), "analyze")
        }
        
        return(list(command = gfBin, args = args, outFile = workFiles[[grp]]$outFile, hash = hashes[grp]))
    }, simplify = FALSE))
}

runGenForm <- function(gfBin, mainArgs, featMZs, groupPeakLists, doMSMS, hashes, cachedSet,
                       workFiles, gNames, adduct, maxProcAmount, maxCmdsPerProc)
{
    cacheDB <- openCacheDBScope()
    
    cachedResults <- pruneList(sapply(hashes, function(h)
    {
        forms <- NULL
        if (!is.null(cachedSet))
            forms <- cachedSet[[h]]
        if (is.null(forms))
            forms <- loadCacheData("formulasGenForm", h, cacheDB)
        return(forms)
    }, simplify = FALSE))
    
    # skip cached results
    doGNames <- setdiff(names(featMZs), names(cachedResults))
    if (length(doGNames) > 0)
    {
        # filter out fgroups without peaklists
        doGNames <- doGNames[sapply(groupPeakLists[doGNames],
                                    function(pl) !is.null(pl[["MS"]]) && (!doMSMS || !is.null(pl[["MSMS"]])))]
    }
    
    for (grp in doGNames)
        writeGenFormFiles(groupPeakLists[[grp]][["MS"]],
                          if (!doMSMS) groupPeakLists[[grp]][["MSMS"]] else NULL,
                          workFiles[[grp]])
    
    cmdQueue <- makeGenFormCmdQueue(gfBin, mainArgs, featMZs[doGNames], groupPeakLists[doGNames],
                                    workFiles[doGNames], hashes, doMSMS)

    ret <- list()
    if (length(cmdQueue) > 0)
    {
        ret <- executeMultiProcess(cmdQueue, function(cmd)
        {
            f <- processGenFormResultFile(cmd$outFile, doMSMS, adduct)
            if (is.null(f))
                f <- data.table()
            saveCacheData("formulasGenForm", f, cmd$hash, cacheDB)
            return(f)
        }, maxProcAmount = maxProcAmount, waitTimeout = 10, maxCmdsPerProc = maxCmdsPerProc)
    }

    if (length(cachedResults) > 0)
    {
        ret <- c(ret, cachedResults)
        ret <- ret[intersect(gNames, names(ret))] # re-order
    }
    
    return(pruneList(ret, checkZeroRows = TRUE))
}

processGenFormMSMSResultFile <- function(file)
{
    formsMSMS <- fread(file, sep = "\t", fill = TRUE, header = FALSE) # fill should be set with strange formatting of file

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

    if (!any(!formsMSMS$isPrecursor)) # no fragments?
        return(NULL)

    fMS <- formsMSMS[isPrecursor == TRUE, ]
    setnames(fMS, 1:7, c("neutral_formula", "dbe", "formula_mz", "error", "MS_match", "MSMS_match", "comb_match"))
    fMS[, isPrecursor := NULL]
    setkey(fMS, "precursorGroup")

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
    setkey(fMSMS, "precursorGroup")

    formsMSMS <- fMS[fMSMS] # re-join
    formsMSMS[, byMSMS := TRUE]
    formsMSMS[, precursorGroup := NULL]

    formsMSMS[, frag_formula := Vectorize(sortFormula)(frag_formula)] # GenForm doesn't seem use Hill sorting

    return(formsMSMS)
}

processGenFormResultFile <- function(file, isMSMS, adduct)
{
    if (file.size(file) == 0)
        return(NULL)

    if (!isMSMS)
    {
        forms <- fread(file, header = FALSE)
        setnames(forms, c("neutral_formula", "dbe", "formula_mz", "error", "MS_match"))
        forms[, neutral_formula := Vectorize(sortFormula)(neutral_formula)] # GenForm doesn't seem use Hill sorting
        forms[, byMSMS := FALSE]
    }
    else
        forms <- processGenFormMSMSResultFile(file)

    if (is.null(forms))
        return(NULL)

    forms[, formula := calculateIonFormula(neutral_formula, adduct)]

    # set correct column types
    numCols <- c("error", "dbe", "MS_match", "frag_mz", "frag_error", "frag_dbe", "MSMS_match", "comb_match")
    for (col in numCols)
        set(forms, j = col, value = as.numeric(forms[[col]]))

    if (isMSMS)
    {
        forms[, neutral_loss := as.character(Vectorize(subtractFormula)(formula, frag_formula))]
        forms <- forms[nzchar(neutral_loss)] # remove fragments that are precursors
    }

    return(if (nrow(forms) > 0) forms else NULL)
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
#'   \Sexpr[results=verbatim,echo=FALSE,stage=build]{cat(patRoon:::readAllFile(system.file("misc",
#'   "genform.txt", package = "patRoon")))}
#'
#' @param hetero Only consider formulae with at least one hetero atom. Sets the
#'   \option{het} commandline option.
#' @param extraOpts An optional character vector with any other commandline
#'   options that will be passed to \command{GenForm}. See the \verb{GenForm
#'   options} section for all available commandline options.
#' @param maxCmdsPerProc Maximum number of commands that should be combined for
#'   each executed process. Combining commands with short runtimes (such as
#'   \command{GenForm}) can significantly increase parallel performance.
#'
#' @note \code{generateFormulasGenForm} always sets the \option{exist} and
#'   \option{oei} \command{GenForm} commandline options.
#'
#' @references \insertRef{Meringer2011}{patRoon}
#'
#' @rdname formula-generation
#' @export
generateFormulasGenForm <- function(fGroups, MSPeakLists, maxMzDev = 5, adduct = "M+H",
                                    elements = "CHNOP", hetero = TRUE, extraOpts = NULL,
                                    calculateBy = "feature", formFeatThreshold = 0.75, MSMode = "both",
                                    maxProcAmount = getOption("patRoon.maxProcAmount"), maxCmdsPerProc = 25)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertNumber(maxMzDev, lower = 0, finite = TRUE, add = ac)
    aapply(checkmate::assertString, . ~ adduct + elements, fixed = list(add = ac))
    checkmate::assertFlag(hetero, add = ac)
    checkmate::assertChoice(calculateBy, c("feature", "group"), add = ac)
    checkmate::assertChoice(MSMode, c("ms", "msms", "both"), add = ac)
    checkmate::assertList(extraOpts, any.missing = FALSE, names = "unique", null.ok = TRUE, add = ac)
    aapply(checkmate::assertCount, . ~ maxProcAmount + maxCmdsPerProc, positive = TRUE, fixed = list(add = ac))
    checkmate::reportAssertions(ac)

    gInfo <- groupInfo(fGroups)
    anaInfo <- analysisInfo(fGroups)
    featIndex <- groupFeatIndex(fGroups)
    fTable <- featureTable(fGroups)
    gCount <- length(fGroups)
    gNames <- names(fGroups)
    pLists <- peakLists(MSPeakLists)

    mainArgs <- c("exist",
                  "oei",
                  "noref",
                  "dbe",
                  "cm",
                  sprintf("ion=%s", adduct),
                  sprintf("ppm=%f", maxMzDev),
                  sprintf("el=%s", elements),
                  extraOpts)

    if (hetero)
        mainArgs <- c(mainArgs, "het")

    formTable <- list()
    cacheDB <- openCacheDBScope() # open manually so caching code doesn't need to on each R/W access
    baseHash <- makeHash(maxMzDev, adduct, elements, hetero, extraOpts)
    setHash <- makeHash(fGroups, MSMode, baseHash)
    cachedSet <- loadCacheSet("formulasGenForm", setHash, cacheDB)
    formHashes <- character(0)
    gfBin <- getGenFormBin()
    
    # ana is optional and not used when only calculating group average formulas
    doGenForm <- function(featMZs, groupPeakLists, ana)
    {
        doGNames <- names(featMZs)
        
        MSHashes <- sapply(doGNames, function(grp) makeHash(featMZs[grp], groupPeakLists[[grp]][["MS"]], baseHash))
        if (MSMode != "ms")
            MSMSHashes <- sapply(doGNames, function(grp) makeHash(MSHashes[grp], groupPeakLists[[grp]][["MSMS"]]))
        else
            MSMSHashes <- NULL
        
        formHashes <<- c(formHashes, MSHashes, MSMSHashes)
        
        workFiles <- sapply(doGNames, function(g) makeGenFormFilesNames(), simplify = FALSE)
        MSForms <- list()
        
        if (!is.null(ana))
        {
            startMsg <- sprintf("Loading all %%s formulas for analysis '%s'...\n", ana)
            endMsg <- "Loaded %d %s formulas for %d features (%.2f%%).\n"
        }
        else
        {
            startMsg <- "Loading all %s formulas...\n"
            endMsg <- "Loaded %d %s formulas for %d feature groups (%.2f%%).\n"
        }
        
        if (MSMode != "msms")
        {
            printf(startMsg, "MS")
            MSForms <- runGenForm(gfBin, mainArgs, featMZs, groupPeakLists, FALSE, MSHashes,
                                  cachedSet, workFiles, gNames, adduct, maxProcAmount, maxCmdsPerProc)

            printf(endMsg, sum(unlist(sapply(MSForms, nrow))), "MS", length(MSForms),
                   if (gCount == 0) 0 else length(MSForms) * 100 / gCount)
        }
        
        if (MSMode != "ms")
        {
            printf(startMsg, "MS/MS")
            
            MSMSForms <- runGenForm(gfBin, mainArgs, featMZs, groupPeakLists[doGNames], TRUE, MSMSHashes,
                                    cachedSet, workFiles, gNames, adduct, maxProcAmount, maxCmdsPerProc)
            
            printf(endMsg, sum(unlist(sapply(MSMSForms, nrow))), "MS/MS", length(MSMSForms),
                   if (gCount == 0) 0 else length(MSMSForms) * 100 / gCount)
            
            MSForms <- sapply(union(names(MSForms), names(MSMSForms)), function(grp)
            {
                ret <- MSForms[[grp]]
                if (!is.null(MSMSForms[[grp]]))
                {
                    if (is.null(ret) || nrow(ret) == 0)
                        ret <- MSMSForms[[grp]]
                    else
                    {
                        # remove any MS formulas also present in MSMS data
                        ret <- ret[!formula %in% MSMSForms[[grp]][["formula"]]]
                        
                        # merge
                        ret <- rbind(ret, MSMSForms[[grp]], fill = TRUE)
                    }
                }
                return(ret)
            }, simplify = FALSE)
        }
        
        return(MSForms)
    }
    
    if (calculateBy == "feature")
    {
        pLists <- peakLists(MSPeakLists)
        
        formTable <- sapply(seq_along(anaInfo$analysis), function(anai)
        {
            ana <- anaInfo$analysis[anai]
            
            if (is.null(pLists[[ana]]))
                return(NULL)
            
            ftinds <- sapply(gNames, function(grp) featIndex[[grp]][anai])
            ftinds <- ftinds[ftinds != 0] # prune missing
            featMZs <- sapply(ftinds, function(fti) fTable[[ana]][["mz"]][fti])
            
            return(doGenForm(featMZs, pLists[[ana]], ana))
        }, simplify = FALSE)
        
        names(formTable) <- anaInfo$analysis
        formTable <- pruneList(formTable, TRUE)
        
        if (is.null(cachedSet))
            saveCacheSet("formulasGenForm", formHashes, setHash, cacheDB)
        
        if (length(formTable) > 0)
            groupFormulas <- generateGroupFormulasByConsensus(formTable, formFeatThreshold)
        else
            groupFormulas <- list()
        
    }
    else
    {
        featMZs <- setNames(gInfo[, "mzs"], gNames)
        groupFormulas <- doGenForm(featMZs, averagedPeakLists(MSPeakLists), NULL)
        formTable <- list()
    }
    
    return(formulas(formulas = formTable, groupFormulas = groupFormulas, algorithm = "GenForm"))
}



generateFormulasGenFormOld <- function(fGroups, MSPeakLists, maxMzDev = 5, adduct = "M+H",
                                       elements = "CHNOP", hetero = TRUE, MSMode = "both", extraOpts = NULL,
                                       maxProcAmount = getOption("patRoon.maxProcAmount"), maxCmdsPerProc = 25)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertNumber(maxMzDev, lower = 0, finite = TRUE, add = ac)
    aapply(checkmate::assertString, . ~ adduct + elements, fixed = list(add = ac))
    checkmate::assertFlag(hetero, add = ac)
    checkmate::assertChoice(MSMode, c("ms", "msms", "both"), add = ac)
    checkmate::assertList(extraOpts, any.missing = FALSE, names = "unique", null.ok = TRUE, add = ac)
    aapply(checkmate::assertCount, . ~ maxProcAmount + maxCmdsPerProc, positive = TRUE, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    gInfo <- groupInfo(fGroups)
    anaInfo <- analysisInfo(fGroups)
    ftind <- groupFeatIndex(fGroups)
    fTable <- featureTable(fGroups)
    gCount <- length(fGroups)
    gNames <- names(fGroups)
    
    mainArgs <- c("exist",
                  "oei",
                  "noref",
                  "dbe",
                  "cm",
                  sprintf("ion=%s", adduct),
                  sprintf("ppm=%f", maxMzDev),
                  sprintf("el=%s", elements),
                  extraOpts)
    
    if (hetero)
        mainArgs <- c(mainArgs, "het")
    
    formTable <- list()
    cacheDB <- openCacheDBScope() # open manually so caching code doesn't need to on each R/W access
    setHash <- makeHash(fGroups, maxMzDev, adduct, elements, hetero, MSMode, extraOpts)
    cachedSet <- loadCacheSet("formulasGenForm", setHash, cacheDB)
    formHashes <- character(0)
    gfBin <- getGenFormBin()
    
    runGenForm <- function(anai, cmdData, doMSMS)
    {
        cmdQueue <- sapply(cmdData, function(cd)
        {
            args <- c(mainArgs, sprintf("ms=%s", cd$MSFile), sprintf("m=%f", cd$precursorMz), sprintf("out=%s", cd$outFile))
            
            if (doMSMS)
            {
                if (is.null(cd$MSMSFile))
                    return(NULL)
                args <- c(args, sprintf("msms=%s", cd$MSMSFile), "analyze")
            }
            
            return(list(command = gfBin, args = args, outFile = cd$outFile,
                        hash = if (doMSMS) cd$MSMSHash else cd$MSHash))
        }, simplify = FALSE)
        cmdQueue <- cmdQueue[!sapply(cmdQueue, is.null)]
        
        formHashes <<- c(formHashes, sapply(cmdQueue, function(cmd) cmd$hash))
        
        cachedResults <- sapply(cmdQueue, function(cmd)
        {
            forms <- NULL
            if (!is.null(cachedSet))
                forms <- cachedSet[[cmd$hash]]
            if (is.null(forms))
                forms <- loadCacheData("formulasGenForm", cmd$hash, cacheDB)
            return(forms)
            
        }, simplify = FALSE)
        cachedResults <- cachedResults[!sapply(cachedResults, is.null)]
        
        cmdQueue <- cmdQueue[setdiff(names(cmdQueue), names(cachedResults))] # remove cached results
        
        ret <- list()
        if (length(cmdQueue) > 0)
        {
            ret <- executeMultiProcess(cmdQueue, function(cmd)
            {
                f <- processGenFormResultFile(cmd$outFile, doMSMS, adduct)
                if (is.null(f))
                    f <- data.table()
                saveCacheData("formulasGenForm", f, cmd$hash, cacheDB)
                return(f)
            }, maxProcAmount = maxProcAmount, waitTimeout = 10, maxCmdsPerProc = maxCmdsPerProc)
        }
        
        if (length(cachedResults) > 0)
        {
            ret <- c(ret, cachedResults)
            ret <- ret[intersect(gNames, names(ret))] # re-order
        }
        
        ngrp <- sum(unlist(sapply(ret, function(ft) nrow(ft) > 0)))
        printf("Loaded %d %s formulas from %d features (%.2f%%).\n", sum(unlist(sapply(ret, nrow))),
               if (doMSMS) "MS/MS" else "MS", ngrp, if (gCount == 0) 0 else ngrp * 100 / gCount)
        
        return(ret)
    }
    
    for (anai in seq_along(anaInfo$analysis))
    {
        ana <- anaInfo$analysis[anai]
        
        cmdData <- sapply(gNames, function(grp)
        {
            if (ftind[[grp]][anai] == 0)
                return(NULL) # feature not present
            
            ana <- anaInfo$analysis[anai]
            
            if (is.null(pLists[[ana]][[grp]]) || is.null(pLists[[ana]][[grp]][["MS"]]))
                return(NULL) # MS or MSMS spectrum probably filtered away
            
            if (MSMode == "msms" && is.null(pLists[[ana]][[grp]][["MSMS"]]))
                return(NULL)
            
            ftmz <- fTable[[ana]][["mz"]][ftind[[grp]][anai]]
            plmz <- getMZFromMSPeakList(ftmz, pLists[[ana]][[grp]][["MS"]])
            
            MSHash <- makeHash(ftmz, pLists[[ana]][[grp]][["MS"]], maxMzDev, adduct, elements, hetero, extraOpts)
            
            if (MSMode != "ms")
                MSMSHash <- makeHash(MSHash, pLists[[ana]][[grp]][["MSMS"]])
            else
                MSMSHash <- NULL
            
            # UNDONE: don't make new files if cached
            files <- initGenFormFiles(pLists[[ana]][[grp]][["MS"]],
                                      if (MSMode != "ms") pLists[[ana]][[grp]][["MSMS"]] else NULL)
            
            return(c(list(precursorMz = plmz, MSHash = MSHash, MSMSHash = MSMSHash), files))
        }, simplify = FALSE)
        cmdData <- cmdData[!sapply(cmdData, is.null)]
        
        formTable[[ana]] <- list()
        
        MSForms <- list()
        if (MSMode != "msms")
        {
            printf("Loading all MS formulas for analysis '%s'...\n", ana)
            MSForms <- runGenForm(anai, cmdData, FALSE)
        }
        
        if (MSMode != "ms")
        {
            printf("Loading all MS/MS formulas for analysis '%s'...\n", ana)
            MSMSForms <- runGenForm(anai, cmdData, TRUE)
            
            MSForms <- sapply(union(names(MSForms), names(MSMSForms)), function(grp)
            {
                ret <- MSForms[[grp]]
                if (!is.null(MSMSForms[[grp]]))
                {
                    if (is.null(ret) || nrow(ret) == 0)
                        ret <- MSMSForms[[grp]]
                    else
                    {
                        # remove any MS formulas also present in MSMS data
                        ret <- ret[!formula %in% MSMSForms[[grp]][["formula"]]]
                        
                        # merge
                        ret <- rbind(ret, MSMSForms[[grp]], fill = TRUE)
                    }
                }
                return(ret)
            }, simplify = FALSE)
        }
        
        formTable[[ana]] <- MSForms
    }
    
    if (is.null(cachedSet))
        saveCacheSet("formulasGenForm", formHashes, setHash, cacheDB)
    
    formTable <- pruneList(sapply(formTable, pruneList, checkZeroRows = TRUE, simplify = FALSE), TRUE)
    
    if (length(formTable) > 0)
        groupFormulas <- generateGroupFormulasByConsensus(formTable, 0.75) # UNDONE
    
    return(formulas(formulas = formTable, groupFormulas = groupFormulas, algorithm = "GenForm"))
}
