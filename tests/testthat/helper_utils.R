testWithSets <- function() T # UNDONE: check environment variable or something

getTestDataPathGeneric <- function() "test_data"

testFile <- function(f, ..., text = FALSE) file.path(getTestDataPath(), paste0(f, ..., if (!text) ".Rds" else ".txt", collapse = ""))
getTestFGroups <- function(anaInfo = getTestAnaInfo(), ...) groupFeatures(getTestFeatures(anaInfo, ...), "openms")
getEmptyFeatures <- function(anaInfo = getTestAnaInfo(), ...) getTestFeatures(anaInfo, noiseThrInt = 1E9, ...)
getEmptyTestFGroups <- function(anaInfo = getTestAnaInfo()) getTestFGroups(anaInfo)[, "none"]

getFormFGroups <- function()
{
    # lower intensity threshold a bit to get benzotriazole in both polarities with sets
    fGroups <- getTestFGroups(getTestAnaInfoAnn(), noiseThrInt = 2E4)
    # convert to screening results to simplify things a bit
    # focus on some +/- suspects and one with a non C fragment
    return(doScreen(fGroups, patRoonData::suspectsPos[c(1:7, 33, 34, 36), ], onlyHits = TRUE))
    
}

getMFTestDBPath <- function() file.path(getTestDataPath(), "test-mf-db.csv")
getCompFGroups <- function()
{
    fGroups <- doScreen(getTestFGroups(getTestAnaInfoAnn(), noiseThrInt = 1E4), patRoonData::suspectsPos, onlyHits = TRUE)
    # just focus on few targets also present in MF test DB
    return(fGroups[, suspects = fread(getMFTestDBPath())$Name])
}

callMF <- function(fGroups, plists, scoreTypes = "fragScore", db = getMFTestDBPath(), to = 300, ...)
{
    doGenComps(fGroups, plists, "metfrag", timeout = to, database = "csv", scoreTypes = scoreTypes,
               extraOpts = list(LocalDatabasePath = db), ...)
}

getDAAnaInfoAll <- function(pat = NULL)
{
    path <- getOption("patRoon.test.DAAnalyses")
    if (is.null(path))
        return(NULL)
    ret <- generateAnalysisInfo(c(file.path(path, "neg"), file.path(path, "pos")),
                                groups = c(rep("solvent-neg", 3), rep("standard-neg", 3),
                                           rep("solvent-pos", 3), rep("standard-pos", 3)),
                                blanks = "solvent", formats = "bruker")
    if (!is.null(pat))
        ret <- ret[grepl(pat, ret$analysis), ]
    return(ret)
}

getMSLibMSPPath <- function() file.path(getTestDataPathGeneric(), "MoNA-export-CASMI_2012.msp")
getMSLibJSONPath <- function() file.path(getTestDataPathGeneric(), "MoNA-export-CASMI_2016.json")

if (testWithSets())
{
    isAnaInfoNeg <- function(anaInfo) grepl("\\-neg", anaInfo$analysis)
    
    getWorkPath <- function(file = "", ...) if (nzchar(file)) file.path("test_temp_sets", file, ...) else "test_temp_sets"
    getTestDataPath <- function() paste0(getTestDataPathGeneric(), "_sets")
    
    getTestAnaInfo <- function()
    {
        return(rbind(patRoonData::exampleAnalysisInfo("positive"), patRoonData::exampleAnalysisInfo("negative")))
    }
    getTestAnaInfoPos <- function(anaInfo = getTestAnaInfo()) anaInfo[!grepl("\\-neg", anaInfo$analysis), ]
    getDAAnaInfo <- function(...) getDAAnaInfoAll(...)
    getTestFeatures <- function(anaInfo = getTestAnaInfo(), noiseThrInt = 3E4, ...)
    {
        an <- isAnaInfoNeg(anaInfo)
        if (any(an))
            ret <- makeSet(findFeatures(anaInfo[!an, ], "openms", noiseThrInt = noiseThrInt, ...),
                           findFeatures(anaInfo[an, ], "openms", noiseThrInt = noiseThrInt, ...),
                           adducts = c("[M+H]+", "[M-H]-"))
        else
            ret <- makeSet(findFeatures(anaInfo, "openms", noiseThrInt = noiseThrInt, ...), adducts = "[M+H]+")
        return(ret)
        
    }
    getTestFGroupsDA <- function(anaInfo)
    {
        an <- isAnaInfoNeg(anaInfo)
        if (any(isAnaInfoNeg(anaInfo)))
            feats <- makeSet(findFeatures(anaInfo[!an, ], "bruker"),
                             findFeatures(anaInfo[an, ], "bruker"),
                             adducts = c("[M+H]+", "[M-H]-"))
        else
            feats <- makeSet(findFeatures(anaInfo[!an, ], "bruker"), adducts = c("[M+H]+", "[M-H]-"))
        return(groupFeatures(feats, "openms"))
    }
    makeOneEmptySetFGroups <- function(fGroups) delete(fGroups, which(analysisInfo(fGroups)$set == "negative"), function(...) TRUE)

    doExportXCMS <- function(x, ...)
    {
        # HACK: if x is a SIRIUS fGroups object, it will not be a set object
        if (is(x, "featureGroupsSIRIUS"))
            getXCMSSet(x, ...)
        else
            getXCMSSet(x, loadRawData = FALSE, set = "positive")
    }
    doExportXCMS3 <- function(x, ...)
    {
        # HACK: if x is a SIRIUS fGroups object, it will not be a set object
        if (is(x, "featureGroupsSIRIUS"))
            getXCMSnExp(x, ...)
        else
            getXCMSnExp(x, loadRawData = FALSE, set = "positive")
    }
    getExpAnaInfo <- function() getTestAnaInfoPos()
    getExpFeats <- function(x) x[, sets = "positive"]
    getExpFG <- function(x) x[, sets = "positive"]
    doExport <- function(x, ...) export(x, ..., set = "positive")
    
    getISTDAssignments <- function(fg) internalStandardAssignments(fg, "positive")

    getTestAnaInfoAnn <- function() getTestAnaInfo()[grepl("standard\\-.+\\-[2-3]", getTestAnaInfo()$analysis), ]
    getTestAnaInfoComponents <- function() getTestAnaInfo()[grepl("(solvent|standard)\\-.+\\-1", getTestAnaInfo()$analysis), ]

    getSIRFormFPsProjPath <- function() file.path(getTestDataPath(), paste0("SIRProjFormFPs", c("-pos", "-neg")))
    getSIRCompProjPath <- function() file.path(getTestDataPath(), paste0("SIRProjComp", c("-pos", "-neg")))
    
    doScreen <- function(fg, susp, ...)
    {
        cols <- !grepl("^mz$", names(susp))
        susp <- if (is.data.table(susp)) susp[, cols, with = FALSE] else susp[, cols, drop = FALSE]
        screenSuspects(fg, susp, ...)
    }
    
    doNormInts <- function(fg, ...) normInts(fg, ..., standards = list(patRoonData::ISTDListPos,
                                                                       patRoonData::ISTDListNeg))
    
    # zero threshold makes comparisons in testing much easier
    doGenForms <- function(..., setThresholdAnn = 0) generateFormulas(..., setThresholdAnn = setThresholdAnn)
    doFormCons <- function(..., setThresholdAnn = 0) consensus(..., setThresholdAnn = setThresholdAnn)
    doGenComps <- function(..., setThresholdAnn = 0) generateCompounds(..., setThresholdAnn = setThresholdAnn)
    doCompCons <- function(..., setThresholdAnn = 0) consensus(..., setThresholdAnn = setThresholdAnn)
    doGenComponents <- function(...) generateComponents(...)
    
    doGenLogicTPs <- function(...) generateTPs("logic", ...)
} else
{
    getWorkPath <- function(file = "", ...) if (nzchar(file)) file.path("test_temp", file, ...) else "test_temp"
    getTestDataPath <- getTestDataPathGeneric
    getTestAnaInfo <- function() patRoonData::exampleAnalysisInfo()
    getDAAnaInfo <- function(...)
    {
        ret <- getDAAnaInfoAll(pat)
        return(ret[grepl("-pos", ret$group, fixed = TRUE), ])
    }
    getTestFeatures <- function(anaInfo = getTestAnaInfo(), noiseThrInt = 3E4, ...)
    {
        ret <- findFeatures(anaInfo, "openms", noiseThrInt = noiseThrInt, ...)
        return(ret)
    }
    
    doExportXCMS <- function(x, ...) getXCMSSet(x, ...)
    doExportXCMS3 <- function(x, ...) getXCMSnExp(x, ...)
    getExpAnaInfo <- function() getTestAnaInfo()
    getExpFeats <- function(x) x
    getExpFG <- function(x) x
    doExport <- function(x, ...) export(x, ...)
    getISTDAssignments <- function(fg) internalStandardAssignments(fg)

    getTestAnaInfoComponents <- function() getTestAnaInfo()[3:4, ]
    getTestAnaInfoAnn <- function() getTestAnaInfo()[4:5, ]
    
    getSIRFormFPsProjPath <- function() file.path(getTestDataPath(), "SIRProjFormFPs")
    getSIRCompProjPath <- function() file.path(getTestDataPath(), "SIRProjComp")
    
    doScreen <- function(fg, susp, ...)
    {
        if ((is.null(susp[["mz"]]) || any(is.na(susp$mz))) && (is.null(susp[["adduct"]]) || any(is.na(susp$adduct))))
            screenSuspects(fg, susp, ..., adduct = "[M+H]+")
        else
            screenSuspects(fg, susp, ...)
    }
    doNormInts <- function(fg, ...) normInts(fg, ..., standards = patRoonData::ISTDListPos)
    doGenForms <- function(...) generateFormulas(..., adduct = "[M+H]+")
    doFormCons <- function(...) consensus(...)
    doGenComps <- function(...) generateCompounds(..., adduct = "[M+H]+")
    doCompCons <- function(...) consensus(...)
    doGenComponents <- function(fGroups, algorithm, ...)
    {
        args <- c(list(fGroups, algorithm), list(...))
        if (!algorithm %in% c("intclust", "specclust"))
            args <- c(args, list(ionization = "positive"))
        do.call(generateComponents, args)
    }
    
    doGenLogicTPs <- function(...) generateTPs("logic", ..., adduct = "[M+H]+")
}


# to make testing a bit easier: precursors don't have to follow filter rules
removePrecursors <- function(plists)
{
    doRmPrecs <- function(obj)
    {
        obj@peakLists <- pruneList(lapply(obj@peakLists,
                                             function(pa) pruneList(lapply(pa, function(pg)
                                                 pruneList(lapply(pg, function(pl) pl[precursor == FALSE]), checkZeroRows = TRUE)),
                                                 checkEmptyElements = TRUE)), checkEmptyElements = TRUE)
        obj@averagedPeakLists <- pruneList(lapply(obj@averagedPeakLists,
                                                     function(pg) pruneList(lapply(pg, function(pl) pl[precursor == FALSE]),
                                                                            checkZeroRows = TRUE)),
                                              checkEmptyElements = TRUE)
        return(obj)
    }
    
    plists <- doRmPrecs(plists)
    
    if (testWithSets())
        plists@setObjects <- lapply(plists@setObjects, doRmPrecs)
    
    return(plists)
}

removeMSPlists <- function(plists, type)
{
    clearpl <- function(pl) { pl[[type]] <- NULL; return(pl) }
    
    doRemove <- function(obj)
    {
        obj@peakLists <- lapply(obj@peakLists,
                                function(pa) pruneList(lapply(pa, clearpl), checkEmptyElements = TRUE))
        obj@averagedPeakLists <- pruneList(lapply(obj@averagedPeakLists, clearpl), checkEmptyElements = TRUE)
        return(obj)
    }
    
    plists <- doRemove(plists)
    
    if (testWithSets())
        plists@setObjects <- lapply(plists@setObjects, doRemove)
    
    return(plists)
}

doDATests <- function() !is.null(getOption("patRoon.test.DAAnalyses"))

makeMZXMLs <- function(anaInfo)
{
    # BUG/HACK: mzR cannot read mzXML files produced by OpenMS. The latter sets the .mzML source file as parentFile,
    # which is currently not supported by pwiz (pwiz keeps the original instrument raw data format). However, pwiz does
    # recognize mzXML as source format so simply convert the files twice...
    convertMSFiles(anaInfo = anaInfo, from = "mzML", outPath = getWorkPath(), to = "mzXML",
                   algorithm = "openms", overWrite = TRUE)
    anaInfo$path <- getWorkPath()
    convertMSFiles(anaInfo = anaInfo, from = "mzXML", outPath = getWorkPath(), to = "mzXML",
                   algorithm = "openms", overWrite = TRUE)
    return(anaInfo)
}

updateSIRIUSAnnProj <- function(SIRPath, clearEmptyFI)
{
    for (pp in normalizePath(SIRPath))
    {
        if (!dir.exists(pp))
            next
        
        if (clearEmptyFI)
        {
            # remove directories without FingerID results, as these will trigger new online searches even with dryRun=TRUE
            allDirs <- dirname(Sys.glob(file.path(pp, "*", "spectrum.ms")))
            dirsWithFI <- dirname(Sys.glob(file.path(pp, "*", "fingerid")))
            unlink(setdiff(allDirs, dirsWithFI), recursive = TRUE)
        }
        
        # zip to save space, .sirius extension seems to be needed. Change dir to fix the root.
        zipf <- paste0(pp, ".sirius")
        unlink(zipf)
        withr::with_dir(pp, zip(zipf, Sys.glob("*")))
    }
}

doGenFormsSIRFPs <- function(fGroups, plists) doGenForms(fGroups, plists, "sirius", dryRun = TRUE, calculateFeatures = FALSE,
                                                         getFingerprints = TRUE,
                                                         projectPath = paste0(getSIRFormFPsProjPath(), ".sirius"))
doGenCompsSIR <- function(fGroups, plists) doGenComps(fGroups, plists, "sirius", dryRun = TRUE,
                                                      projectPath = paste0(getSIRCompProjPath(), ".sirius"))

updateSIRIUSFormFPsProj <- function(...)
{
    unlink(getSIRFormFPsProjPath(), recursive = TRUE)
    withOpt(cache.mode = "none", doGenForms(..., algorithm = "sirius", projectPath = getSIRFormFPsProjPath(),
                                            calculateFeatures = FALSE, getFingerprints = TRUE))
    updateSIRIUSAnnProj(getSIRFormFPsProjPath(), clearEmptyFI = FALSE)
}

updateSIRIUSCompProj <- function(...)
{
    unlink(getSIRCompProjPath(), recursive = TRUE)
    withOpt(cache.mode = "none", doGenComps(..., algorithm = "sirius", projectPath = getSIRCompProjPath()))
    updateSIRIUSAnnProj(getSIRCompProjPath(), clearEmptyFI = TRUE)
}

expect_file <- function(object, file, removeIfExists = TRUE)
{
    if (removeIfExists && file.exists(file))
        file.remove(file)

    act <- quasi_label(rlang::enquo(object))
    expect(file.exists(file), sprintf("failed to generate %s", file))
    invisible(act$val)
}

expect_range <- function(object, r)
{
    act <- quasi_label(rlang::enquo(object))
    act$r <- range(act$val)
    expect(numGTE(act$r[1], r[1]) && numLTE(act$r[2], r[2]),
           sprintf("range of %s is %f - %f which is outside %f - %f",
                   act$lab, act$r[1], act$r[2], r[1], r[2]))
    invisible(act$val)
}

expect_gt_or_zero <- function(object, expected)
{
    act <- quasi_label(rlang::enquo(object))
    expect(all(object == 0 | object > expected), "object not greater than or zero")
    invisible(act$val)
}

expect_known_show <- function(object, file)
{
    act <- quasi_label(rlang::enquo(object))

    text <- capture_output_lines(show(act$val))

    # remove object size as it may vary even if object remain the same
    text <- text[!grepl("Object size", text)]

    text <- paste0(text, collapse = "")

    # based on simplified expect_known_output
    if (!file.exists(file))
    {
        warning("Creating reference output", call. = FALSE)
        cat(text, file = file)
        succeed()
    }
    else
    {
        ref <- patRoon:::readAllFile(file)
        cat(text, file = file)

        cmp <- compare(text, enc2native(ref))
        expect(cmp$equal, sprintf("show reference of %s has changed\n%s", act$lab, cmp$message))
    }

    invisible(act$val)
}

expect_plot <- function(object)
{
    tf <- tempfile()
    withr::with_png(tf, act <- quasi_label(rlang::enquo(object)))
    expect(file.exists(tf), "failed to generate plot")
    invisible(act$val)
}

expect_ggplot <- function(object)
{
    act <- quasi_label(rlang::enquo(object))
    expect(!is.null(object), "NULL plot")

    if (!is.null(object))
    {
        tf <- tempfile()
        withr::with_png(tf, print(object))
        expect(file.exists(tf), "failed to generate plot")
    }

    invisible(act$val)
}

makeReportHTML <- function(fGroups, path = getWorkPath(), overrideSettings = list(), ...)
{
    overrideSettings <- modifyList(
        list(
            general = list(noDate = TRUE), # for reproducibility
            # limits to speed up a bit
            formulas = list(topMost = 5),
            compounds = list(topMost = 5)
        ),
        overrideSettings, keep.null = TRUE
    )
    
    report(fGroups, path = path, openReport = FALSE, parallel = FALSE, overrideSettings = overrideSettings, ...)
}

expect_reportHTML <- function(object)
{
    # generate report twice: without and with cached results
    
    rpFile <- getWorkPath("report.html")
    unlink(rpFile) # in case it already exists
    unlink(getWorkPath("report_files"), recursive = TRUE)
    
    clearCache("reportHTML") # reset cached plots
    act <- quasi_label(rlang::enquo(object))
    
    expect(file.exists(rpFile), "failed to generate report")

    uniqueLines <- function(path)
    {
        # HACK: bslib sets random IDs on creation --> remove IDs to allow reproducible report generation
        lines <- readLines(path)
        lines <- gsub("data-tabsetid=\"[[:digit:]]+\"", "", lines)
        return(gsub("\"[#]?tab\\-[[:digit:]]+\\-[[:digit:]]+\"", "", lines))
    }
    
    if (file.exists(rpFile))
    {
        rpLines <- uniqueLines(rpFile)
        act <- quasi_label(rlang::enquo(object))
        expect(all.equal(rpLines, uniqueLines(rpFile)), "cached report differs")
    }

    invisible(act$val)
}

expect_HTML <- function(object)
{
    # UNDONE: use vdiffr when it supports it (https://github.com/r-lib/vdiffr/issues/60)
    
    out <- tempfile(fileext = ".html")
    htmlwidgets::saveWidget({
        act <- quasi_label(rlang::enquo(object))
        act$val
    }, out)
    expect(file.exists(out), "failed to plot HTML widget")
    invisible(act$val)
}

expect_equal_scr <- function(object, expected, ...)
{
    act <- quasi_label(rlang::enquo(object))
    rmCols <- c("InChI", "InChIKey", "SMILES", "neutralMass", "formula", "adduct")
    expect(isTRUE(all.equal(object[, setdiff(names(object), rmCols), with = FALSE], expected[, setdiff(names(expected), rmCols), with = FALSE], ...)),
           "screening results differ")
    invisible(act$val)
}

expect_doppel <- function(title, fig)
{
    # use different path for sets
    if (testWithSets())
    {
        # from vdiffr::expect_doppelganger(): get current context
        # context <- get(".context", envir = testthat::get_reporter())
        # return(vdiffr::expect_doppelganger(..., path = paste0(context, "_set")))
        # UNDONE: for now modify title, hopefully we can separate nicely again with the new snapshot system some day...
        title <- paste0("sets-", title)
    }
    return(vdiffr::expect_doppelganger(title, fig))
}

# HACK: workaround for non imported checkmate namespace
makeExpectation <- checkmate::makeExpectation
vname <- checkmate::vname

expect_csv_file <- checkmate::makeExpectationFunction(patRoon:::checkCSVFile)

# call dollar operator with string
callDollar <- function(x, name) eval(substitute(x$NAME_ARG, list(NAME_ARG = name)))

initXCMS <- function()
{
    # temporary workaround, fixed in https://github.com/sneumann/xcms/pull/462
    if (utils::packageVersion("xcms") < "3.10")
        library(xcms)
}

testFeatAnnADT <- function(obj)
{
    expect_setequal(as.data.table(obj)$group, groupNames(obj))
    expect_equal(nrow(as.data.table(obj)), length(obj))
    
    checkmate::expect_names(names(as.data.table(obj, fragments = TRUE)),
                            must.include = c("frag_ion_formula", "frag_mz"))
    expect_gt(nrow(as.data.table(obj, fragments = TRUE)), length(obj))
    
    checkmate::expect_names(names(as.data.table(obj, countElements = c("C", "H"))),
                            must.include = c("C", "H"))
    checkmate::expect_names(names(as.data.table(obj, countFragElements = c("C", "H"))),
                            must.include = c("frag_C", "frag_H"))

    OMTab <- as.data.table(obj, OM = TRUE)
    checkmate::qexpectr(OMTab[, c(unlist(strsplit("CHNOPS", "")), paste0(unlist(strsplit("HNOPS", "")), "C"),
                                  "DBE_AI", "AI")], "N+")
    checkmate::expect_character(OMTab[["classification"]], min.chars = 1, any.missing = FALSE, len = nrow(OMTab))
}

testSIRFPSubset <- function(obj)
{
    grp <- names(which.max(sapply(obj@fingerprints, ncol)))[1] # get fg with FP results for multiple candidates
    candidates <- intersect(names(obj@fingerprints[[grp]]), obj[[grp]]$neutral_formula)
    delOne <- delete(obj, i = grp, j = function(ann, ...) ann$neutral_formula == candidates[1])
    delAll <- delete(obj, i = grp)
    if (length(candidates) == 1)
        # UNDONE: unfortunately SIRIUS compounds test results have only one formula candidate for each group, so we
        # cannot test this. But for formulas we do, and the code used should be the same so we stick with the current
        # situation for now...
        expect_null(delAll@fingerprints[[grp]])
    else
        checkmate::expect_names(names(delOne@fingerprints[[grp]]), disjunct.from = candidates[1])
}

getOldNewRefs <- function(ref)
{
    oldref <- tempfile(fileext = ".Rds")
    ref <- file.path("tests", "testthat", testFile(ref))
    executeCommand("git", c("show", paste0("HEAD:", ref)), stdout = oldref)
    return(list(old = readRDS(oldref), new = readRDS(ref)))
}

compareRef <- function(ref)
{
    on <- getOldNewRefs(ref)
    waldo::compare(on$old, on$new, tolerance = 0.00001)
}
