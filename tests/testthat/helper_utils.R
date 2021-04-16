testWithSets <- function() T # UNDONE: check environment variable or something

testFile <- function(f, ..., text = FALSE) file.path(getTestDataPath(), paste0(f, ..., if (!text) ".Rds" else ".txt", collapse = ""))
getTestFGroups <- function(anaInfo = getTestAnaInfo(), ...) groupFeatures(getTestFeatures(anaInfo, ...), "openms")
getEmptyFeatures <- function(anaInfo = getTestAnaInfo(), ...) getTestFeatures(anaInfo, noiseThrInt = 1E9, ...)
getEmptyTestFGroups <- function(anaInfo = getTestAnaInfo()) getTestFGroups(anaInfo)[, "none"]

getMFTestDBPath <- function() file.path(getTestDataPath(), "test-mf-db.csv")
getCompFGroups <- function()
{
    fGroups <- doScreen(getTestFGroups(getTestAnaInfoAnn()), patRoonData::suspectsPos, onlyHits = TRUE)
    # just focus on few targets also present in MF test DB
    return(fGroups[, suspects = fread(getMFTestDBPath())$Name])
}

callMF <- function(fGroups, plists, scoreTypes = "fragScore", db = getMFTestDBPath(), to = 300, ...)
{
    doGenComps(fGroups, plists, "metfrag", timeout = to, database = "csv", scoreTypes = scoreTypes,
               extraOpts = list(LocalDatabasePath = db), ...)
}

if (testWithSets())
{
    isAnaInfoNeg <- function(anaInfo) grepl("\\-neg", anaInfo$analysis)
    
    getWorkPath <- function(file = "", ...) if (nzchar(file)) file.path("test_temp_sets", file, ...) else "test_temp_sets"
    getTestDataPath <- function() "test_data_sets"
    getTestAnaInfo <- function()
    {
        return(rbind(patRoonData::exampleAnalysisInfo("positive"), patRoonData::exampleAnalysisInfo("negative")))
    }
    getTestAnaInfoPos <- function(anaInfo = getTestAnaInfo()) anaInfo[!grepl("\\-neg", anaInfo$analysis), ]
    getTestFeatures <- function(anaInfo = getTestAnaInfo(), ...)
    {
        an <- isAnaInfoNeg(anaInfo)
        if (any(an))
            ret <- makeSet(findFeatures(anaInfo[!an, ], "openms", ...),
                           findFeatures(anaInfo[an, ], "openms", ...),
                           adducts = c("[M+H]+", "[M-H]-"))
        else
            ret <- makeSet(findFeatures(anaInfo, "openms", ...), adducts = "[M+H]+")
        return(filter(ret, absMinIntensity = 5E4)) # reduce a bit as we don't need all of them for testing...
        
    }
    makeOneEmptySetFGroups <- function(fGroups) delete(fGroups, which(analysisInfo(fGroups)$set == "negative"), function(...) TRUE)

    doExportXCMS <- function(x, ...) getXCMSSet(x, exportedData = FALSE, set = "positive")
    doExportXCMS3 <- function(x, ...) getXCMSnExp(x, exportedData = FALSE, set = "positive")
    getExpAnaInfo <- function() getTestAnaInfoPos()
    getExpFeats <- function(x) x[, sets = "positive"]
    getExpFG <- function(x) x[, sets = "positive"]
    doExport <- function(x, ...) export(x, ..., set = "positive")

    getTestAnaInfoAnn <- function() getTestAnaInfo()[grepl("standard\\-.+\\-[2-3]", getTestAnaInfo()$analysis), ]
    getTestAnaInfoComponents <- function() getTestAnaInfo()[grepl("(solvent|standard)\\-.+\\-1", getTestAnaInfo()$analysis), ]
    
    doScreen <- function(fg, susp, ...) screenSuspects(fg, susp[, !grepl("mz", names(susp), fixed = TRUE), drop = FALSE], ...)
        
    # zero threshold makes comparisons in testing much easier
    doGenForms <- function(..., setThresholdAnn = 0) generateFormulas(..., setThresholdAnn = setThresholdAnn)
    doFormCons <- function(..., setThresholdAnn = 0) consensus(..., setThresholdAnn = setThresholdAnn)
    doGenComps <- function(..., setThresholdAnn = 0) generateCompounds(..., setThresholdAnn = setThresholdAnn)
    doCompCons <- function(..., setThresholdAnn = 0) consensus(..., setThresholdAnn = setThresholdAnn)
    doGenComponents <- function(...) generateComponents(...)
} else
{
    getWorkPath <- function(file = "", ...) if (nzchar(file)) file.path("test_temp", file, ...) else "test_temp"
    getTestDataPath <- function() "test_data"
    getTestAnaInfo <- function() patRoonData::exampleAnalysisInfo()
    getTestFeatures <- function(anaInfo = getTestAnaInfo(), ...)
    {
        ret <- findFeatures(anaInfo, "openms", ...)
        return(filter(ret, absMinIntensity = 5E4)) # reduce a bit as we don't need all of them for testing...
    }
    
    doExportXCMS <- function(x, ...) getXCMSSet(x, ...)
    doExportXCMS3 <- function(x, ...) getXCMSnExp(x, ...)
    getExpAnaInfo <- function() getTestAnaInfo()
    getExpFeats <- function(x) x
    getExpFG <- function(x) x
    doExport <- function(x, ...) export(x, ...)

    getTestAnaInfoComponents <- function() getTestAnaInfo()[3:4, ]
    getTestAnaInfoAnn <- function() getTestAnaInfo()[4:5, ]
    
    doScreen <- function(...) screenSuspects(...)
    doGenForms <- function(...) generateFormulas(..., adduct = "[M+H]+")
    doFormCons <- function(...) consensus(...)
    doGenComps <- function(...) generateCompounds(..., adduct = "[M+H]+")
    doCompCons <- function(...) consensus(...)
    doGenComponents <- function(fGroups, algorithm, ...)
    {
        args <- c(list(fGroups, algorithm), list(...))
        if (algorithm != "intclust")
            args <- c(args, list(ionization = "positive"))
        do.call(generateComponents, args)
    }
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

getDAAnaInfo <- function()
{
    path <- getOption("patRoon.test.DAAnalyses")
    if (is.null(path))
        return(NULL)
    return(generateAnalysisInfo(path, groups = c(rep("standard", 3), rep("solvent", 3)),
                                blanks = "solvent"))
}

doDATests <- function() !is.null(getOption("patRoon.test.DAAnalyses"))

makeMZXMLs <- function(anaInfo)
{
    convertMSFiles(anaInfo = anaInfo, from = "mzML", outPath = getWorkPath(), to = "mzXML",
                   algorithm = "openms")
    anaInfo$path <- getWorkPath()
    return(anaInfo)
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
    expect(act$r[1] >= r[1] && act$r[2] <= r[2],
           sprintf("range of %s is %.1f - %.1f which is outside %.1f - %.1f",
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

# noDate should be TRUE for consistent report generation (ie when verifying cached report)
makeReportHTML <- function(fGroups, ...) reportHTML(fGroups, getWorkPath(), openReport = FALSE, noDate = TRUE, ...)

expect_reportHTML <- function(object)
{
    # generate report twice: without and with cache
    
    rpFile <- getWorkPath("report.html")
    unlink(getWorkPath("report.html")) # in case it already exists
    
    withr::with_options(list(patRoon.cache.mode = "save"), act <- quasi_label(rlang::enquo(object)))
    
    expect(file.exists(rpFile), "failed to generate report")

    if (file.exists(rpFile))
    {
        rpHash <- makeFileHash(rpFile)
        act <- quasi_label(rlang::enquo(object))
        expect(rpHash == makeFileHash(rpFile), "cached report differs")
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

expect_doppel <- function(...)
{
    # use different path for sets
    if (testWithSets())
    {
        # from vdiffr::expect_doppelganger(): get current context
        context <- get(".context", envir = testthat::get_reporter())
        return(vdiffr::expect_doppelganger(..., path = paste0(context, "_set")))
    }
    return(vdiffr::expect_doppelganger(...))
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

