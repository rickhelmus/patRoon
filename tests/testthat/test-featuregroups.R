context("feature groups")

fList <- findFeatures(getTestAnaInfo(), "openms", logPath = NULL)
fgOpenMS <- groupFeatures(fList, "openms")
fgXCMS <- groupFeatures(fList, "xcms")

fListEmpty <- findFeatures(getTestAnaInfo(), "openms", thr = 1E9, logPath = NULL)
fgOpenMSEmpty <- groupFeatures(fListEmpty, "openms")
fgXCMSEmpty <- groupFeatures(fListEmpty, "xcms")

test_that("verify feature grouping output", {
    expect_known_value(groups(fgOpenMS), testFile("fg-openms"))
    expect_known_value(groups(fgXCMS), testFile("fg-xcms"))
})

test_that("verify show output", {
    expect_known_show(fgOpenMS, testFile("fg-show-openms", text = TRUE))
    expect_known_show(fgXCMS, testFile("fg-show-xcms", text = TRUE))
})

test_that("empty objects work", {
    expect_length(fgOpenMSEmpty, 0)
    expect_length(fgXCMSEmpty, 0)
})

test_that("basic subsetting", {
    expect_length(fgOpenMS[, 1:50], 50)
    expect_length(fgOpenMS[, "nope"], 0)
    expect_length(fgOpenMS["nope"], 0)
    expect_equivalent(analysisInfo(fgOpenMS[1:3]), getTestAnaInfo()[1:3, ])
    expect_equivalent(analysisInfo(fgOpenMS[c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE)]),
                     getTestAnaInfo()[c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE), ])
    expect_equivalent(analysisInfo(fgOpenMS[getTestAnaInfo()$analysis[4:6]]), getTestAnaInfo()[4:6, ])
    expect_equal(length(fgOpenMS[FALSE]), 0)
    expect_length(fgOpenMSEmpty[, 1:50], 0)
})

expfile <- file.path(getWorkPath(), "export.csv") # NOTE: will be removed prior to each test automatically
test_that("exporting works", {
    expect_file(export(fgOpenMS, "brukerpa", expfile), expfile)
    expect_file(export(fgOpenMS, "brukertasq", expfile), expfile)
    expect_file(export(fgOpenMS, "mzmine", expfile), expfile)
    expect_error(export(fgOpenMSEmpty, "brukerpa", expfile))
    expect_file(export(fgOpenMSEmpty, "brukertasq", expfile), expfile)
    expect_file(export(fgOpenMSEmpty, "mzmine", expfile), expfile)
})

test_that("groupTable works", {
    expect_equal(nrow(groupTable(fgOpenMS)), length(fgOpenMS))
    # first 3 cols contain general info, then rep group ints
    expect_equal(ncol(groupTable(fgOpenMS, average = TRUE)), 3 + length(unique(getTestAnaInfo()$group)))
    expect_equal(nrow(groupTable(fgOpenMSEmpty)), 0)
})

test_that("unique works", {
    # note: only have two rep groups
    
    expect_equivalent(unique(fgOpenMS, which = "standard"),
                      unique(fgOpenMS, which = "standard", relativeTo = "solvent"))
    expect_equivalent(unique(fgOpenMS, which = "standard"), unique(fgOpenMS, which = "standard", outer = TRUE))
    expect_lt(length(unique(fgOpenMS, which = "standard")), length(fgOpenMS))
    expect_equal(length(unique(fgOpenMS, which = c("standard", "solvent"))), length(fgOpenMS))
    expect_lt(length(unique(fgOpenMS, which = c("standard", "solvent"), outer = TRUE)), length(fgOpenMS))
    expect_length(unique(fgOpenMSEmpty, which = c("standard", "solvent")), 0)
})

test_that("overlap works", {
    # note: only have two rep groups
    
    expect_lt(length(overlap(fgOpenMS, which = c("standard", "solvent"))), length(fgOpenMS))
    expect_length(overlap(fgOpenMSEmpty, which = c("standard", "solvent")), 0)
})

minInt <- function(fg)
{
    # collapse to vector with use.names = FALSE: https://stackoverflow.com/a/12796124/9264518
    g <- unlist(groups(fg), use.names = FALSE)
    min(g[g != 0])
}

test_that("basic filtering", {
    expect_gte(minInt(filter(fgOpenMS, intensityThreshold = 500)), 500)
    
    expect_range(groupInfo(filter(fgOpenMS, retentionRange = c(120, 200)))$rts, range(120, 200))
    expect_equivalent(filter(fgOpenMS, retentionRange = c(0, -1)), fgOpenMS)
    expect_range(groupInfo(filter(fgOpenMS, mzRange = c(200, 300)))$mzs, range(200, 300))
    expect_equivalent(filter(fgOpenMS, mzRange = c(0, -1)), fgOpenMS)
    expect_lt(length(filter(fgOpenMS, chromWidthRange = c(0, 30))), length(fgOpenMS))
    expect_equivalent(filter(fgOpenMS, chromWidthRange = c(0, -1)), fgOpenMS)
    
    expect_identical(unique(analysisInfo(filter(fgOpenMS, rGroups = "standard"))$group), "standard")
    expect_identical(unique(analysisInfo(filter(fgOpenMS, removeRefAnalyses = TRUE))$group), "standard")
    expect_known_output(filter(fgOpenMS, relAbundance = 0.5), testFile("fgf-relabu", text = TRUE))
    expect_known_output(filter(fgOpenMS, absAbundance = 3), testFile("fgf-absabu", text = TRUE))
    expect_known_output(filter(fgOpenMS, interRelRGroupAbundance = 1), testFile("fgf-inter_rel_rg", text = TRUE))
    expect_known_output(filter(fgOpenMS, interAbsRGroupAbundance = 2), testFile("fgf-inter_abs_rg", text = TRUE))
    expect_known_output(filter(fgOpenMS, intraRGroupAbundance = 1), testFile("fgf-intra_rg", text = TRUE))
    expect_known_output(filter(fgOpenMS, minBlankThreshold = 5), testFile("fgf-bl", text = TRUE))
    expect_known_output(filter(fgOpenMS, intensityThreshold = 500, minBlankThreshold = 5,
                               retentionRange = c(120, -1), intraRGroupAbundance = 1),
                        testFile("fgf-combi", text = TRUE))
    expect_known_output(filter(fgOpenMS, intensityThreshold = 500, minBlankThreshold = 5,
                               retentionRange = c(120, -1), intraRGroupAbundance = 1, negate = TRUE),
                        testFile("fgf-combi-neg", text = TRUE))
    expect_length(filter(fgOpenMSEmpty, intensityThreshold = 500, minBlankThreshold = 5,
                         retentionRange = c(120, -1), intraRGroupAbundance = 1), 0)
})

test_that("replicate group subtraction", {
    # should be as these are the only two rep groups
    expect_setequal(names(replicateGroupSubtract(fgOpenMS, "solvent")),
                    names(unique(fgOpenMS, which = "standard")))
    expect_length(replicateGroupSubtract(fgOpenMSEmpty, "solvent"), 0)
})

fGCompOpenMS <- comparison(openms = fgOpenMS, xcms = fgXCMS, groupAlgo = "openms")
fGCons <- consensus(fGCompOpenMS)
fGCompXCMS <- comparison(openms = fgOpenMS, xcms = fgXCMS, groupAlgo = "xcms")
fgCompOneEmpty <- comparison(openms = fgOpenMS, xcms = fgXCMSEmpty, groupAlgo = "openms")
fGConsOneEmpty <- consensus(fgCompOneEmpty)
fgCompBothEmpty <- comparison(openms = fgOpenMSEmpty, xcms = fgXCMSEmpty, groupAlgo = "openms")
fGConsBothEmpty <- consensus(fgCompBothEmpty)

test_that("verify feature group comparison", {
    expect_known_value(groups(fGCompOpenMS@comparedFGroups), testFile("fg-comp-openms"))
    expect_known_value(groups(fGCompXCMS@comparedFGroups), testFile("fg-comp-xcms"))
    expect_named(fGCompOpenMS, c("openms", "xcms"))
    expect_named(fGCompXCMS, c("openms", "xcms"))
    expect_named(fGCompOpenMS[1], "openms")
    expect_named(fGCompOpenMS[2], "xcms")
    expect_known_value(groups(fGCons), testFile("fg-comp-cons"))
    expect_lt(length(consensus(fGCompOpenMS, relAbundance = 1)), length(fGCons))
    expect_lt((length(unique(fGCompOpenMS, which = "openms")) +
                  length(unique(fGCompOpenMS, which = "xcms"))), length(fGCons))
    expect_length(fgCompOneEmpty, 2)
    expect_length(fGConsOneEmpty, length(fgOpenMS))
    expect_length(fgCompBothEmpty, 2)
    expect_length(fGConsBothEmpty, 0)
})

subFGroups <- fgOpenMS[, 1:25]
test_that("reporting works", {
    expect_file(reportCSV(subFGroups, getWorkPath(), reportFeatures = TRUE),
                file.path(getWorkPath(), sprintf("%s.csv", class(subFGroups))))
    for (ana in getTestAnaInfo()$analysis)
        checkmate::expect_file_exists(file.path(getWorkPath(), "features",
                                                sprintf("%s-%s.csv", class(getFeatures(subFGroups)), ana)))
    
    expect_file(reportPDF(subFGroups, getWorkPath()), getWorkPath(sprintf("%s.pdf", class(subFGroups))))
    
    expect_file(reportMD(subFGroups, getWorkPath()), getWorkPath("report.html"))
    
    # skip if pngquant is not specified and not in PATH
    # assign condition to variable as expression seems to be to complicated for skip...
    havePngQuant <- (!is.null(getOption("patRoon.path.pngquant")) && nzchar(getOption("patRoon.path.pngquant"))) ||
        nzchar(Sys.which(sprintf("pngquant%s", if (Sys.info()[["sysname"]] == "Windows") ".exe" else "")))
    skip_if_not(havePngQuant)
    expect_error(reportMD(subFGroups, getWorkPath("pngquant"), optimizePng = TRUE), NA)
    expect_lt(file.size(getWorkPath("pngquant", "report.html")), file.size(getWorkPath("report.html")))
})

test_that("reporting with empty object works", {
    expect_error(reportCSV(fgOpenMSEmpty, getWorkPath(), reportFeatures = TRUE), NA)
    expect_error(reportPDF(fgOpenMSEmpty, getWorkPath()), NA)
    expect_file(reportMD(fgOpenMSEmpty, getWorkPath()), getWorkPath("report.html"))
})

test_that("plotting works", {
    expect_doppel("retmz", function() plot(fgOpenMS))
    expect_doppel("retmz-comp", function() plot(fGCompOpenMS))
    
    expect_doppel("intensity-def", function() plotInt(fgOpenMS))
    expect_doppel("intensity-avg", function() plotInt(fgOpenMS, TRUE))
    
    expect_doppel("chord-def", function() plotChord(fgOpenMS))
    expect_doppel("chord-selflinks", function() plotChord(fgOpenMS, addSelfLinks = TRUE))
    expect_doppel("chord-nortmz", function() plotChord(fgOpenMS, addRetMzPlots = FALSE))
    expect_doppel("chord-comp", function() plotChord(fGCompOpenMS))
    
    expect_doppel("eic-def", function() plotEIC(subFGroups))
    expect_doppel("eic-rtmin", function() plotEIC(subFGroups, retMin = TRUE))
    expect_doppel("eic-tm1", function() plotEIC(subFGroups, topMost = 1))
    expect_doppel("eic-area", function() plotEIC(subFGroups, showPeakArea = TRUE))
    expect_doppel("eic-cbr", function() plotEIC(subFGroups, colourBy = "rGroups"))
    expect_doppel("eic-cbf", function() plotEIC(subFGroups, colourBy = "fGroups"))
    expect_doppel("eic-ann", function() plotEIC(subFGroups, annotate = "mz"))
    
    expect_doppel("venn", function() plotVenn(fgOpenMS))
    expect_doppel("venn-comp", function() plotVenn(fGCompOpenMS))
})

test_that("plotting empty objects works", {
    expect_doppel("retmz-empty", function() plot(fgOpenMSEmpty))
    expect_doppel("retmz", function() plot(fGConsOneEmpty)) # should be same as fgOpenMS
    expect_doppel("retmz-comp-empty", function() plot(fgCompBothEmpty))
    
    expect_doppel("intensity-def-empty", function() plotInt(fgOpenMSEmpty))
    expect_doppel("intensity-avg-empty", function() plotInt(fgOpenMSEmpty, TRUE))
    
    expect_error(plotChord(fgOpenMSEmpty))
    expect_doppel("chord-def", function() plotChord(fGConsOneEmpty)) # should be same as fgOpenMS
    expect_error(plotChord(fgCompBothEmpty))
    
    expect_doppel("eic-def-empty", function() plotEIC(fgOpenMSEmpty))

    expect_error(plotVenn(fgOpenMSEmpty))
    expect_doppel("venn", function() plotVenn(fGConsOneEmpty)) # should be same as fgOpenMS
    expect_error(plotVenn(fgCompBothEmpty))
})
