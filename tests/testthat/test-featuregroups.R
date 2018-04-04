context("feature groups")

fList <- findFeatures(getTestAnaInfo(), "openms", logPath = NULL)

fgOpenMS <- groupFeatures(fList, "openms")
fgXCMS <- groupFeatures(fList, "xcms")

test_that("verify feature grouping output", {
    expect_known_value(groups(fgOpenMS), testFile("fg-openms"))
    expect_known_value(groups(fgXCMS), testFile("fg-xcms"))
})

test_that("verify show output", {
    expect_known_show(fgOpenMS, testFile("fg-show-openms", text = TRUE))
    expect_known_show(fgXCMS, testFile("fg-show-xcms", text = TRUE))
})

test_that("basic subsetting", {
    expect_equal(length(fgOpenMS[, 1:50]), 50)
    expect_equivalent(analysisInfo(fgOpenMS[1:3]), getTestAnaInfo()[1:3, ])
    expect_equivalent(analysisInfo(fgOpenMS[c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE)]),
                     getTestAnaInfo()[c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE), ])
    expect_equivalent(analysisInfo(fgOpenMS[getTestAnaInfo()$analysis[4:6]]), getTestAnaInfo()[4:6, ])
    expect_equal(length(fgOpenMS[FALSE]), 0)
})

expfile <- file.path(getWorkPath(), "export.csv") # NOTE: will be removed prior to each test automatically
test_that("exporting works", {
    expect_file(export(fgOpenMS, "brukerpa", expfile), expfile)
    expect_file(export(fgOpenMS, "brukertasq", expfile), expfile)
    expect_file(export(fgOpenMS, "mzmine", expfile), expfile)
})

test_that("groupTable works", {
    expect_equal(nrow(groupTable(fgOpenMS)), length(fgOpenMS))
    # first 3 cols contain general info, then rep group ints
    expect_equal(ncol(groupTable(fgOpenMS, average = TRUE)), 3 + length(unique(getTestAnaInfo()$group)))
})

test_that("unique works", {
    # note: only have two rep groups
    
    expect_equivalent(unique(fgOpenMS, which = "standard"),
                      unique(fgOpenMS, which = "standard", relativeTo = "solvent"))
    expect_equivalent(unique(fgOpenMS, which = "standard"), unique(fgOpenMS, which = "standard", outer = TRUE))
    expect_lt(length(unique(fgOpenMS, which = "standard")), length(fgOpenMS))
    expect_equal(length(unique(fgOpenMS, which = c("standard", "solvent"))), length(fgOpenMS))
    expect_equal(length(unique(fgOpenMS, which = c("standard", "solvent"), outer = TRUE)), 0)
})

test_that("unique works", {
    # note: only have two rep groups
    
    expect_lt(length(overlap(fgOpenMS, which = c("standard", "solvent"))), length(fgOpenMS))
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
    expect_identical(unique(analysisInfo(filter(fgOpenMS, rGroups = "standard"))$group), "standard")
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
})

test_that("replicate group subtraction", {
    # should be as these are the only two rep groups
    expect_setequal(names(replicateGroupSubtract(fgOpenMS, "solvent")),
                    names(unique(fgOpenMS, which = "standard")))
})

fGCompOpenMS <- comparison(openms = fgOpenMS, xcms = fgXCMS, groupAlgo = "openms")
fGCons <- consensus(fGCompOpenMS)
fGCompXCMS <- comparison(openms = fgOpenMS, xcms = fgXCMS, groupAlgo = "xcms")

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
    skip_if_not((!is.null(getOption("patRoon.path.pngquant")) && nzchar(getOption("patRoon.path.pngquant"))) ||
                nzchar(Sys.which(sprintf("pngquant%s", if (Sys.info()[["sysname"]] == "Windows") ".exe" else ""))))
    expect_error(reportMD(subFGroups, getWorkPath("pngquant"), optimizePng = TRUE), NA)
    expect_lt(file.size(getWorkPath("pngquant", "report.html")), file.size(getWorkPath("report.html")))
})

test_that("plotting works", {
    vdiffr::expect_doppelganger("retmz", function() plot(fgOpenMS))
    vdiffr::expect_doppelganger("retmz-comp", function() plot(fGCompOpenMS))
    
    vdiffr::expect_doppelganger("intensity-def", function() plotInt(fgOpenMS))
    vdiffr::expect_doppelganger("intensity-avg", function() plotInt(fgOpenMS, TRUE))
    
    vdiffr::expect_doppelganger("chord-def", function() plotChord(fgOpenMS))
    vdiffr::expect_doppelganger("chord-selflinks", function() plotChord(fgOpenMS, addSelfLinks = TRUE))
    vdiffr::expect_doppelganger("chord-nortmz", function() plotChord(fgOpenMS, addRetMzPlots = FALSE))
    vdiffr::expect_doppelganger("chord-comp", function() plotChord(fGCompOpenMS))
    
    vdiffr::expect_doppelganger("eic-def", function() plotEIC(subFGroups))
    vdiffr::expect_doppelganger("eic-rtmin", function() plotEIC(subFGroups, retMin = TRUE))
    vdiffr::expect_doppelganger("eic-tm1", function() plotEIC(subFGroups, topMost = 1))
    vdiffr::expect_doppelganger("eic-area", function() plotEIC(subFGroups, showPeakArea = TRUE))
    vdiffr::expect_doppelganger("eic-cbr", function() plotEIC(subFGroups, colourBy = "rGroups"))
    vdiffr::expect_doppelganger("eic-cbf", function() plotEIC(subFGroups, colourBy = "fGroups"))
    vdiffr::expect_doppelganger("eic-ann", function() plotEIC(subFGroups, annotate = "mz"))
    
    vdiffr::expect_doppelganger("venn", function() plotVenn(fgOpenMS))
    vdiffr::expect_doppelganger("venn-comp", function() plotVenn(fGCompOpenMS))
})
