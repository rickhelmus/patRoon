context("feature groups")

fList <- findFeatures(testAnaInfo, "openms", logPath = NULL)

fgOpenMS <- groupFeatures(fList, "openms")
fgXCMS <- groupFeatures(fList, "xcms")

test_that("verify feature grouping output", {
    expect_known_value(groups(fgOpenMS), testFile("fg-openms"))
    expect_known_value(groups(fgXCMS), testFile("fg-xcms"))
})

test_that("verify show output", {
    expect_known_output(show(fgOpenMS), testFile("fg-show-openms", text = TRUE))
    expect_known_output(show(fgXCMS), testFile("fg-show-xcms", text = TRUE))
})

test_that("basic subsetting", {
    expect_equal(length(fgOpenMS[, 1:50]), 50)
    expect_equivalent(analysisInfo(fgOpenMS[1:3]), testAnaInfo[1:3, ])
    expect_equivalent(analysisInfo(fgOpenMS[c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE)]),
                     testAnaInfo[c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE), ])
    expect_equivalent(analysisInfo(fgOpenMS[testAnaInfo$analysis[4:6]]), testAnaInfo[4:6, ])
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
    expect_equal(ncol(groupTable(fgOpenMS, average = TRUE)), 3 + length(unique(testAnaInfo$group)))
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
