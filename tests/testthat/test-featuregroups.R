context("feature groups")

initXCMS()

fList <- findFeatures(getTestAnaInfo(), "openms", logPath = NULL)
fgOpenMS <- groupFeatures(fList, "openms")
fgXCMS <- groupFeatures(fList, "xcms")
fgXCMS3 <- groupFeatures(fList, "xcms3")

# fList with dummy concs
anaInfoConc <- cbind(getTestAnaInfo(), list(conc = c(NA, NA, NA, 1, 2, 3)))
# modify replicate groups so we can test averaging
anaInfoConc$group[grepl("standard", anaInfoConc$group)] <- c("standard-1", "standard-2", "standard-2")
fListConc <- findFeatures(anaInfoConc, "openms", logPath = NULL)
fgOpenMSConc <- groupFeatures(fListConc, "openms")

fListEmpty <- findFeatures(getTestAnaInfo(), "openms", noiseThrInt = 1E9, logPath = NULL)
fgOpenMSEmpty <- groupFeatures(fListEmpty, "openms")
fgXCMSEmpty <- groupFeatures(fListEmpty, "xcms")
fgXCMS3Empty <- groupFeatures(fListEmpty, "xcms3")

test_that("verify feature grouping output", {
    expect_known_value(groups(fgOpenMS), testFile("fg-openms"))
    expect_known_value(groups(fgXCMS), testFile("fg-xcms"))
    expect_known_value(groups(fgXCMS3), testFile("fg-xcms3"))
    
    # extraOpts
    expect_equal(fgOpenMS, groupFeatures(fList, "openms",
                                         extraOptsRT = list("-algorithm:pairfinder:distance_RT:max_difference" = 30)))
    expect_equal(fgOpenMS, groupFeatures(fList, "openms",
                                         extraOptsGroup = list("-algorithm:distance_RT:max_difference" = 12)))
})

test_that("verify show output", {
    expect_known_show(fgOpenMS, testFile("fg-show-openms", text = TRUE))
    expect_known_show(fgXCMS, testFile("fg-show-xcms", text = TRUE))
    expect_known_show(fgXCMS3, testFile("fg-show-xcms3", text = TRUE))
})

test_that("empty objects work", {
    expect_length(fgOpenMSEmpty, 0)
    expect_length(fgXCMSEmpty, 0)
    expect_length(fgXCMS3Empty, 0)
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

    expect_length(fgOpenMS[[1]], length(analyses(fgOpenMS)))
    expect_length(fgOpenMS[[1, 1]], 1)
    expect_equivalent(fgOpenMS[[4, 50]], groups(fgOpenMS)[[4, 50]])
    expect_equivalent(fgOpenMS[[analyses(fgOpenMS)[4], names(fgOpenMS)[50]]], groups(fgOpenMS)[[4, 50]])
    expect_equivalent(fgOpenMS[[4]], groups(fgOpenMS)[[4]])
    expect_equivalent(callDollar(fgOpenMS, names(fgOpenMS)[4]), fgOpenMS[[4]])
})

expfile <- file.path(getWorkPath(), "export.csv") # NOTE: will be removed prior to each test automatically
test_that("exporting works", {
    expect_file(export(fgOpenMS, "brukerpa", expfile), expfile)
    expect_file(export(fgOpenMS, "brukertasq", expfile), expfile)
    expect_file(export(fgOpenMS, "mzmine", expfile), expfile)
    expect_error(export(fgOpenMSEmpty, "brukerpa", expfile))
})

XCMSImpXCMS <- getXCMSSet(fgXCMS)
XCMSImpXCMS3 <- getXCMSSet(fgXCMS3, exportedData = FALSE)
XCMSImpOpenMS <- getXCMSSet(fgOpenMS, exportedData = FALSE)
test_that("XCMS conversion", {
    expect_equal(nrow(xcms::groups(XCMSImpXCMS)), length(fgXCMS))
    expect_equal(nrow(xcms::groups(XCMSImpXCMS3)), length(fgXCMS3))
    expect_equal(nrow(xcms::groups(XCMSImpOpenMS)), length(fgOpenMS))
    
    expect_known_value(xcms::groups(XCMSImpXCMS), testFile("fg-xcms_import_xcms"))
    expect_known_value(xcms::groups(XCMSImpXCMS3), testFile("fg-xcms_import_xcms3"))
    expect_known_value(xcms::groups(XCMSImpOpenMS), testFile("fg-xcms_import_openms"))
    
    expect_equal(unname(groups(importFeatureGroupsXCMS(XCMSImpXCMS, getTestAnaInfo()))), unname(groups(fgXCMS)))
    expect_equal(unname(groups(importFeatureGroupsXCMS(XCMSImpXCMS3, getTestAnaInfo()))), unname(groups(fgXCMS3)))
    expect_equal(unname(groups(importFeatureGroupsXCMS(XCMSImpOpenMS, getTestAnaInfo()))), unname(groups(fgOpenMS)))
})

XCMS3ImpXCMS <- getXCMSnExp(fgXCMS)
XCMS3ImpXCMS3 <- getXCMSnExp(fgXCMS3)
XCMS3ImpOpenMS <- getXCMSnExp(fgOpenMS)
test_that("XCMS3 conversion", {
    expect_equal(nrow(xcms::featureDefinitions(XCMS3ImpXCMS)), length(fgXCMS))
    expect_equal(nrow(xcms::featureDefinitions(XCMS3ImpXCMS3)), length(fgXCMS3))
    expect_equal(nrow(xcms::featureDefinitions(XCMS3ImpOpenMS)), length(fgOpenMS))
    
    expect_known_value(xcms::featureDefinitions(XCMS3ImpXCMS), testFile("fg-xcms3_import_xcms"))
    expect_known_value(xcms::featureDefinitions(XCMS3ImpXCMS3), testFile("fg-xcms3_import_xcms3"))
    expect_known_value(xcms::featureDefinitions(XCMS3ImpOpenMS), testFile("fg-xcms3_import_openms"))
    
    expect_equal(unname(groups(importFeatureGroupsXCMS3(XCMS3ImpXCMS, getTestAnaInfo()))), unname(groups(fgXCMS)))
    expect_equal(unname(groups(importFeatureGroupsXCMS3(XCMS3ImpXCMS3, getTestAnaInfo()))), unname(groups(fgXCMS3)))
    expect_equal(unname(groups(importFeatureGroupsXCMS3(XCMS3ImpOpenMS, getTestAnaInfo()))), unname(groups(fgOpenMS)))
})

regr <- as.data.table(fgOpenMSConc, features = TRUE, regression = TRUE)
test_that("as.data.table works", {
    expect_equal(nrow(as.data.table(fgOpenMS)), length(fgOpenMS))

    # first 3 cols contain general info, then rep group ints
    expect_equal(ncol(as.data.table(fgOpenMS, average = TRUE)), 3 + length(unique(getTestAnaInfo()$group)))

    # UNDONE: intensities are sometimes higher than areas?
    # expect_gt_or_zero(as.data.table(fgOpenMS, areas = TRUE), as.data.table(fgOpenMS, areas = FALSE))
    # check if area from first group of first analysis corresponds to its feature data
    expect_equal(as.data.table(fgOpenMS, areas = TRUE)[[analyses(fgOpenMS)[1]]][1],
                 featureTable(fgOpenMS)[[analyses(fgOpenMS)[1]]][["area"]][groupFeatIndex(fgOpenMS)[[c(1, 1)]]])
    
    expect_range(nrow(as.data.table(fgOpenMS, features = TRUE)), length(fgOpenMS) * c(1, length(analyses(fgOpenMS))))

    expect_warning(as.data.table(fgOpenMS, regression = TRUE)) # no conc specified
    checkmate::expect_names(names(regr), must.include = "RSQ")
    checkmate::expect_names(names(as.data.table(fgOpenMSConc, features = FALSE, regression = TRUE)),
                            must.include = "RSQ")
    checkmate::expect_names(names(as.data.table(fgOpenMSConc, features = FALSE, average = TRUE,
                                                regression = TRUE)), must.include = "RSQ")
    expect_true(all(is.na(regr$conc) | is.na(regr$conc_reg) | regr$RSQ < 0.9 |
                        abs(regr$conc - regr$conc_reg) < 0.5)) # calculated concentrations should be somewhat close


    expect_equal(nrow(as.data.table(fgOpenMSEmpty, average = TRUE)), 0)
    expect_equal(nrow(as.data.table(fgOpenMSEmpty, features = TRUE)), 0)
    expect_equal(nrow(as.data.table(fgOpenMSEmpty, average = TRUE, features = TRUE)), 0)
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

minInt <- function(fg, rel)
{
    # collapse to vector with use.names = FALSE: https://stackoverflow.com/a/12796124/9264518
    g <- unlist(groups(fg), use.names = FALSE)
    if (rel)
        return(min(g[g != 0]) / max(g))
    return(min(g[g != 0]))
}

test_that("basic filtering", {
    expect_gte(minInt(filter(fgOpenMS, absMinIntensity = 1500), FALSE), 1500)
    expect_gte(minInt(filter(fgOpenMS, relMinIntensity = 0.2), TRUE), 0.2)
    expect_gte(minInt(filter(fgOpenMS, preAbsMinIntensity = 1500), FALSE), 1500)
    expect_gte(minInt(filter(fgOpenMS, preRelMinIntensity = 0.2), TRUE), 0.2)

    expect_range(groupInfo(filter(fgOpenMS, retentionRange = c(120, 200)))$rts, c(120, 200))
    expect_equivalent(filter(fgOpenMS, retentionRange = c(0, Inf)), fgOpenMS)
    expect_range(groupInfo(filter(fgOpenMS, mzRange = c(200, 300)))$mzs, c(200, 300))
    expect_equivalent(filter(fgOpenMS, mzRange = c(0, Inf)), fgOpenMS)
    expect_range(groupInfo(filter(fgOpenMS, mzDefectRange = c(0.1, 0.2)))$mzs %% 1, c(0.1, 0.2))
    expect_equivalent(filter(fgOpenMS, mzDefectRange = c(0, 1)), fgOpenMS)
    expect_lt(length(filter(fgOpenMS, chromWidthRange = c(0, 30))), length(fgOpenMS))
    expect_equivalent(filter(fgOpenMS, chromWidthRange = c(0, Inf)), fgOpenMS)

    expect_identical(replicateGroups(filter(fgOpenMS, rGroups = "standard")), "standard")
    expect_identical(replicateGroups(fgOpenMS[, rGroups = "standard"]), "standard")
    expect_identical(replicateGroups(filter(fgOpenMS, removeBlanks = TRUE)), "standard")
    expect_identical(replicateGroups(removeEmptyAnalyses(filter(fgOpenMS, relMinFeatures = 0.7))), "standard")
    expect_identical(replicateGroups(removeEmptyAnalyses(filter(fgOpenMS, absMinFeatures = 400))), "standard")
    expect_identical(replicateGroups(removeEmptyAnalyses(filter(fgOpenMS, blankThreshold = 1E6))), "standard")

    expect_known_output(filter(fgOpenMS, relMinAnalyses = 0.5), testFile("fgf-minana-rel", text = TRUE))
    expect_known_output(filter(fgOpenMS, absMinAnalyses = 3), testFile("fgf-minana-abs", text = TRUE))
    expect_known_output(filter(fgOpenMS, relMinReplicates = 1), testFile("fgf-minrep-rel", text = TRUE))
    expect_known_output(filter(fgOpenMS, absMinReplicates = 2), testFile("fgf-minrep-abs", text = TRUE))
    expect_known_output(filter(fgOpenMS, relMinFeatures = 0.75), testFile("fgf-minfeat-rel", text = TRUE))
    expect_known_output(filter(fgOpenMS, absMinFeatures = 450), testFile("fgf-minfeat-abs", text = TRUE))
    expect_known_output(filter(fgOpenMS, relMinReplicateAbundance = 1), testFile("fgf-minrepabu-rel", text = TRUE))
    expect_known_output(filter(fgOpenMS, absMinReplicateAbundance = 3), testFile("fgf-minrepabu-abs", text = TRUE))
    expect_known_output(filter(fgOpenMS, maxReplicateIntRSD = 0.5), testFile("fgf-reprsd", text = TRUE))
    expect_known_output(filter(fgOpenMS, blankThreshold = 5), testFile("fgf-bl", text = TRUE))
    expect_known_output(filter(fgOpenMS, absMinIntensity = 1500, blankThreshold = 5,
                               retentionRange = c(120, Inf), relMinReplicateAbundance = 1),
                        testFile("fgf-combi", text = TRUE))
    expect_known_output(filter(fgOpenMS, absMinIntensity = 1500, blankThreshold = 5,
                               retentionRange = c(120, Inf), relMinReplicateAbundance = 1, negate = TRUE),
                        testFile("fgf-combi-neg", text = TRUE))
    expect_length(filter(fgOpenMSEmpty, absMinIntensity = 1500, blankThreshold = 5,
                         retentionRange = c(120, Inf), relMinReplicateAbundance = 1), 0)
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
# fGCompXCMS3 <- comparison(openms = fgOpenMS, xcms = fgXCMS, groupAlgo = "xcms3")
fgCompOneEmpty <- comparison(openms = fgOpenMS, xcms = fgXCMSEmpty, groupAlgo = "openms")
fGConsOneEmpty <- consensus(fgCompOneEmpty)
fgCompBothEmpty <- comparison(openms = fgOpenMSEmpty, xcms = fgXCMSEmpty, groupAlgo = "openms")
fGConsBothEmpty <- consensus(fgCompBothEmpty)

test_that("verify feature group comparison", {
    expect_known_value(groups(fGCompOpenMS@comparedFGroups), testFile("fg-comp-openms"))
    expect_known_value(groups(fGCompXCMS@comparedFGroups), testFile("fg-comp-xcms"))

    expect_named(fGCompOpenMS, c("openms", "xcms"))
    expect_named(fGCompXCMS, c("openms", "xcms"))
    # expect_named(fGCompXCMS3, c("openms", "xcms"))
    expect_named(fGCompOpenMS[1], "openms")
    expect_named(fGCompOpenMS[2], "xcms")

    expect_equivalent(fGCompOpenMS[[1]], fgOpenMS)
    expect_equivalent(fGCompOpenMS[[names(fGCompOpenMS)[2]]], fgXCMS)
    expect_equivalent(fGCompOpenMS[[names(fGCompOpenMS)[2]]], fgXCMS)
    expect_equivalent(callDollar(fGCompOpenMS, names(fGCompOpenMS)[1]), fgOpenMS)

    expect_known_value(groups(fGCons), testFile("fg-comp-cons"))

    expect_lt(length(consensus(fGCompOpenMS, relMinAbundance = 1)), length(fGCons))
    expect_length(fgCompOneEmpty, 2)
    expect_length(fGConsOneEmpty, length(fgOpenMS))
    expect_length(fgCompBothEmpty, 2)
    expect_length(fGConsBothEmpty, 0)

    expect_equal(length(consensus(fGCompOpenMS, uniqueFrom = 1)) +
                 length(consensus(fGCompOpenMS, uniqueFrom = 2)) +
                 length(consensus(fGCompOpenMS, relMinAbundance = 1)), length(fGCons))
    expect_equal(length(consensus(fGCompOpenMS, uniqueFrom = 1:2, uniqueOuter = TRUE)) +
                 length(consensus(fGCompOpenMS, relMinAbundance = 1)), length(fGCons))
    expect_length(consensus(fGCompOpenMS, uniqueFrom = 1:2), length(fGCons))
    expect_lt(length(consensus(fGCompOpenMS, uniqueFrom = 1:2, uniqueOuter = TRUE)), length(fGCons))
    expect_length(consensus(fgCompBothEmpty, uniqueFrom = 1), 0)
    expect_length(consensus(fgCompBothEmpty, uniqueFrom = 1, uniqueOuter = TRUE), 0)
})

subFGroups <- fgOpenMS[, 1:25]
test_that("reporting works", {
    expect_file(reportCSV(subFGroups, getWorkPath(), reportFeatures = TRUE),
                file.path(getWorkPath(), sprintf("%s.csv", class(subFGroups))))
    for (ana in getTestAnaInfo()$analysis)
        checkmate::expect_file_exists(file.path(getWorkPath(), "features",
                                                sprintf("%s-%s.csv", class(getFeatures(subFGroups)), ana)))

    expect_file(reportPDF(subFGroups, getWorkPath()), getWorkPath(sprintf("%s.pdf", class(subFGroups))))

    expect_reportHTML(makeReportHTML(subFGroups))

    # skip if pngquant is not specified and not in PATH
    # assign condition to variable as expression seems to be to complicated for skip...
    havePngQuant <- (!is.null(getOption("patRoon.path.pngquant")) && nzchar(getOption("patRoon.path.pngquant"))) ||
        nzchar(Sys.which(sprintf("pngquant%s", if (Sys.info()[["sysname"]] == "Windows") ".exe" else "")))
    skip_if_not(havePngQuant)
    expect_error(reportHTML(subFGroups, getWorkPath("pngquant"), optimizePng = TRUE, openReport = FALSE), NA)
    expect_lt(file.size(getWorkPath("pngquant", "report.html")), file.size(getWorkPath("report.html")))
})

test_that("reporting with empty object works", {
    expect_error(reportCSV(fgOpenMSEmpty, getWorkPath(), reportFeatures = TRUE), NA)
    expect_error(reportPDF(fgOpenMSEmpty, getWorkPath()), NA)
    expect_error(makeReportHTML(fgOpenMSEmpty), NA)
})

test_that("plotting works", {
    expect_doppel("retmz", function() plot(fgOpenMS, colourBy = "fGroups", showLegend = FALSE))
    expect_doppel("retmz-singlec", function() plot(fgOpenMS, colourBy = "none", col = "blue"))
    expect_doppel("retmz-rgroups", function() plot(fgOpenMS, colourBy = "rGroups"))
    expect_doppel("retmz-comp", function() plot(fGCompOpenMS, colourBy = "fGroups", showLegend = FALSE))

    expect_doppel("intensity-def", function() plotInt(fgOpenMS))
    expect_doppel("intensity-avg", function() plotInt(fgOpenMS, TRUE))

    expect_doppel("chord-def", function() plotChord(fgOpenMS))
    expect_doppel("chord-selflinks", function() plotChord(fgOpenMS, addSelfLinks = TRUE))
    expect_doppel("chord-nortmz", function() plotChord(fgOpenMS, addRetMzPlots = FALSE))
    expect_doppel("chord-outer", function() plotChord(fgOpenMS,
                                                      outerGroups = c("standard-1" = "grp1",
                                                                      "standard-2" = "grp2",
                                                                      "standard-3" = "grp2",
                                                                      "solvent-1" = "grp3",
                                                                      "solvent-2" = "grp4",
                                                                      "solvent-3" = "grp5")))
    expect_doppel("chord-comp", function() plotChord(fGCompOpenMS))
    expect_error(plotChord(unique(fgOpenMS, which = replicateGroups(fgOpenMS), outer = TRUE),
                           average = TRUE)) # stops with nothing to plot: no overlap
    expect_plot(plotChord(unique(fgOpenMS, which = replicateGroups(fgOpenMS), outer = TRUE),
                          average = TRUE, addSelfLinks = TRUE)) # unless there are self links

    expect_doppel("eic-def", function() plotEIC(subFGroups))
    expect_doppel("eic-rtmin", function() plotEIC(subFGroups, retMin = TRUE))
    expect_doppel("eic-tm1", function() plotEIC(subFGroups, topMost = 1))
    expect_doppel("eic-area", function() plotEIC(subFGroups, showPeakArea = TRUE))
    expect_doppel("eic-cbr", function() plotEIC(subFGroups, colourBy = "rGroups"))
    expect_doppel("eic-cbf", function() plotEIC(subFGroups, colourBy = "fGroups"))
    expect_doppel("eic-ann", function() plotEIC(subFGroups, annotate = "mz"))

    expect_doppel("venn", function() plotVenn(fgOpenMS))
    expect_doppel("venn-comp", function() plotVenn(fGCompOpenMS))
    expect_equal(expect_plot(plotVenn(fgOpenMS, which = c("solvent", "standard")))$areas[2],
                 length(filter(fgOpenMS, rGroups = "standard")))
    expect_equal(expect_plot(plotVenn(fGCompOpenMS))$areas[2], length(fgXCMS))
    expect_equal(expect_plot(plotVenn(fGCompOpenMS))$intersectionCounts,
                 length(consensus(fGCompOpenMS, relMinAbundance = 1)))

    # vdiffr doesn't work with UpSet
    expect_ggplot(plotUpSet(fgOpenMS))
    expect_ggplot(plotUpSet(fGCompOpenMS))
})

test_that("plotting empty objects works", {
    expect_doppel("retmz-empty", function() plot(fgOpenMSEmpty))
    expect_doppel("retmz-empty", function() plot(fgOpenMSEmpty, colourBy = "rGroups"))
    expect_doppel("retmz", function() plot(fGConsOneEmpty, colourBy = "fGroups", showLegend = FALSE)) # should be same as fgOpenMS
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

    expect_error(plotUpSet(fgOpenMSEmpty))
    expect_ggplot(plotUpSet(fGConsOneEmpty))
    expect_error(plotUpSet(fgCompBothEmpty))
})
