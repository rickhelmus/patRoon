context("feature groups")

initXCMS()

fList <- getTestFeatures()
fgOpenMS <- groupFeatures(fList, "openms")
fgXCMS <- groupFeatures(fList, "xcms")
fgXCMS3 <- groupFeatures(fList, "xcms3")

# fList with dummy concs
anaInfoConc <- cbind(getTestAnaInfo(), list(conc = c(NA, NA, NA, 1, 2, 3)))
# modify replicate groups so we can test averaging
anaInfoConc$group[grepl("standard", anaInfoConc$group)] <- c("standard-1", "standard-2", "standard-2")
fListConc <- findFeatures(anaInfoConc, "openms")
fgOpenMSConc <- groupFeatures(fListConc, "openms")

fListEmpty <- getEmptyFeatures()
fgOpenMSEmpty <- groupFeatures(fListEmpty, "openms")
fgXCMSEmpty <- groupFeatures(fListEmpty, "xcms")
fgXCMS3Empty <- groupFeatures(fListEmpty, "xcms3")

test_that("verify feature grouping output", {
    expect_known_value(groupTable(fgOpenMS), testFile("fg-openms"))
    expect_known_value(groupTable(fgXCMS), testFile("fg-xcms"))
    
    # extraOpts
    expect_equal(groupTable(fgOpenMS),
                 groupTable(groupFeatures(fList, "openms",
                                          extraOptsRT = list("-algorithm:pairfinder:distance_RT:max_difference" = 30))))
    expect_equal(groupTable(fgOpenMS),
                 groupTable(groupFeatures(fList, "openms",
                                          extraOptsGroup = list("-algorithm:distance_RT:max_difference" = 12))))
    expect_known_value(groupTable(fgXCMS3), testFile("fg-xcms3"))
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

if (!testWithSets())
{
    fgOpenMSAnn <- fgOpenMS
    adducts(fgOpenMSAnn) <- rep("[M+H]+", length(fgOpenMSAnn))
    fgOpenMSAnn2 <- fgOpenMSAnn
    adducts(fgOpenMSAnn2)[3] <- "[M+K]+"
}

test_that("adducts setting", {
    skip_if(testWithSets()) # sets tests done below
    
    expect_length(adducts(fgOpenMS), 0)
    expect_length(adducts(fgOpenMSAnn), length(fgOpenMSAnn))
    expect_length(adducts(fgOpenMSAnn2), length(fgOpenMSAnn2))
    expect_setequal(adducts(fgOpenMSAnn), "[M+H]+")
    expect_equal(adducts(fgOpenMSAnn2)[3], "[M+K]+", check.attributes = FALSE)
    
    # verify neutral masses
    expect_true(all(sapply(seq_len(nrow(annotations(fgOpenMSAnn2))), function(i)
    {
        ann <- annotations(fgOpenMSAnn2)[i]
        return(isTRUE(all.equal(ann$neutralMass + adductMZDelta(as.adduct(ann$adduct)),
                                groupInfo(fgOpenMSAnn2)[ann$group, "mzs"])))
    })))
})

# to compare original anaInfo: ignore extra columns that may have been added afterwards
getAnaInfo <- function(fg) analysisInfo(fg)[, c("path", "analysis", "group", "blank")]

test_that("basic subsetting", {
    expect_length(fgOpenMS[, 1:50], 50)
    expect_length(fgOpenMS[, "nope"], 0)
    expect_length(fgOpenMS["nope"], 0)
    expect_equivalent(getAnaInfo(fgOpenMS[1:3]), getTestAnaInfo()[1:3, ])
    expect_equivalent(getAnaInfo(fgOpenMS[c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE)]),
                      getTestAnaInfo()[c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE), ])
    expect_equivalent(getAnaInfo(fgOpenMS[getTestAnaInfo()$analysis[4:6]]), getTestAnaInfo()[4:6, ])
    expect_equal(length(fgOpenMS[FALSE]), 0)
    expect_length(fgOpenMSEmpty[, 1:50], 0)

    expect_length(fgOpenMS[[1]], length(analyses(fgOpenMS)))
    expect_length(fgOpenMS[[1, 1]], 1)
    expect_equivalent(fgOpenMS[[4, 50]], groupTable(fgOpenMS)[[4, 50]])
    expect_equivalent(fgOpenMS[[analyses(fgOpenMS)[4], names(fgOpenMS)[50]]], groupTable(fgOpenMS)[[4, 50]])
    expect_equivalent(fgOpenMS[[4]], groupTable(fgOpenMS)[[4]])
    expect_equivalent(callDollar(fgOpenMS, names(fgOpenMS)[4]), fgOpenMS[[4]])
})

expfile <- file.path(getWorkPath(), "export.csv") # NOTE: will be removed prior to each test automatically
test_that("exporting works", {
    expect_file(doExport(fgOpenMS, "brukerpa", expfile), expfile)
    expect_file(doExport(fgOpenMS, "brukertasq", expfile), expfile)
    expect_file(doExport(fgOpenMS, "mzmine", expfile), expfile)
    expect_error(doExport(fgOpenMSEmpty, "brukerpa", expfile))
})

XCMSImpXCMS <- doExportXCMS(fgXCMS)
XCMSImpXCMS3 <- doExportXCMS(fgXCMS3, exportedData = FALSE)
XCMSImpOpenMS <- doExportXCMS(fgOpenMS, exportedData = FALSE)
test_that("XCMS conversion", {
    expect_equal(nrow(xcms::groups(XCMSImpXCMS)), length(getExpFG(fgXCMS)))
    expect_equal(nrow(xcms::groups(XCMSImpXCMS3)), length(getExpFG(fgXCMS3)))
    expect_equal(nrow(xcms::groups(XCMSImpOpenMS)), length(getExpFG(fgOpenMS)))
    
    expect_known_value(xcms::groups(XCMSImpXCMS), testFile("fg-xcms_import_xcms"))
    expect_known_value(xcms::groups(XCMSImpXCMS3), testFile("fg-xcms_import_xcms3"))
    expect_known_value(xcms::groups(XCMSImpOpenMS), testFile("fg-xcms_import_openms"))
    
    expect_equal(unname(groupTable(importFeatureGroupsXCMS(XCMSImpXCMS, getExpAnaInfo()))),
                 unname(groupTable(getExpFG(fgXCMS))))
    expect_equal(unname(groupTable(importFeatureGroupsXCMS(XCMSImpXCMS3, getExpAnaInfo()))),
                 unname(groupTable(getExpFG(fgXCMS3))))
    expect_equal(unname(groupTable(importFeatureGroupsXCMS(XCMSImpOpenMS, getExpAnaInfo()))),
                 unname(groupTable(getExpFG(fgOpenMS))))
})

XCMS3ImpXCMS <- doExportXCMS3(fgXCMS, exportedData = FALSE)
XCMS3ImpXCMS3 <- doExportXCMS3(fgXCMS3)
XCMS3ImpOpenMS <- doExportXCMS3(fgOpenMS, exportedData = FALSE)

test_that("XCMS3 conversion", {
    expect_equal(nrow(xcms::featureDefinitions(XCMS3ImpXCMS)), length(getExpFG(fgXCMS)))
    expect_equal(nrow(xcms::featureDefinitions(XCMS3ImpXCMS3)), length(getExpFG(fgXCMS3)))
    expect_equal(nrow(xcms::featureDefinitions(XCMS3ImpOpenMS)), length(getExpFG(fgOpenMS)))
    
    expect_known_value(xcms::featureDefinitions(XCMS3ImpXCMS), testFile("fg-xcms3_import_xcms"))
    expect_known_value(xcms::featureDefinitions(XCMS3ImpXCMS3), testFile("fg-xcms3_import_xcms3"))
    expect_known_value(xcms::featureDefinitions(XCMS3ImpOpenMS), testFile("fg-xcms3_import_openms"))
    
    expect_equal(unname(groupTable(importFeatureGroupsXCMS3(XCMS3ImpXCMS, getExpAnaInfo()))),
                 unname(groupTable(getExpFG(fgXCMS))))
    expect_equal(unname(groupTable(importFeatureGroupsXCMS3(XCMS3ImpXCMS3, getExpAnaInfo()))),
                 unname(groupTable(getExpFG(fgXCMS3))))
    expect_equal(unname(groupTable(importFeatureGroupsXCMS3(XCMS3ImpOpenMS, getExpAnaInfo()))),
                 unname(groupTable(getExpFG(fgOpenMS))))
})

regr <- as.data.table(fgOpenMSConc, features = TRUE, regression = TRUE)
test_that("as.data.table works", {
    expect_equal(nrow(as.data.table(fgOpenMS)), length(fgOpenMS))

    checkmate::expect_names(names(as.data.table(fgOpenMS, average = TRUE)), must.include = replicateGroups(fgOpenMS))
    
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
                      unique(fgOpenMS, which = "standard",
                             relativeTo = setdiff(replicateGroups(fgOpenMS), "standard")))
    expect_equivalent(unique(fgOpenMS, which = "standard"), unique(fgOpenMS, which = "standard", outer = TRUE))
    expect_lt(length(unique(fgOpenMS, which = "standard")), length(fgOpenMS))
    expect_equal(length(unique(fgOpenMS, which = replicateGroups(fgOpenMS))), length(fgOpenMS))
    expect_lt(length(unique(fgOpenMS, which = replicateGroups(fgOpenMS), outer = TRUE)), length(fgOpenMS))
    expect_length(unique(fgOpenMSEmpty, which = replicateGroups(fgOpenMS)), 0)
})

test_that("overlap works", {
    # note: only have two rep groups

    expect_lt(length(overlap(fgOpenMS, which = replicateGroups(fgOpenMS))), length(fgOpenMS))
    expect_length(overlap(fgOpenMSEmpty, which = replicateGroups(fgOpenMS)), 0)
})

minInt <- function(fg, rel)
{
    # collapse to vector with use.names = FALSE: https://stackoverflow.com/a/12796124/9264518
    g <- unlist(groupTable(fg), use.names = FALSE)
    if (rel)
        return(min(g[g != 0]) / max(g))
    return(min(g[g != 0]))
}

stdRGs <- if (testWithSets()) c("standard", "standard-set2") else "standard"

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
    expect_identical(replicateGroups(filter(fgOpenMS, removeBlanks = TRUE)), stdRGs)
    expect_identical(replicateGroups(removeEmptyAnalyses(filter(fgOpenMS, relMinFeatures = 0.7))), stdRGs)
    expect_identical(replicateGroups(removeEmptyAnalyses(filter(fgOpenMS, absMinFeatures = 400))), stdRGs)
    expect_identical(replicateGroups(removeEmptyAnalyses(filter(fgOpenMS, blankThreshold = 1E6))), stdRGs)

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
    expect_setequal(names(replicateGroupSubtract(fgOpenMS, "solvent")),
                    names(unique(fgOpenMS, which = setdiff(replicateGroups(fgOpenMS), "solvent"))))
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
    expect_known_value(groupTable(fGCompOpenMS@comparedFGroups), testFile("fg-comp-openms"))
    expect_known_value(groupTable(fGCompXCMS@comparedFGroups), testFile("fg-comp-xcms"))

    expect_named(fGCompOpenMS, c("openms", "xcms"))
    expect_named(fGCompXCMS, c("openms", "xcms"))
    # expect_named(fGCompXCMS3, c("openms", "xcms"))
    expect_named(fGCompOpenMS[1], "openms")
    expect_named(fGCompOpenMS[2], "xcms")

    expect_equivalent(fGCompOpenMS[[1]], fgOpenMS)
    expect_equivalent(fGCompOpenMS[[names(fGCompOpenMS)[2]]], fgXCMS)
    expect_equivalent(fGCompOpenMS[[names(fGCompOpenMS)[2]]], fgXCMS)
    expect_equivalent(callDollar(fGCompOpenMS, names(fGCompOpenMS)[1]), fgOpenMS)

    expect_known_value(groupTable(fGCons), testFile("fg-comp-cons"))

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
    for (ana in analyses(subFGroups))
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

chordGroups <- c("standard-1" = "grp1",
                 "standard-2" = "grp2",
                 "standard-3" = "grp2",
                 "solvent-1" = "grp3",
                 "solvent-2" = "grp4",
                 "solvent-3" = "grp5")
if (testWithSets())
    chordGroups <- c(chordGroups, setNames(chordGroups, c(paste0("set2-", names(chordGroups)))))

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
    expect_doppel("chord-outer", function() plotChord(fgOpenMS, outerGroups = chordGroups))
    expect_doppel("chord-comp", function() plotChord(fGCompOpenMS))
    expect_error(plotChord(unique(fgOpenMS, which = replicateGroups(fgOpenMS), outer = TRUE),
                           average = TRUE)) # stops with nothing to plot: no overlap
    expect_plot(plotChord(unique(fgOpenMS, which = replicateGroups(fgOpenMS), outer = TRUE),
                          average = TRUE, addSelfLinks = TRUE)) # unless there are self links

    expect_doppel("eic-def", function() plotChroms(subFGroups))
    expect_doppel("eic-rtmin", function() plotChroms(subFGroups, retMin = TRUE))
    expect_doppel("eic-tm1", function() plotChroms(subFGroups, topMost = 1))
    expect_doppel("eic-area", function() plotChroms(subFGroups, showPeakArea = TRUE))
    expect_doppel("eic-cbr", function() plotChroms(subFGroups, colourBy = "rGroups"))
    expect_doppel("eic-cbf", function() plotChroms(subFGroups, colourBy = "fGroups"))
    expect_doppel("eic-ann", function() plotChroms(subFGroups, annotate = "mz"))

    expect_doppel("venn", function() plotVenn(fgOpenMS))
    # use conc fGroups as it has >2 rGroups
    expect_doppel("venn-multiple", function() plotVenn(fgOpenMSConc, which = list(standards = paste0("standard-", 1:3), solvents = "solvent")))
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

    expect_doppel("eic-def-empty", function() plotChroms(fgOpenMSEmpty))

    expect_error(plotVenn(fgOpenMSEmpty))
    expect_doppel("venn", function() plotVenn(fGConsOneEmpty)) # should be same as fgOpenMS
    expect_error(plotVenn(fgCompBothEmpty))

    expect_error(plotUpSet(fgOpenMSEmpty))
    expect_ggplot(plotUpSet(fGConsOneEmpty))
    expect_error(plotUpSet(fgCompBothEmpty))
})

if (testWithSets())
{
    fgOpenMSOneEmptySet <- makeSet(unset(fgOpenMS, "set1"), unset(fgOpenMSEmpty, "set2"), groupAlgo = "openms",
                                   adducts = NULL, labels = c("set1", "set2"))
    
    fgUn1 <- unset(fgOpenMS, "set1"); fgUn2 <- unset(fgOpenMS, "set2")
    fgUn2 <- filter(fgUn2, absMinIntensity = 1E5)
    fgOpenMSDiffSet <- makeSet(fgUn1, fgUn2, groupAlgo = "openms", adducts = NULL, labels = c("set1", "set2"))
    
    fgUn2 <- unset(fgOpenMS, "set2")
    adducts(fgUn2) <- rep("[M-H]-", length(fgUn2))
    fgOpenMSDiffAdductSet <- makeSet(fgUn1, fgUn2, groupAlgo = "openms", adducts = NULL, labels = c("set1", "set2"))
    
    fgOpenMSDiffAdduct <- fgOpenMS
    adducts(fgOpenMSDiffAdduct, reGroup = FALSE, set = "set1")[names(fgOpenMSDiffAdduct)[3]] <- "[M+K]+"
    fgOpenMSDiffAdductRG <- fgOpenMS
    adducts(fgOpenMSDiffAdductRG, reGroup = TRUE, set = "set1")[names(fgOpenMSDiffAdduct)[3]] <- "[M+K]+"
    
    fgUniqueSet2 <- unique(fgOpenMS, which = "set2", sets = TRUE)
}

test_that("sets functionality", {
    skip_if_not(testWithSets())
    
    # proper (de)neutralization
    expect_equal(mean(groupInfo(unset(fgOpenMS, "set1"))$mzs) - mean(groupInfo(fgOpenMS[, sets = "set1"])$mzs),
                 patRoon:::adductMZDelta(as.adduct("[M+H]+")))
    expect_equal(analysisInfo(unset(fgOpenMS, "set1")), getTestAnaInfoSet1())
    expect_equal(analysisInfo(fgOpenMS[, sets = "set1"])[, 1:4], getTestAnaInfoSet1())
    expect_equal(unique(annotations(fgOpenMS)$adduct), "[M+H]+")
    expect_equal(fgOpenMS, fgOpenMS[, sets = sets(fgOpenMS)])
    expect_length(fgOpenMS[, sets = character()], 0)
    expect_length(fgOpenMS[, sets = "set1"], length(fgOpenMS) - length(fgUniqueSet2))
    expect_length(unset(fgOpenMS, set = "set1"), length(fgOpenMS) - length(fgUniqueSet2))
    expect_equal(sets(filter(fgOpenMS, sets = "set1", negate = TRUE)), "set2")
    
    # can't make empty sets from fGroups atm
    expect_error(makeSet(unset(fgOpenMSEmpty, "set1"), unset(fgOpenMSEmpty, "set2"),
                         adducts = NULL, labels = c("set1", "set2")))
    
    expect_length(fgOpenMSOneEmptySet, length(fgOpenMS) - length(fgUniqueSet2))
    expect_length(filter(fgOpenMSOneEmptySet, absMinSets = 1), length(fgOpenMS) - length(fgUniqueSet2))
    expect_length(filter(fgOpenMSOneEmptySet, absMinSets = 2), 0)
    expect_length(filter(fgOpenMSOneEmptySet, relMinSets = 1), 0)
    
    expect_length(filter(fgOpenMSDiffSet, absMinSets = 1), length(fgOpenMSDiffSet))
    expect_lt(length(filter(fgOpenMSDiffSet, absMinSets = 2)), length(fgOpenMSDiffSet))
    expect_lt(length(filter(fgOpenMSDiffSet, relMinSets = 1)), length(fgOpenMSDiffSet))
    expect_setequal(adducts(fgOpenMSDiffAdductSet, set = "set2"), "[M-H]-")
    expect_gt(length(fgOpenMSDiffAdductSet), length(fgOpenMS)) # different adducts: less grouping
    
    expect_equal(adducts(fgOpenMSDiffAdduct, set = "set1")[names(fgOpenMSDiffAdduct)[3]], "[M+K]+",
                 check.attributes = FALSE)
    expect_true(any(adducts(fgOpenMSDiffAdductRG, set = "set1") == "[M+K]+"))
    expect_setequal(names(fgOpenMSDiffAdduct), names(fgOpenMS))
    # should re-group --> new group names
    expect_false(checkmate::testSubset(names(fgOpenMSDiffAdductRG), names(fgOpenMS)))
    expect_true(all(sapply(seq_len(nrow(annotations(fgOpenMSDiffAdductRG))), function(i)
    {
        ann <- annotations(fgOpenMSDiffAdductRG)[i]
        # neutral mass equals neutralized group mass
        return(isTRUE(all.equal(ann$neutralMass, groupInfo(fgOpenMSDiffAdductRG)[ann$group, "mzs"])))
    })))
    
    expect_length(filter(fgOpenMSEmpty, relMinSets = 1), 0)
    
    expect_lt(length(fgUniqueSet2), length(fgOpenMS))
    expect_equal(length(unique(fgOpenMS, which = sets(fgOpenMS), sets = TRUE)), length(fgOpenMS))
    expect_length(unique(fgOpenMSEmpty, which = sets(fgOpenMSEmpty), sets = TRUE), 0)
    
    expect_lt(length(overlap(fgOpenMSDiffSet, which = sets(fgOpenMSDiffSet), sets = TRUE)), length(fgOpenMSDiffSet))
    expect_length(overlap(fgOpenMSEmpty, which = sets(fgOpenMSDiffSet), sets = TRUE), 0)
    
    expect_doppel("venn-sets", function() plotVenn(fgOpenMS, sets = TRUE))
})
