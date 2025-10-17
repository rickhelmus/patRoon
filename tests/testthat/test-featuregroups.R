# SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

context("feature groups")

initXCMS()

fList <- getTestFeatures(noiseThrInt = 1E5)

fgOpenMS <- groupFeatures(fList, "openms")
# BUG: old XCMS obiwarp doesn' work with non sets test data
fgXCMS <- groupFeatures(fList, "xcms")
fgXCMS3 <- groupFeatures(fList, "xcms3")
fgKPIC2 <- groupFeatures(fList, "kpic2")
fgSIRIUS <- groupFeatures(analysisInfo(fList)[1,], "sirius") # only do first analysis to avoid long run times

fgOpenMSQ <- calculatePeakQualities(fgOpenMS)

fListEmpty <- getEmptyFeatures()
fgOpenMSEmpty <- groupFeatures(fListEmpty, "openms")
fgXCMSEmpty <- groupFeatures(fListEmpty, "xcms")
fgXCMS3Empty <- groupFeatures(fListEmpty, "xcms3")
fgKPIC2Empty <- groupFeatures(fListEmpty, "kpic2")
fgOpenMSEmptyQ <- calculatePeakQualities(fgOpenMSEmpty)

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
    expect_known_value(groupTable(fgKPIC2), testFile("fg-kpic2"))
    expect_known_value(groupTable(fgSIRIUS), testFile("fg-sirius"))
    expect_known_value(groupTable(fgOpenMSQ), testFile("fg-openms-qual"))
})

test_that("verify show output", {
    expect_known_show(fgOpenMS, testFile("fg-show-openms", text = TRUE))
    expect_known_show(fgXCMS, testFile("fg-show-xcms", text = TRUE))
    expect_known_show(fgXCMS3, testFile("fg-show-xcms3", text = TRUE))
    expect_known_show(fgKPIC2, testFile("fg-show-kpic2", text = TRUE))
    expect_known_show(fgSIRIUS, testFile("fg-show-sirius", text = TRUE))
    expect_known_show(fgOpenMSQ, testFile("fg-show-openms-qual", text = TRUE))
})

test_that("empty objects work", {
    expect_length(fgOpenMSEmpty, 0)
    expect_length(fgXCMSEmpty, 0)
    expect_length(fgXCMS3Empty, 0)
    expect_length(fgKPIC2Empty, 0)
    expect_length(fgOpenMSEmptyQ, 0)
})

# to compare original anaInfo: ignore extra columns that may have been added afterwards
getAnaInfo <- function(fg) analysisInfo(fg)[, -"set"]

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
    
    expect_equal(replicates(fgOpenMS[, replicates = "standard-pos"]), "standard-pos")
    expect_equal(replicates(fgOpenMS[, ni = replicate == "standard-pos"]), "standard-pos")
})

revSeq <- seq(length(analyses(fgOpenMS)), 1)
fgOpenMSRev <- fgOpenMS[revSeq, reorder = TRUE]
test_that("anaInfo modification", {
    expect_equal(fgOpenMS[revSeq, reorder = FALSE], fgOpenMS)
    expect_equal(fgOpenMSRev, { fg <- fgOpenMS; analysisInfo(fg) <- analysisInfo(fg)[revSeq, ]; fg })
    expect_equal(groupTable(fgOpenMSRev), groupTable(fgOpenMS)[revSeq])
    expect_equal(groupFeatIndex(fgOpenMSRev), groupFeatIndex(fgOpenMS)[revSeq])
})


expfile <- file.path(getWorkPath(), "export.csv") # NOTE: will be removed prior to each test automatically
test_that("exporting works", {
    expect_file(doExport(fgOpenMS, "brukerpa", expfile), expfile)
    expect_file(doExport(fgOpenMS, "brukertasq", expfile), expfile)
    expect_file(doExport(fgOpenMS, "mzmine", expfile), expfile)
    expect_error(doExport(fgOpenMSEmpty, "brukerpa", expfile))
})

XCMSImpXCMS <- doExportXCMS(fgXCMS)
XCMSImpXCMS3 <- doExportXCMS(fgXCMS3, loadRawData = FALSE)
XCMSImpOpenMS <- doExportXCMS(fgOpenMS, loadRawData = FALSE)
XCMSImpKPIC2 <- doExportXCMS(fgKPIC2, loadRawData = FALSE)
XCMSImpSIRIUS <- doExportXCMS(fgSIRIUS, loadRawData = FALSE)
test_that("XCMS conversion", {
    expect_equal(nrow(xcms::groups(XCMSImpXCMS)), length(getExpFG(fgXCMS)))
    expect_equal(nrow(xcms::groups(XCMSImpXCMS3)), length(getExpFG(fgXCMS3)))
    expect_equal(nrow(xcms::groups(XCMSImpOpenMS)), length(getExpFG(fgOpenMS)))
    expect_equal(nrow(xcms::groups(XCMSImpKPIC2)), length(getExpFG(fgKPIC2)))
    expect_equal(nrow(xcms::groups(XCMSImpSIRIUS)), length(getExpFG(fgSIRIUS)))
    
    expect_known_value(xcms::groups(XCMSImpXCMS), testFile("fg-xcms_import_xcms"))
    expect_known_value(xcms::groups(XCMSImpXCMS3), testFile("fg-xcms_import_xcms3"))
    expect_known_value(xcms::groups(XCMSImpOpenMS), testFile("fg-xcms_import_openms"))
    expect_known_value(xcms::groups(XCMSImpKPIC2), testFile("fg-xcms_import_kpic2"))
    expect_known_value(xcms::groups(XCMSImpSIRIUS), testFile("fg-xcms_import_sirius"))
    
    expect_equal(unname(groupTable(importFeatureGroups(XCMSImpXCMS, "xcms", getExpAnaInfo()))),
                 unname(groupTable(getExpFG(fgXCMS))))
    expect_equal(unname(groupTable(importFeatureGroups(XCMSImpXCMS3, "xcms", getExpAnaInfo()))),
                 unname(groupTable(getExpFG(fgXCMS3))))
    expect_equal(unname(groupTable(importFeatureGroups(XCMSImpOpenMS, "xcms", getExpAnaInfo()))),
                 unname(groupTable(getExpFG(fgOpenMS))))
    expect_equal(unname(groupTable(importFeatureGroups(XCMSImpKPIC2, "xcms", getExpAnaInfo()))),
                 unname(groupTable(getExpFG(fgKPIC2))))
    expect_equal(unname(groupTable(importFeatureGroups(XCMSImpSIRIUS, "xcms", getExpAnaInfo()[1, ]))),
                 unname(groupTable(getExpFG(fgSIRIUS))))
})

XCMS3ImpXCMS <- doExportXCMS3(fgXCMS, loadRawData = FALSE)
XCMS3ImpXCMS3 <- doExportXCMS3(fgXCMS3)
XCMS3ImpOpenMS <- doExportXCMS3(fgOpenMS, loadRawData = FALSE)
XCMS3ImpKPIC2 <- doExportXCMS3(fgKPIC2, loadRawData = FALSE)
XCMS3ImpSIRIUS <- doExportXCMS3(fgSIRIUS, loadRawData = FALSE)

test_that("XCMS3 conversion", {
    expect_equal(nrow(xcms::featureDefinitions(XCMS3ImpXCMS)), length(getExpFG(fgXCMS)))
    expect_equal(nrow(xcms::featureDefinitions(XCMS3ImpXCMS3)), length(getExpFG(fgXCMS3)))
    expect_equal(nrow(xcms::featureDefinitions(XCMS3ImpOpenMS)), length(getExpFG(fgOpenMS)))
    expect_equal(nrow(xcms::featureDefinitions(XCMS3ImpKPIC2)), length(getExpFG(fgKPIC2)))
    expect_equal(nrow(xcms::featureDefinitions(XCMS3ImpSIRIUS)), length(getExpFG(fgSIRIUS)))
    
    expect_known_value(xcms::featureDefinitions(XCMS3ImpXCMS), testFile("fg-xcms3_import_xcms"))
    expect_known_value(xcms::featureDefinitions(XCMS3ImpXCMS3), testFile("fg-xcms3_import_xcms3"))
    expect_known_value(xcms::featureDefinitions(XCMS3ImpOpenMS), testFile("fg-xcms3_import_openms"))
    expect_known_value(xcms::featureDefinitions(XCMS3ImpKPIC2), testFile("fg-xcms3_import_kpic2"))
    expect_known_value(xcms::featureDefinitions(XCMS3ImpSIRIUS)[, names(xcms::featureDefinitions(XCMS3ImpSIRIUS)) != "peakidx"],
                       testFile("fg-xcms3_import_sirius")) # NOTE: peakidx not consistent
    
    expect_equal(unname(groupTable(importFeatureGroups(XCMS3ImpXCMS, "xcms3", getExpAnaInfo()))),
                 unname(groupTable(getExpFG(fgXCMS))))
    expect_equal(unname(groupTable(importFeatureGroups(XCMS3ImpXCMS3, "xcms3", getExpAnaInfo()))),
                 unname(groupTable(getExpFG(fgXCMS3))))
    expect_equal(unname(groupTable(importFeatureGroups(XCMS3ImpOpenMS, "xcms3", getExpAnaInfo()))),
                 unname(groupTable(getExpFG(fgOpenMS))))
    expect_equal(unname(groupTable(importFeatureGroups(XCMS3ImpKPIC2, "xcms3", getExpAnaInfo()))),
                 unname(groupTable(getExpFG(fgKPIC2))))
    expect_equal(unname(groupTable(importFeatureGroups(XCMS3ImpSIRIUS, "xcms3", getExpAnaInfo()[1, ]))),
                 unname(groupTable(getExpFG(fgSIRIUS))))
})

fGroupsRegr <- fgOpenMS
analysisInfo(fGroupsRegr)$conc <- seq_along(analyses(fGroupsRegr))
fGroupsRegrNA <- fGroupsRegr
analysisInfo(fGroupsRegrNA)$conc[c(TRUE, FALSE)] <- NA
regr <- as.data.table(fGroupsRegr, regression = TRUE)
regrF <- as.data.table(fGroupsRegr, regression = TRUE, features = TRUE)
regrFBA <- as.data.table(fGroupsRegr, regression = TRUE, features = TRUE, regressionBy = "set", average = TRUE)
regrFNA <- as.data.table(fGroupsRegrNA, regression = TRUE, features = TRUE, anaInfoCols = "conc")

FCParams <- getFCParams(c("solvent-pos", "standard-pos"))
fctbl <- as.data.table(fgOpenMS, FCParams = FCParams)
fctblF <- as.data.table(fgOpenMS, FCParams = FCParams, features = TRUE)
test_that("as.data.table works", {
    expect_equal(nrow(as.data.table(fgOpenMS)), length(fgOpenMS))

    checkmate::expect_names(names(as.data.table(fgOpenMS, average = TRUE)),
                            must.include = getADTIntCols(replicates(fgOpenMS)))
    checkmate::expect_names(names(as.data.table(fgOpenMS, average = "set")),
                            must.include = getADTIntCols(sets(fgOpenMS)))
    checkmate::expect_names(names(as.data.table(fgOpenMS, average = "fGroups")), must.include = "intensity")
    
    
    # UNDONE: intensities are sometimes higher than areas?
    # expect_gt_or_zero(as.data.table(fgOpenMS, areas = TRUE), as.data.table(fgOpenMS, areas = FALSE))
    # check if area from first group of first analysis corresponds to its feature data
    expect_equal(as.data.table(fgOpenMS, areas = TRUE)[[getADTIntCols(analyses(fgOpenMS)[1])]][11],
                 featureTable(fgOpenMS)[[analyses(fgOpenMS)[1]]][["area"]][groupFeatIndex(fgOpenMS)[[c(11, 1)]]])
    
    expect_equal(nrow(as.data.table(fgOpenMS, features = TRUE)), length(getFeatures(fgOpenMS)))
    expect_equal(nrow(as.data.table(fgOpenMS, features = TRUE, average = "analysis")), length(getFeatures(fgOpenMS)))
    expect_equal(nrow(as.data.table(fgOpenMS, features = TRUE, average = "fGroups")), length(fgOpenMS))
    expect_setequal(as.data.table(fgOpenMS, features = TRUE, average = "set")$average_group, sets(fgOpenMS))

    checkmate::expect_names(names(as.data.table(fGroupsRegr, features = TRUE, anaInfoCols = "conc")),
                            must.include = "anaInfo_conc")
    
    expect_error(as.data.table(fgOpenMS, regression = TRUE)) # no conc specified
    expect_equal(regr, as.data.table(fGroupsRegr, regression = "conc"))
    expect_equal(regrF, as.data.table(fGroupsRegr, features = TRUE, regression = "conc"))
    testRegrTab(fGroupsRegr, FALSE, FALSE, FALSE)
    testRegrTab(fGroupsRegr, FALSE, TRUE, FALSE)
    testRegrTab(fGroupsRegr, FALSE, TRUE, TRUE)
    testRegrTab(fGroupsRegr, TRUE, FALSE, FALSE)
    testRegrTab(fGroupsRegr, TRUE, TRUE, FALSE)
    testRegrTab(fGroupsRegr, TRUE, TRUE, TRUE)
    expect_true(all(is.na(regrF$x_reg) | is.na(regrF$x_reg) | is.na(regrF$RSQ) | regrF$RSQ < 0.9 |
                        abs(regrF$x_reg - regrF$anaInfo_conc) < 0.5)) # calculated x values should be somewhat close
    expect_true(any(is.na(regrFNA$anaInfo_conc) & !is.na(regrFNA$x_reg))) # at least some x estimations should be there
    # regression properties should be equal per fGroup+regression_group
    expect_equal(uniqueN(regrFBA[, .(group, regression_group, slope, intercept, RSQ, p)]),
                 uniqueN(regrFBA[, .(group, regression_group)]))

    checkmate::expect_names(names(fctbl), must.include = c("FC", "FC_log", "PV", "PV_log", "classification"))
    checkmate::expect_subset(fctbl$classification, c("insignificant", "FC", "increase", "decrease", "significant"))
    checkmate::expect_names(names(fctblF), must.include = c("FC", "FC_log", "PV", "PV_log", "classification"))
    checkmate::expect_subset(fctblF$classification, c("insignificant", "FC", "increase", "decrease", "significant"))
    
    expect_identical(as.data.table(fgOpenMS), as.data.table(fgOpenMSQ, qualities = FALSE)) # nothing extra reported
    expect_identical(as.data.table(fgOpenMS), as.data.table(fgOpenMS, qualities = "both")) # nothing to report
    checkmate::expect_names(names(as.data.table(fgOpenMSQ, qualities = "quality")), must.include = featureQualityNames())
    checkmate::expect_names(names(as.data.table(fgOpenMSQ, qualities = "score")),
                            must.include = featureQualityNames(scores = TRUE))
    checkmate::expect_names(names(as.data.table(fgOpenMSQ, qualities = "both")),
                            must.include = c(featureQualityNames(), featureQualityNames(scores = TRUE)))
    
    expect_equal(nrow(as.data.table(fgOpenMSEmpty, average = TRUE)), 0)
    expect_equal(nrow(as.data.table(fgOpenMSEmpty, features = TRUE)), 0)
    expect_equal(nrow(as.data.table(fgOpenMSEmpty, average = TRUE, features = TRUE)), 0)
    expect_equal(nrow(as.data.table(fgOpenMSEmpty, FCParams = FCParams)), 0)
    expect_equal(nrow(as.data.table(fgOpenMSEmptyQ, qualities = "both")), 0)
})

test_that("unique works", {
    # note: only have two rep groups

    expect_equal(unique(fgOpenMS, aggregate = FALSE,
                        which = analysisInfo(fgOpenMS)[replicate == "standard-pos"]$analysis),
                 unique(fgOpenMS, aggregate = TRUE, which = "standard-pos"))
    expect_equivalent(unique(fgOpenMS, which = "standard-pos"),
                      unique(fgOpenMS, which = "standard-pos",
                             relativeTo = setdiff(replicates(fgOpenMS), "standard-pos")))
    expect_equivalent(unique(fgOpenMS, which = "standard-pos"), unique(fgOpenMS, which = "standard-pos", outer = TRUE))
    expect_lt(length(unique(fgOpenMS, which = "standard-pos")), length(fgOpenMS))
    expect_equal(length(unique(fgOpenMS, which = replicates(fgOpenMS))), length(fgOpenMS))
    expect_lt(length(unique(fgOpenMS, which = replicates(fgOpenMS), outer = TRUE)), length(fgOpenMS))
    expect_length(unique(fgOpenMSEmpty, which = replicates(fgOpenMS)), 0)
})

test_that("overlap works", {
    # note: only have two rep groups

    expect_equal(overlap(fgOpenMS, aggregate = FALSE, which = analyses(fgOpenMS)),
                 overlap(fgOpenMS, aggregate = TRUE, which = replicates(fgOpenMS)))
    expect_lt(length(overlap(fgOpenMS, which = replicates(fgOpenMS))), length(fgOpenMS))
    expect_length(overlap(fgOpenMSEmpty, which = replicates(fgOpenMS)), 0)
})

minInt <- function(fg, rel)
{
    # collapse to vector with use.names = FALSE: https://stackoverflow.com/a/12796124/9264518
    g <- unlist(groupTable(fg), use.names = FALSE)
    if (rel)
        return(min(g[g != 0]) / max(g))
    return(min(g[g != 0]))
}
featCounts <- function(fg, rel)
{
    counts <- unlist(as.data.table(fg)[, lapply(.SD, function(x) sum(x > 0)), .SDcols = getADTIntCols(analyses(fg))])
    return(if (rel) counts / length(fg) else counts)
}

stdReps <- c("standard-pos", "standard-neg")

qr <- list(ZigZag = c(0.2, 0.6), TPASRScore = c(0.5, 0.9))
fgOpenMSQFF <- filter(fgOpenMSQ, featQualityRange = qr)
fgOpenMSQFFTab <- as.data.table(fgOpenMSQFF, features = TRUE, qualities = "both")
fgOpenMSQFG <- filter(fgOpenMSQ, groupQualityRange = qr)
fgOpenMSQFFN <- filter(fgOpenMSQ, featQualityRange = qr, negate = TRUE)
fgOpenMSQFFTabN <- as.data.table(fgOpenMSQFFN, features = TRUE, qualities = "both")
fgOpenMSQFGN <- filter(fgOpenMSQ, groupQualityRange = qr, negate = TRUE)

test_that("delete and filter", {
    checkmate::expect_names(analyses(delete(fgOpenMS, i = 1)), disjunct.from = analyses(fgOpenMS)[1])
    checkmate::expect_names(analyses(delete(fgOpenMS, i = analyses(fgOpenMS)[1])), disjunct.from = analyses(fgOpenMS)[1])
    expect_length(delete(fgOpenMS, i = analyses(fgOpenMS)), 0)
    checkmate::expect_names(names(delete(fgOpenMS, j = 1:5)), disjunct.from = names(fgOpenMS)[1:5])
    checkmate::expect_names(names(delete(fgOpenMS, j = names(fgOpenMS)[1:5])), disjunct.from = names(fgOpenMS)[1:5])
    expect_true(all(groupTable(delete(fgOpenMS, j = function(...) 3))[3, ] == 0))
    expect_length(delete(fgOpenMS, j = function(...) TRUE), 0)
    expect_equal(delete(fgOpenMS, i = character()), fgOpenMS)
    expect_equal(delete(fgOpenMS, j = integer()), fgOpenMS)
    expect_length(delete(fgOpenMS), 0)
    expect_length(delete(fgOpenMSEmpty), 0)
    
    expect_gte(minInt(filter(fgOpenMS, absMinIntensity = 1500), FALSE), 1500)
    expect_gte(minInt(filter(fgOpenMS, relMinIntensity = 0.2), TRUE), 0.2)
    expect_gte(minInt(filter(fgOpenMS, preAbsMinIntensity = 1500), FALSE), 1500)
    expect_gte(minInt(filter(fgOpenMS, preRelMinIntensity = 0.2), TRUE), 0.2)
    expect_true(all(sapply(groupTable(filter(fgOpenMS, absMinMaxIntensity = 1E6)), function(x) max(x) >= 1E6)))
    expect_true(all(sapply(groupTable(filter(fgOpenMS, relMinMaxIntensity = 0.5)),
                           function(x) max(x) >= (0.5 * max(groupTable(fgOpenMS))))))

    expect_range(groupInfo(filter(fgOpenMS, retentionRange = c(120, 200)))$ret, c(120, 200))
    # expect_equivalent(filter(fgOpenMS, retentionRange = c(0, Inf)), fgOpenMS)
    expect_equivalent(filter(fgXCMS3, retentionRange = c(0, Inf)), fgXCMS3) # NOTE: cannot use OpenMS as it may  yield negative RTs...
    expect_range(groupInfo(filter(fgOpenMS, mzRange = c(200, 300)))$mz, c(200, 300))
    expect_equivalent(filter(fgOpenMS, mzRange = c(0, Inf)), fgOpenMS)
    expect_range(groupInfo(filter(fgOpenMS, mzDefectRange = c(0.1, 0.2)))$mz %% 1, c(0.1, 0.2))
    expect_equivalent(filter(fgOpenMS, mzDefectRange = c(0, 1)), fgOpenMS)
    expect_lt(length(filter(fgOpenMS, chromWidthRange = c(0, 30))), length(fgOpenMS))
    expect_equivalent(filter(fgOpenMS, chromWidthRange = c(0, Inf)), fgOpenMS)

    expect_identical(replicates(filter(fgOpenMS, replicates = "standard-pos")), "standard-pos")
    expect_identical(replicates(fgOpenMS[, replicates = "standard-pos"]), "standard-pos")
    expect_identical(replicates(filter(fgOpenMS, removeBlanks = TRUE)), stdReps)
    expect_identical(replicates(removeEmptyAnalyses(filter(fgOpenMS, blankThreshold = 1E6))), stdReps)
    
    expect_gte(min(featCounts(filter(fgOpenMS, relMinFeatures = 0.3), TRUE)), 0.3)
    expect_gte(min(featCounts(filter(fgOpenMS, absMinFeatures = length(fgOpenMS) * 0.3), FALSE)),
               length(fgOpenMS) * 0.3)
    
    expect_range(fgOpenMSQFFTab[[names(qr)[1]]], qr[[1]])
    expect_range(fgOpenMSQFFTab[[names(qr)[2]]], qr[[2]])
    expect_true(all(fgOpenMSQFFTabN[[names(qr)[1]]] < qr[[c(1, 1)]] | fgOpenMSQFFTabN[[names(qr)[1]]] > qr[[c(1, 2)]]))
    expect_true(all(fgOpenMSQFFTabN[[names(qr)[2]]] < qr[[c(2, 1)]] | fgOpenMSQFFTabN[[names(qr)[2]]] > qr[[c(2, 2)]]))
    
    expect_range(groupQualities(fgOpenMSQFG)[[names(qr)[1]]], qr[[1]])
    expect_range(groupScores(fgOpenMSQFG)[[names(qr)[2]]], qr[[2]])
    expect_true(all(groupQualities(fgOpenMSQFGN)[[names(qr)[1]]] < qr[[c(1, 1)]] |
                        groupQualities(fgOpenMSQFGN)[[names(qr)[1]]] > qr[[c(1, 2)]]))
    expect_true(all(groupScores(fgOpenMSQFGN)[[names(qr)[2]]] < qr[[c(2, 1)]] |
                        groupScores(fgOpenMSQFGN)[[names(qr)[2]]] > qr[[c(2, 2)]]))

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

test_that("replicate subtraction", {
    expect_setequal(names(replicateSubtract(fgOpenMS, "solvent-pos")),
                    names(unique(fgOpenMS, which = setdiff(replicates(fgOpenMS), "solvent-pos"))))
    expect_length(replicateSubtract(fgOpenMSEmpty, "solvent-pos"), 0)
})

fgISTD <- fgOpenMS
fgISTD@features@analysisInfo <- copy(fgISTD@features@analysisInfo)
fgISTD@features@analysisInfo[, norm_conc := rep(c(NA, NA, NA, 1, 2, 1), length.out = nrow(analysisInfo(fgISTD)))]
fgNormISTDMin1 <- doNormInts(fgISTD, "istd", ISTDRTWindow = 120, ISTDMZWindow = 300, minISTDs = 1)
fgNormISTDMin2 <- doNormInts(fgISTD, "istd", ISTDRTWindow = 120, ISTDMZWindow = 300, minISTDs = 2)
fgNormISTDEmpty <- doNormInts(fgOpenMSEmpty, "istd")
fgNormTIC <- doNormInts(fgISTD, "tic")
fgNormConc <- doNormInts(fgISTD, "conc")
fgNormGroup <- doNormInts(fgISTD, groupNorm = TRUE)

checkISTDNormWindows <- function(fg, minISTD, RTWin, MZWin)
{
    ia <- getISTDAssignments(fg)
    ia <- ia[lengths(ia) > minISTD]
    gInfo <- groupInfo(fg)
    calcDev <- function(grp, istdGrps, what) abs(gInfo[group == grp][[what]] - gInfo[match(istdGrps, group)][[what]])
    maxRTDev <- max(unlist(Map(names(ia), ia, f = calcDev, MoreArgs = list(what = "ret"))))
    maxMZDev <- max(unlist(Map(names(ia), ia, f = calcDev, MoreArgs = list(what = "mz"))))
    expect_lte(maxRTDev, RTWin)
    expect_lte(maxMZDev, MZWin)
}

test_that("Normalization", {
    # first three analyses are without IS, should give zero intensities
    expect_true(all(groupTable(fgNormISTDMin1[1:3], normalized = TRUE) == 0))
    expect_true(all(groupTable(fgNormConc[1:3], normalized = TRUE) == 0))
    
    expect_true(all(lengths(getISTDAssignments(fgNormISTDMin1)) >= 1))
    expect_true(all(lengths(getISTDAssignments(fgNormISTDMin2)) >= 2))
    
    checkISTDNormWindows(fgNormISTDMin1, 1, 120, 300)
    
    expect_equal(groupTable(fgNormTIC, normalized = FALSE), groupTable(fgISTD))
    expect_range(groupTable(fgNormTIC, normalized = TRUE), c(0, 1))
    expect_range(groupTable(fgNormGroup, normalized = TRUE), c(0, 1))
    # sample 5 in fgNormConc has halved intensities
    expect_equal(groupTable(fgNormConc[5], normalized = TRUE) * 2, groupTable(fgISTD[5]))
    
    expect_length(fgNormISTDEmpty, 0)
    expect_length(doNormInts(fgOpenMSEmpty, "tic"), 0)
    expect_length(doNormInts(fgOpenMSEmpty, "conc"), 0)
    expect_length(doNormInts(fgOpenMSEmpty, groupNorm = TRUE), 0)
})

fGCompOpenMS <- comparison(openms = fgOpenMS, xcms = fgXCMS, groupAlgo = "openms")
fGCompXCMS <- comparison(openms = fgOpenMS, xcms = fgXCMS, groupAlgo = "xcms")
fGCompXCMS3 <- comparison(openms = fgOpenMS, xcms = fgXCMS, groupAlgo = "xcms3")
fGCompKPIC2 <- comparison(openms = fgOpenMS, xcms = fgXCMS, groupAlgo = "kpic2")
fgCompOneEmpty <- comparison(openms = fgOpenMS, xcms = fgXCMSEmpty, groupAlgo = "openms")
fgCompBothEmpty <- comparison(openms = fgOpenMSEmpty, xcms = fgXCMSEmpty, groupAlgo = "openms")

test_that("verify feature group comparison", {
    expect_known_value(groupTable(fGCompOpenMS@comparedFGroups), testFile("fg-comp-openms"))
    expect_known_value(groupTable(fGCompXCMS@comparedFGroups), testFile("fg-comp-xcms"))
    expect_known_value(groupTable(fGCompXCMS3@comparedFGroups), testFile("fg-comp-xcms3"))
    expect_known_value(groupTable(fGCompKPIC2@comparedFGroups), testFile("fg-comp-kpic2"))

    expect_named(fGCompOpenMS, c("openms", "xcms"))
    expect_named(fGCompXCMS, c("openms", "xcms"))
    expect_named(fGCompXCMS3, c("openms", "xcms"))
    expect_named(fGCompKPIC2, c("openms", "xcms"))
    expect_named(fGCompOpenMS[1], "openms")
    expect_named(fGCompOpenMS[2], "xcms")

    expect_equivalent(fGCompOpenMS[[1]], fgOpenMS)
    expect_equivalent(fGCompOpenMS[[names(fGCompOpenMS)[2]]], fgXCMS)
    expect_equivalent(fGCompOpenMS[[names(fGCompOpenMS)[2]]], fgXCMS)
    expect_equivalent(callDollar(fGCompOpenMS, names(fGCompOpenMS)[1]), fgOpenMS)

    expect_length(fgCompOneEmpty, 2)
    expect_length(fgCompBothEmpty, 2)
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
})

test_that("reporting with empty object works", {
    expect_error(reportCSV(fgOpenMSEmpty, getWorkPath(), reportFeatures = TRUE), NA)
    expect_error(reportPDF(fgOpenMSEmpty, getWorkPath()), NA)
    expect_error(makeReportHTML(fgOpenMSEmpty), NA)
})

test_that("plotting works", {
    expect_doppel("retmz", function() plot(fgOpenMS, groupBy = "fGroups", showLegend = FALSE))
    expect_doppel("retmz-singlec", function() plot(fgOpenMS, groupBy = NULL, col = "blue"))
    expect_doppel("retmz-rep", function() plot(fgOpenMS, groupBy = "replicate"))
    expect_doppel("retmz-comp", function() plot(fGCompOpenMS, groupBy = "fGroups", showLegend = FALSE))

    expect_doppel("intensity-def", function() plotInt(fgOpenMS))
    expect_doppel("intensity-avg", function() plotInt(fgOpenMS, average = TRUE))
    expect_doppel("intensity-area", function() plotInt(fgOpenMS, areas = TRUE))
    expect_doppel("intensity-norm", function() plotInt(fgOpenMS, normalized = TRUE))
    expect_doppel("intensity-group_rep", function() plotInt(fgOpenMS, groupBy = "replicate"))
    expect_doppel("intensity-xBy", function() plotInt(fGroupsRegr, xBy = "conc"))
    expect_doppel("intensity-regr", function() plotInt(fGroupsRegr, regression = TRUE, xBy = "conc"))
    expect_doppel("intensity-regr-group", function() plotInt(fGroupsRegr, regression = TRUE, xBy = "conc", groupBy = "replicate"))
    expect_doppel("intensity-regr-group-leg", function() plotInt(fGroupsRegr[, 1:25], regression = TRUE, xBy = "conc",
                                                                 groupBy = "replicate", showLegend = TRUE))
    expect_doppel("intensity-regr-group-leg1", function() plotInt(fGroupsRegr[, 1], regression = TRUE, xBy = "conc",
                                                                 groupBy = "replicate", showLegend = TRUE))
    
    expect_doppel("chord-def", function() plotChord(fgOpenMS))
    expect_doppel("chord-selflinks", function() plotChord(fgOpenMS, addSelfLinks = TRUE))
    expect_doppel("chord-nortmz", function() plotChord(fgOpenMS, addRetMzPlots = FALSE))
    expect_doppel("chord-aggr", function() plotChord(fgOpenMS, aggregate = TRUE))
    expect_doppel("chord-outer", function() plotChord(fgOpenMS, groupBy = "replicate"))
    expect_doppel("chord-comp", function() plotChord(fGCompOpenMS))
    expect_error(plotChord(unique(fgOpenMS, which = replicates(fgOpenMS), outer = TRUE),
                           aggregate = TRUE)) # stops with nothing to plot: no overlap
    expect_plot(plotChord(unique(fgOpenMS, which = replicates(fgOpenMS), outer = TRUE),
                          aggregate = TRUE, addSelfLinks = TRUE)) # unless there are self links

    expect_doppel("eic-def", function() plotChroms(subFGroups))
    expect_doppel("eic-rtmin", function() plotChroms(subFGroups, retMin = TRUE))
    expect_doppel("eic-tm1", function() plotChroms(subFGroups, EICParams = getDefEICParams(topMost = 1)))
    expect_doppel("eic-tm1rg", function() plotChroms(subFGroups,
                                                     EICParams = getDefEICParams(topMost = 1, topMostByReplicate = TRUE)))
    expect_doppel("eic-area", function() plotChroms(subFGroups, showPeakArea = TRUE))
    expect_doppel("eic-gbr", function() plotChroms(subFGroups, groupBy = "replicate"))
    expect_doppel("eic-gbf", function() plotChroms(subFGroups, groupBy = "fGroups"))
    expect_doppel("eic-ann", function() plotChroms(subFGroups, annotate = "mz"))
    # below two should be mostly the same, but xlim and group rect will be slightly different since subsetting removes
    # some of the feature data that is used to determine the limits for these. For now just compare the two figures
    # manually.
    expect_doppel("eic-sub", function() plotChroms(subFGroups[1, 1]))
    expect_doppel("eic-sub2", function() plotChroms(subFGroups, analysis = analyses(subFGroups)[1],
                                                    groupName = names(subFGroups)[1]))

    expect_doppel("eic-3d-def", function() plotChroms3D(subFGroups[1, 1]))
    expect_doppel("eic-3d-rtmin", function() plotChroms3D(subFGroups[1, 1], retMin = TRUE))
    expect_doppel("eic-3d-nolim", function() plotChroms3D(subFGroups[1, 1], showLimits = FALSE))

    expect_doppel("venn", function() plotVenn(fgOpenMS))
    # use conc fGroups as it has >2 replicates
    expect_doppel("venn-multiple", function() plotVenn(fgOpenMS, which = c("standard-pos", "solvent-pos"),
                                                       aggregate = TRUE))
    expect_doppel("venn-comp", function() plotVenn(fGCompOpenMS))
    expect_equal(expect_plot(plotVenn(fgOpenMS, which = c("solvent-pos", "standard-pos")))$areas[2],
                 length(filter(fgOpenMS, replicates = "standard-pos")))
    expect_equal(expect_plot(plotVenn(fGCompOpenMS))$areas[2], length(fgXCMS))

    # vdiffr doesn't work with UpSet
    expect_ggplot(plotUpSet(fgOpenMS, aggregate = FALSE))
    expect_ggplot(plotUpSet(fgOpenMS, aggregate = TRUE))
    expect_ggplot(plotUpSet(fgOpenMS, which = c("standard-pos", "solvent-pos"), aggregate = TRUE))
    expect_ggplot(plotUpSet(fGCompOpenMS))
    
    expect_doppel("volcano", function() plotVolcano(fgOpenMS, FCParams))
})

test_that("plotting empty objects works", {
    expect_doppel("retmz-empty", function() plot(fgOpenMSEmpty))
    expect_doppel("retmz-empty", function() plot(fgOpenMSEmpty, groupBy = "replicate"))
    expect_doppel("retmz-comp-empty", function() plot(fgCompBothEmpty))

    expect_doppel("intensity-def-empty", function() plotInt(fgOpenMSEmpty))
    expect_doppel("intensity-avg-empty", function() plotInt(fgOpenMSEmpty, TRUE))

    expect_error(plotChord(fgOpenMSEmpty))
    expect_error(plotChord(fgCompBothEmpty))

    expect_doppel("eic-def-empty", function() plotChroms(fgOpenMSEmpty))
    expect_doppel("eic-3d-def-empty", function() plotChroms3D(fgOpenMSEmpty[1]))

    expect_error(plotVenn(fgOpenMSEmpty))
    expect_error(plotVenn(fgCompBothEmpty))

    expect_error(plotUpSet(fgOpenMSEmpty))
    expect_error(plotUpSet(fgCompBothEmpty))
})

fgOpenMSOneEmptySet <- makeSet(unset(fgOpenMS, "positive"), unset(fgOpenMSEmpty, "negative"), groupAlgo = "openms",
                               adducts = NULL, labels = c("positive", "negative"))

fgOpenMSDiffAdduct <- fgOpenMS
adducts(fgOpenMSDiffAdduct, reGroup = FALSE, set = "positive")[names(fgOpenMSDiffAdduct[, sets = "positive"])[3]] <- "[M+K]+"
fgOpenMSDiffAdductRG <- fgOpenMS
adducts(fgOpenMSDiffAdductRG, reGroup = TRUE, set = "positive")[names(fgOpenMSDiffAdductRG[, sets = "positive"])] <- "[M+K]+"

fgUniqueSet2 <- unique(fgOpenMS, which = "negative", aggregate = "set")

test_that("sets functionality", {
    # proper (de)neutralization
    expect_equal(patRoon:::calculateMasses(groupInfo(unset(fgOpenMS, "positive"))$mz, as.adduct("[M+H]+"), "neutral"),
                 groupInfo(fgOpenMS[, sets = "positive"])$mz)
    expect_equal(analysisInfo(unset(fgOpenMS, "positive"), TRUE), getTestAnaInfoPos())
    expect_equivalent(analysisInfo(fgOpenMS[, sets = "positive"], FALSE)[, -"set"], getTestAnaInfoPos())
    expect_setequal(annotations(fgOpenMS)$adduct, c("[M+H]+", "[M-H]-"))
    expect_equal(fgOpenMS, fgOpenMS[, sets = sets(fgOpenMS)])
    expect_equal(sets(fgOpenMS[, sets = "positive"]), "positive")
    expect_equal(sets(fgOpenMS[, ni = set == "positive"]), "positive")
    expect_length(fgOpenMS[, sets = character()], 0)
    expect_length(fgOpenMS[, sets = "positive"], length(fgOpenMS) - length(fgUniqueSet2))
    expect_equal(sets(fgOpenMS[, sets = rev(sets(fgOpenMS)), reorder = TRUE]), rev(sets(fgOpenMS)))
    expect_length(unset(fgOpenMS, set = "positive"), length(fgOpenMS) - length(fgUniqueSet2))
    expect_equal(sets(filter(fgOpenMS, sets = "positive", negate = TRUE)), "negative")
    
    # can't make empty sets from fGroups atm
    expect_error(makeSet(unset(fgOpenMSEmpty, "positive"), unset(fgOpenMSEmpty, "negative"),
                         adducts = NULL, labels = c("positive", "negative")))
    
    expect_length(fgOpenMSOneEmptySet, length(fgOpenMS) - length(fgUniqueSet2))
    expect_length(filter(fgOpenMSOneEmptySet, absMinSets = 1), length(fgOpenMS) - length(fgUniqueSet2))
    expect_length(filter(fgOpenMSOneEmptySet, absMinSets = 2), 0)
    expect_length(filter(fgOpenMSOneEmptySet, relMinSets = 1), 0)
    
    expect_length(filter(fgOpenMS, absMinSets = 1), length(fgOpenMS))
    expect_lt(length(filter(fgOpenMS, absMinSets = 2)), length(fgOpenMS))
    expect_lt(length(filter(fgOpenMS, relMinSets = 1)), length(fgOpenMS))
    expect_setequal(adducts(fgOpenMS, set = "negative"), "[M-H]-")
    
    expect_equal(fgOpenMSDiffAdduct@annotationsChanged[set == "positive" & group == names(fgOpenMSDiffAdduct[, sets = "positive"])[3]]$adduct, "[M+K]+")
    # no grouping, annotations should be the same
    expect_equal(annotations(fgOpenMSDiffAdduct), annotations(fgOpenMS))
    expect_false(length(fgOpenMSDiffAdductRG) == length(fgOpenMS)) # different adducts: different grouping
    expect_setequal(adducts(fgOpenMSDiffAdductRG, set = "positive"), "[M+K]+")
    expect_setequal(names(fgOpenMSDiffAdduct), names(fgOpenMS))
    # should re-group --> new group names
    expect_false(checkmate::testSubset(names(fgOpenMSDiffAdductRG), names(fgOpenMS)))
    expect_true(all(sapply(seq_len(nrow(annotations(fgOpenMSDiffAdductRG))), function(i)
    {
        ann <- annotations(fgOpenMSDiffAdductRG)[i]
        # neutral mass equals neutralized group mass
        return(isTRUE(all.equal(ann$neutralMass, groupInfo(fgOpenMSDiffAdductRG)[match(ann$group, group)]$mz)))
    })))
    
    expect_length(filter(fgOpenMSEmpty, relMinSets = 1), 0)
    
    expect_lt(length(fgUniqueSet2), length(fgOpenMS))
    expect_equal(length(unique(fgOpenMS, which = sets(fgOpenMS), aggregate = "set")), length(fgOpenMS))
    expect_length(unique(fgOpenMSEmpty, which = sets(fgOpenMSEmpty), aggregate = "set"), 0)
    
    expect_lt(length(overlap(fgOpenMS, which = sets(fgOpenMS), aggregate = "set")), length(fgOpenMS))
    expect_length(overlap(fgOpenMSEmpty, which = sets(fgOpenMS), aggregate = "set"), 0)
    
    expect_doppel("intensity-avg_set", function() plotInt(fgOpenMS, average = "set"))
    expect_doppel("intensity-group_set", function() plotInt(fgOpenMS, groupBy = "set"))
    expect_doppel("intensity-regr-set", function() plotInt(fGroupsRegr[, 1:25], regression = TRUE, xBy = "conc",
                                                           groupBy = "set", showLegend = TRUE))
    
    expect_doppel("venn-sets", function() plotVenn(fgOpenMS, aggregate = "set"))
    
    expect_ggplot(plotUpSet(fgOpenMS, aggregate = "set"))
    
    expect_doppel("chord-aggr-set", function() plotChord(fgOpenMS, aggregate = "set"))
    expect_doppel("chord-outer-set", function() plotChord(fgOpenMS, aggregate = TRUE, groupBy = "set"))
    
    expect_doppel("eic-gby_set", function() plotChroms(subFGroups, groupBy = "set"))
    
    expect_HTML(plotGraph(fgNormISTDMin1, onlyPresent = FALSE, set = "positive"))
    expect_HTML(plotGraph(fgNormISTDMin1, onlyPresent = TRUE, set = "positive"))
    expect_HTML(plotGraph(fgNormISTDEmpty, set = "positive"))
})

fListUS <- unset(fList, "positive")
fgNS <- groupFeatures(fListUS, "openms")
fgNSAnn <- fgNS
adducts(fgNSAnn) <- rep("[M+H]+", length(fgNSAnn))
fgNSAnn2 <- fgNSAnn
adducts(fgNSAnn2)[3] <- "[M+K]+"

fGCompNS <- comparison(openms = fgNS, xcms = groupFeatures(fListUS, "xcms", rtalign = FALSE), groupAlgo = "openms")
fGConsNS <- consensus(fGCompNS)

fgOpenMSNSEmpty <- groupFeatures(getEmptyFeaturesNS(), "openms")
fgXCMSNSEmpty <- groupFeatures(getEmptyFeaturesNS(), "xcms")
fgCompNSOneEmpty <- comparison(openms = fgNS, xcms = fgXCMSNSEmpty, groupAlgo = "openms")
fGConsNSOneEmpty <- consensus(fgCompNSOneEmpty)
fgCompNSBothEmpty <- comparison(openms = fgOpenMSNSEmpty, xcms = fgXCMSNSEmpty, groupAlgo = "openms")
fGConsNSBothEmpty <- consensus(fgCompNSBothEmpty)

test_that("set unsupported functionality", {
    expect_length(adducts(fgNS), 0)
    expect_length(adducts(fgNSAnn), length(fgNSAnn))
    expect_length(adducts(fgNSAnn2), length(fgNSAnn2))
    expect_setequal(adducts(fgNSAnn), "[M+H]+")
    expect_equal(adducts(fgNSAnn2)[3], "[M+K]+", check.attributes = FALSE)
    # verify neutral masses
    expect_true(all(sapply(seq_len(nrow(annotations(fgNSAnn2))), function(i)
    {
        ann <- annotations(fgNSAnn2)[i]
        return(isTRUE(all.equal(patRoon:::calculateMasses(ann$neutralMass, as.adduct(ann$adduct), "mz"),
                                groupInfo(fgNSAnn2)[match(ann$group, group)]$mz)))
    })))

    expect_known_value(groupTable(fGConsNS), testFile("fg-comp-cons"))
    
    expect_lt(length(consensus(fGCompNS, relMinAbundance = 1)), length(fGConsNS))
    
    expect_length(fGConsNSOneEmpty, length(fgNS))
    expect_length(fGConsNSBothEmpty, 0)
    
    expect_equal(length(consensus(fGCompNS, uniqueFrom = 1)) +
                     length(consensus(fGCompNS, uniqueFrom = 2)) +
                     length(consensus(fGCompNS, relMinAbundance = 1)), length(fGConsNS))
    expect_equal(length(consensus(fGCompNS, uniqueFrom = 1:2, uniqueOuter = TRUE)) +
                     length(consensus(fGCompNS, relMinAbundance = 1)), length(fGConsNS))
    expect_length(consensus(fGCompNS, uniqueFrom = 1:2), length(fGConsNS))
    expect_lt(length(consensus(fGCompNS, uniqueFrom = 1:2, uniqueOuter = TRUE)), length(fGConsNS))
    expect_length(consensus(fgCompNSBothEmpty, uniqueFrom = 1), 0)
    expect_length(consensus(fgCompNSBothEmpty, uniqueFrom = 1, uniqueOuter = TRUE), 0)
    
    expect_equal(expect_plot(plotVenn(fGCompNS))$intersectionCounts, length(consensus(fGCompNS, relMinAbundance = 1)))
})
