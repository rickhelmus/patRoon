context("features")

initXCMS()

anaInfo <- getTestAnaInfo()
anaInfoOne <- getTestAnaInfo()[4, ]

ffOpenMS <- getTestFeatures(anaInfo)
# ffXCMS <- findFeatures(anaInfoOne, "xcms", noise = 3E4)
ffXCMS3 <- findFeatures(anaInfoOne, "xcms3", xcms::CentWaveParam(noise = 3E4))
# UNDONE: ignore warnings about clusters...
ffKPIC2 <- withCallingHandlers(findFeatures(anaInfoOne, "kpic2", level = 1E5),
                               warning = function(w) if (grepl("number of clusters", w, fixed = TRUE)) invokeRestart("muffleWarning"))
ffSIRIUS <- findFeatures(anaInfoOne, "sirius")

# generate mzXML files for enviPick
exDataFiles <- list.files(patRoonData::exampleDataPath(), "\\.mzML$", full.names = TRUE)
epAnaInfo <- makeMZXMLs(anaInfoOne)
ffEP <- findFeatures(epAnaInfo, "envipick", minpeak = 25)

ffOpenMSQ <- calculatePeakQualities(ffOpenMS)

ffEmpty <- getTestFeatures(anaInfoOne, noiseThrInt = 1E9)
ffEmptyQ <- calculatePeakQualities(ffEmpty)


if (doDATests())
{
    ffDA <- findFeatures(getDAAnaInfo("std1")[1, ], "bruker")

    # NOTE: use 2nd analysis here so first can be re-used for MS peaklists/formulas...
    ffDAEmpty <- findFeatures(getDAAnaInfo("std2")[1, ], "bruker", endRange = 0.01, doFMF = "force")
}

# Remove ID column: not reproducible
OpenMSFTable <- function(ff) sapply(featureTable(ff), function(fts) fts[, -"ID"], simplify = FALSE)

test_that("verify feature finder output", {
    expect_known_value(OpenMSFTable(ffOpenMS), testFile("ff-openms"), tolerance = 1E-5) # increased tolerance value for win/lin deviations
    # expect_known_value(featureTable(ffXCMS), testFile("ff-xcms"))
    expect_known_value(featureTable(ffXCMS3), testFile("ff-xcms3"))
    expect_known_value(featureTable(ffEP), testFile("ff-envipick"))
    expect_known_value(featureTable(ffKPIC2), testFile("ff-kpic2"))
    expect_known_value(featureTable(ffSIRIUS), testFile("ff-sirius"))
    
    expect_known_value(OpenMSFTable(ffOpenMSQ), testFile("ff-openms-qual"))
    
    # extraOpts
    expect_equal(OpenMSFTable(ffOpenMS),
                 OpenMSFTable(getTestFeatures(anaInfo,
                                              extraOpts = list("-algorithm:common:noise_threshold_int" = 30000))))

    skip_if_not(doDATests())
    expect_known_value(featureTable(ffDA), testFile("ff-DA"))
})

test_that("verify show output", {
    expect_known_show(ffOpenMS, testFile("ff-show-openms", text = TRUE))
    # expect_known_show(ffXCMS, testFile("ff-show-xcms", text = TRUE))
    expect_known_show(ffXCMS3, testFile("ff-show-xcms3", text = TRUE))
    expect_known_show(ffEP, testFile("ff-show-envipick", text = TRUE))
    expect_known_show(ffKPIC2, testFile("ff-show-kpic2", text = TRUE))
    expect_known_show(ffSIRIUS, testFile("ff-show-sirius", text = TRUE))

    expect_known_show(ffOpenMSQ, testFile("ff-show-openms-qual", text = TRUE))
    
    skip_if_not(doDATests())
    expect_known_show(ffDA, testFile("ff-DA", text = TRUE))
})

test_that("verify empty object can be generated", {
    expect_length(ffEmpty, 0)
    #expect_length(suppressWarnings(findFeatures(anaInfoOne, "xcms", noise = 1E9)), 0)
    expect_length(findFeatures(epAnaInfo, "envipick", minint = 1E8, maxint = 1E9, minpeak = 100), 0) # add minpeak to speed-up
    expect_length(ffEmptyQ, 0)

    skip_if_not(doDATests())
    expect_length(ffDAEmpty, 0)
})

test_that("basic subsetting", {
    expect_length(ffOpenMS["nope"], 0)
    expect_equivalent(analyses(ffOpenMS[1:2]), anaInfo$analysis[1:2])
    expect_equivalent(analyses(ffOpenMS[anaInfo$analysis[2:3]]), anaInfo$analysis[2:3])
    expect_equivalent(analyses(ffOpenMS[c(TRUE, FALSE)]), anaInfo$analysis[c(TRUE, FALSE)])
    expect_equal(length(ffOpenMS[FALSE]), 0)
    expect_length(ffEmpty[1:5], 0)

    expect_equivalent(ffOpenMS[[2]], featureTable(ffOpenMS)[[2]])
    expect_equivalent(ffOpenMS[[analyses(ffOpenMS)[2]]], featureTable(ffOpenMS)[[2]])
    expect_equivalent(callDollar(ffOpenMS, analyses(ffOpenMS)[2]), ffOpenMS[[2]])
})

qr <- list(ZigZag = c(0.2, 0.6), TPASRScore = c(0.5, 0.9))
ffOpenMSQF <- filter(ffOpenMSQ, qualityRange = qr)
ffOpenMSQFTab <- as.data.table(ffOpenMSQF)
ffOpenMSQFN <- filter(ffOpenMSQ, qualityRange = qr, negate = TRUE)
ffOpenMSQFNTab <- as.data.table(ffOpenMSQFN)

test_that("delete and filter", {
    checkmate::expect_names(analyses(delete(ffOpenMS, i = 1)), disjunct.from = analyses(ffOpenMS)[1])
    checkmate::expect_names(analyses(delete(ffOpenMS, i = analyses(ffOpenMS)[1])), disjunct.from = analyses(ffOpenMS)[1])
    expect_length(delete(ffOpenMS, i = analyses(ffOpenMS)), 0)
    expect_false(delete(ffOpenMS, j = 1)[[1]]$ID[1] == ffOpenMS[[1]]$ID[1])
    expect_length(delete(ffOpenMS, j = 3:4), length(ffOpenMS) - (length(analyses(ffOpenMS)) * 2))
    expect_false(delete(ffOpenMS, j = function(...) 1)[[1]]$ID[1] == ffOpenMS[[1]]$ID[1])
    expect_length(delete(ffOpenMS, j = function(...) 3:4), length(ffOpenMS) - (length(analyses(ffOpenMS)) * 2))
    expect_length(delete(ffOpenMS, j = function(...) TRUE), 0)
    expect_equal(delete(ffOpenMS, i = character()), ffOpenMS)
    expect_equal(delete(ffOpenMS, j = integer()), ffOpenMS)
    expect_length(delete(ffOpenMS), 0)
    expect_length(delete(ffEmpty), 0)
    
    expect_gte(min(filter(ffOpenMS, absMinIntensity = 500)[[1]]$intensity), 500)
    expect_gte(min(filter(ffOpenMS, relMinIntensity = 0.2)[[1]]$intensity), 0.2 * max(ffOpenMS[[1]]$intensity))

    expect_range(filter(ffOpenMS, retentionRange = c(120, 300))[[1]]$ret, c(120, 300))
    expect_equivalent(filter(ffOpenMS, retentionRange = c(0, Inf)), ffOpenMS)
    expect_range(filter(ffOpenMS, mzRange = c(200, 300))[[1]]$mz, c(200, 300))
    expect_equivalent(filter(ffOpenMS, mzRange = c(0, Inf)), ffOpenMS)
    expect_range(filter(ffOpenMS, mzDefectRange = c(0.1, 0.2))[[1]]$mz %% 1, c(0.1, 0.2))
    expect_equivalent(filter(ffOpenMS, mzDefectRange = c(0, 1)), ffOpenMS)
    expect_lt(length(filter(ffOpenMS, chromWidthRange = c(0, 30))), length(ffOpenMS))
    expect_equivalent(filter(ffOpenMS, chromWidthRange = c(0, Inf)), ffOpenMS)

    expect_range(ffOpenMSQFTab[[names(qr)[1]]], qr[[1]])
    expect_range(ffOpenMSQFTab[[names(qr)[2]]], qr[[2]])
    expect_true(all(ffOpenMSQFNTab[[names(qr)[1]]] < qr[[c(1, 1)]] | ffOpenMSQFNTab[[names(qr)[1]]] > qr[[c(1, 2)]]))
    expect_true(all(ffOpenMSQFNTab[[names(qr)[2]]] < qr[[c(2, 1)]] | ffOpenMSQFNTab[[names(qr)[2]]] > qr[[c(2, 2)]]))
    
    expect_known_output(filter(ffOpenMS, absMinIntensity = 500, retentionRange = c(120, Inf),
                               mzRange = c(100, 400)),
                        testFile("ff-combi", text = TRUE))
    expect_known_output(filter(ffOpenMS, absMinIntensity = 500, retentionRange = c(120, Inf),
                               mzRange = c(100, 400), negate = TRUE),
                        testFile("ff-combi-neg", text = TRUE))
    expect_length(filter(ffEmpty, absMinIntensity = 500, retentionRange = c(120, Inf),
                         mzRange = c(100, 400)), 0)
})

test_that("basic usage", {
    expect_equal(nrow(as.data.table(ffOpenMS)), length(ffOpenMS))
    checkmate::expect_names(names(as.data.table(ffOpenMSQ)),
                                  must.include = c(featureQualityNames(group = FALSE),
                                                   featureQualityNames(group = FALSE, scores = TRUE)))
})

# XCMSImpXCMS <- getXCMSSet(ffXCMS)
XCMSImpXCMS3 <- getXCMSSet(ffXCMS3, loadRawData = FALSE)
XCMSImpOpenMS <- doExportXCMS(ffOpenMS, loadRawData = FALSE)
XCMSImpEP <- getXCMSSet(ffEP, loadRawData = FALSE)
XCMSImpKPIC2 <- getXCMSSet(ffKPIC2, loadRawData = FALSE)
XCMSImpSIRIUS <- getXCMSSet(ffSIRIUS, loadRawData = FALSE)
featMZs <- function(f) lapply(featureTable(f), "[[", "mz")
test_that("XCMS conversion", {
    # expect_equal(nrow(xcms::peaks(XCMSImpXCMS)), length(ffXCMS))
    expect_equal(nrow(xcms::peaks(XCMSImpXCMS3)), length(ffXCMS3))
    expect_equal(nrow(xcms::peaks(XCMSImpOpenMS)), length(getExpFeats(ffOpenMS)))
    expect_equal(nrow(xcms::peaks(XCMSImpEP)), length(ffEP))
    expect_equal(nrow(xcms::peaks(XCMSImpKPIC2)), length(ffKPIC2))
    expect_equal(nrow(xcms::peaks(XCMSImpSIRIUS)), length(ffSIRIUS))
    
    # expect_known_value(xcms::peaks(XCMSImpXCMS), testFile("ff-xcms_import_xcms"))
    expect_known_value(xcms::peaks(XCMSImpXCMS3), testFile("ff-xcms_import_xcms3"))
    expect_known_value(xcms::peaks(XCMSImpOpenMS), testFile("ff-xcms_import_openms"))
    expect_known_value(xcms::peaks(XCMSImpEP), testFile("ff-xcms_import_ep"))
    expect_known_value(xcms::peaks(XCMSImpKPIC2), testFile("ff-xcms_import_kpic2"))
    expect_known_value(xcms::peaks(XCMSImpSIRIUS), testFile("ff-xcms_import_sirius"))
    
    # expect_equal(featMZs(importFeatures(anaInfoOne, "xcms", XCMSImpXCMS)), featMZs(ffXCMS))
    expect_equal(featMZs(importFeatures(anaInfoOne, "xcms", XCMSImpXCMS3)), featMZs(ffXCMS3))
    expect_equal(featMZs(importFeatures(epAnaInfo, "xcms", XCMSImpEP)), featMZs(ffEP))
    expect_equal(featMZs(importFeatures(anaInfoOne, "xcms", XCMSImpKPIC2)), featMZs(ffKPIC2))
    expect_equal(featMZs(importFeatures(anaInfoOne, "xcms", XCMSImpSIRIUS)), featMZs(ffSIRIUS))
    
    skip_if(testWithSets())
    expect_equal(featMZs(importFeatures(anaInfo, "xcms", XCMSImpOpenMS)), featMZs(ffOpenMS))    
})

# XCMS3ImpXCMS <- getXCMSnExp(ffXCMS, loadRawData = FALSE)
XCMS3ImpXCMS3 <- getXCMSnExp(ffXCMS3)
XCMS3ImpOpenMS <- doExportXCMS3(ffOpenMS, loadRawData = FALSE)
XCMS3ImpEP <- getXCMSnExp(ffEP, loadRawData = FALSE)
XCMS3ImpKPIC2 <- getXCMSnExp(ffKPIC2, loadRawData = FALSE)
XCMS3ImpSIRIUS <- getXCMSnExp(ffSIRIUS, loadRawData = FALSE)
test_that("XCMS3 conversion", {
    # expect_equal(nrow(xcms::chromPeaks(XCMS3ImpXCMS)), length(ffXCMS))
    expect_equal(nrow(xcms::chromPeaks(XCMS3ImpXCMS3)), length(ffXCMS3))
    expect_equal(nrow(xcms::chromPeaks(XCMS3ImpOpenMS)), length(getExpFeats(ffOpenMS)))
    expect_equal(nrow(xcms::chromPeaks(XCMS3ImpEP)), length(ffEP))
    expect_equal(nrow(xcms::chromPeaks(XCMS3ImpKPIC2)), length(ffKPIC2))
    expect_equal(nrow(xcms::chromPeaks(XCMS3ImpSIRIUS)), length(ffSIRIUS))
    
    # expect_known_value(xcms::chromPeaks(XCMS3ImpXCMS), testFile("ff-xcms3_import_xcms"))
    expect_known_value(xcms::chromPeaks(XCMS3ImpXCMS3), testFile("ff-xcms3_import_xcms3"))
    expect_known_value(xcms::chromPeaks(XCMS3ImpOpenMS), testFile("ff-xcms3_import_openms"))
    expect_known_value(xcms::chromPeaks(XCMS3ImpEP), testFile("ff-xcms3_import_ep"))
    expect_known_value(xcms::chromPeaks(XCMS3ImpKPIC2), testFile("ff-xcms3_import_kpic2"))
    expect_known_value(xcms::chromPeaks(XCMS3ImpSIRIUS), testFile("ff-xcms3_import_sirius"))
    
    # expect_equal(featMZs(importFeatures(anaInfoOne, "xcms3", XCMS3ImpXCMS)), featMZs(ffXCMS))
    expect_equal(featMZs(importFeatures(anaInfoOne, "xcms3", XCMS3ImpXCMS3)), featMZs(ffXCMS3))
    expect_equal(featMZs(importFeatures(epAnaInfo, "xcms3", XCMS3ImpEP)), featMZs(ffEP))
    expect_equal(featMZs(importFeatures(anaInfoOne, "xcms3", XCMS3ImpKPIC2)), featMZs(ffKPIC2))
    expect_equal(featMZs(importFeatures(anaInfoOne, "xcms3", XCMS3ImpSIRIUS)), featMZs(ffSIRIUS))
    
    skip_if(testWithSets())
    expect_equal(featMZs(importFeatures(anaInfo, "xcms3", XCMS3ImpOpenMS)), featMZs(ffOpenMS))
})

test_that("Sets functionality", {
    skip_if_not(testWithSets())
    
    # proper (de)neutralization
    expect_equal(patRoon:::calculateMasses(unset(ffOpenMS, "positive")[[1]]$mz, as.adduct("[M+H]+"), "neutral"),
                 ffOpenMS[[1]]$mz)
    expect_equal(analysisInfo(unset(ffOpenMS, "positive"), TRUE), getTestAnaInfoPos())
    expect_equal(analysisInfo(ffOpenMS[, sets = "positive"], TRUE)[, 1:4], getTestAnaInfoPos())
    expect_equal(unique(ffOpenMS[[1]]$adduct), "[M+H]+")
    expect_equal(sets(filter(ffOpenMS, sets = "positive", negate = TRUE)), "negative")
    expect_length(ffOpenMS[, sets = character()], 0)
    expect_length(makeSet(ffXCMS3, ffXCMS3[FALSE], adducts = "[M+H]+"), length(ffXCMS3))
    expect_length(makeSet(ffXCMS3[FALSE], ffXCMS3[FALSE], adducts = "[M+H]+"), 0)
})
