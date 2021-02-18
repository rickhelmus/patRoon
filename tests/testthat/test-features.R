context("features")

initXCMS()

anaInfo <- getTestAnaInfo()
anaInfoOne <- getTestAnaInfo()[4, ]

ffOpenMS <- getTestFeatures(anaInfo)
ffXCMS <- findFeatures(anaInfoOne, "xcms")
ffXCMS3 <- findFeatures(anaInfoOne, "xcms3")

# generate mzXML files for enviPick
exDataFiles <- list.files(patRoonData::exampleDataPath(), "\\.mzML$", full.names = TRUE)
# convertMSFiles(file.path(anaInfoOne$path, paste0(anaInfoOne$analysis, ".mzML")),
#                outPath = getWorkPath(), to = "mzXML", algorithm = "openms",
#                logPath = NULL)
epAnaInfo <- makeMZXMLs(anaInfoOne)
ffEP <- findFeatures(epAnaInfo, "envipick")

ffEmpty <- getTestFeatures(anaInfoOne, noiseThrInt = 1E9)

if (doDATests())
{
    anaInfoDA <- getDAAnaInfo()[1, ]
    ffDA <- findFeatures(getDAAnaInfo()[1, ], "bruker")

    # NOTE: use 2nd analysis here so first can be re-used for MS peaklists/formulas...
    ffDAEmpty <- findFeatures(getDAAnaInfo()[2, ], "bruker", endRange = 0.01, doFMF = "force")
}

# Remove ID column: not reproducible
OpenMSFTable <- function(ff) sapply(featureTable(ff), function(fts) fts[, -"ID"], simplify = FALSE)

test_that("verify feature finder output", {
    expect_known_value(OpenMSFTable(ffOpenMS), testFile("ff-openms"), tolerance = 1E-5) # increased tolerance value for win/lin deviations
    expect_known_value(featureTable(ffXCMS), testFile("ff-xcms"))
    expect_known_value(featureTable(ffXCMS3), testFile("ff-xcms3"))
    expect_known_value(featureTable(ffEP), testFile("ff-envipick"))
    
    # extraOpts
    expect_equal(OpenMSFTable(ffOpenMS),
                 OpenMSFTable(getTestFeatures(anaInfo,
                                              extraOpts = list("-algorithm:common:noise_threshold_int" = 1000))))

    skip_if_not(doDATests())
    expect_known_value(featureTable(ffDA), testFile("ff-DA"))
})

test_that("verify show output", {
    expect_known_show(ffOpenMS, testFile("ff-show-openms", text = TRUE))
    expect_known_show(ffXCMS, testFile("ff-show-xcms", text = TRUE))
    expect_known_show(ffXCMS3, testFile("ff-show-xcms3", text = TRUE))
    expect_known_show(ffEP, testFile("ff-show-envipick", text = TRUE))

    skip_if_not(doDATests())
    expect_known_show(ffDA, testFile("ff-DA", text = TRUE))
})

test_that("verify empty object can be generated", {
    expect_length(ffEmpty, 0)
    expect_length(suppressWarnings(findFeatures(anaInfoOne, "xcms", snthresh = 1E9)), 0)
    expect_length(findFeatures(epAnaInfo, "envipick", minint = 1E8, maxint = 1E9), 0)

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

test_that("basic filtering", {
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
})

XCMSImpXCMS <- getXCMSSet(ffXCMS)
XCMSImpXCMS3 <- getXCMSSet(ffXCMS3, exportedData = FALSE)
XCMSImpOpenMS <- doExportXCMS(ffOpenMS, exportedData = FALSE)
XCMSImpEP <- getXCMSSet(ffEP, exportedData = FALSE)
featMZs <- function(f) lapply(featureTable(f), "[[", "mz")
test_that("XCMS conversion", {
    expect_equal(nrow(xcms::peaks(XCMSImpXCMS)), length(ffXCMS))
    expect_equal(nrow(xcms::peaks(XCMSImpXCMS3)), length(ffXCMS3))
    expect_equal(nrow(xcms::peaks(XCMSImpOpenMS)), length(getExpFeats(ffOpenMS)))
    expect_equal(nrow(xcms::peaks(XCMSImpEP)), length(ffEP))
    
    expect_known_value(xcms::peaks(XCMSImpXCMS), testFile("ff-xcms_import_xcms"))
    expect_known_value(xcms::peaks(XCMSImpXCMS3), testFile("ff-xcms_import_xcms3"))
    expect_known_value(xcms::peaks(XCMSImpOpenMS), testFile("ff-xcms_import_openms"))
    expect_known_value(xcms::peaks(XCMSImpEP), testFile("ff-xcms_import_ep"))
    
    expect_equal(featMZs(importFeatures(anaInfoOne, "xcms", XCMSImpXCMS)), featMZs(ffXCMS))
    expect_equal(featMZs(importFeatures(anaInfoOne, "xcms", XCMSImpXCMS3)), featMZs(ffXCMS3))
    expect_equal(featMZs(importFeatures(epAnaInfo, "xcms", XCMSImpEP)), featMZs(ffEP))
    
    skip_if(testWithSets())
    expect_equal(featMZs(importFeatures(anaInfo, "xcms", XCMSImpOpenMS)), featMZs(ffOpenMS))    
})

XCMS3ImpXCMS <- getXCMSnExp(ffXCMS, exportedData = FALSE)
XCMS3ImpXCMS3 <- getXCMSnExp(ffXCMS3)
XCMS3ImpOpenMS <- doExportXCMS3(ffOpenMS, exportedData = FALSE)
XCMS3ImpEP <- getXCMSnExp(ffEP, exportedData = FALSE)
test_that("XCMS3 conversion", {
    expect_equal(nrow(xcms::chromPeaks(XCMS3ImpXCMS)), length(ffXCMS))
    expect_equal(nrow(xcms::chromPeaks(XCMS3ImpXCMS3)), length(ffXCMS3))
    expect_equal(nrow(xcms::chromPeaks(XCMS3ImpOpenMS)), length(getExpFeats(ffOpenMS)))
    expect_equal(nrow(xcms::chromPeaks(XCMS3ImpEP)), length(ffEP))
    
    expect_known_value(xcms::chromPeaks(XCMS3ImpXCMS), testFile("ff-xcms3_import_xcms"))
    expect_known_value(xcms::chromPeaks(XCMS3ImpXCMS3), testFile("ff-xcms3_import_xcms3"))
    expect_known_value(xcms::chromPeaks(XCMS3ImpOpenMS), testFile("ff-xcms3_import_openms"))
    expect_known_value(xcms::chromPeaks(XCMS3ImpEP), testFile("ff-xcms3_import_ep"))
    
    expect_equal(featMZs(importFeatures(anaInfoOne, "xcms3", XCMS3ImpXCMS)), featMZs(ffXCMS))
    expect_equal(featMZs(importFeatures(anaInfoOne, "xcms3", XCMS3ImpXCMS3)), featMZs(ffXCMS3))
    expect_equal(featMZs(importFeatures(epAnaInfo, "xcms3", XCMS3ImpEP)), featMZs(ffEP))
    
    skip_if(testWithSets())
    expect_equal(featMZs(importFeatures(anaInfo, "xcms3", XCMS3ImpOpenMS)), featMZs(ffOpenMS))
})

test_that("Sets functionality", {
    skip_if_not(testWithSets())
    
    # proper (de)neutralization
    expect_equal(mean(unset(ffOpenMS, "set1")[[1]]$mz) - mean(ffOpenMS[[1]]$mz),
                 patRoon:::adductMZDelta(as.adduct("[M+H]+")))
    expect_equal(analysisInfo(unset(ffOpenMS, "set1")), getTestAnaInfoSet1())
    expect_equal(analysisInfo(ffOpenMS[, sets = "set1"])[, 1:4], getTestAnaInfoSet1())
    expect_equal(unique(ffOpenMS[[1]]$adduct), "[M+H]+")
})
