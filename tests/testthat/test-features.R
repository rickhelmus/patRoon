context("features")

anaInfo <- getTestAnaInfo()[4, ]

ffOpenMS <- findFeatures(anaInfo, "openms", logPath = NULL)
ffXCMS <- findFeatures(anaInfo, "xcms")

# generate mzXML files for enviPick
exDataFiles <- list.files(patRoonData::exampleDataPath(), "\\.mzML$", full.names = TRUE)
convertMSFiles(anaInfo, "mzML", "mzXML", getWorkPath())
epAnaInfo <- getTestAnaInfo(getWorkPath())
ffEP <- findFeatures(epAnaInfo, "envipick")

test_that("verify feature finder output", {
    # Don't store ID column: not reproducible
    expect_known_value(sapply(featureTable(ffOpenMS), function(fts) fts[, -"ID"], simplify = FALSE),
                       testFile("ff-openms"), tolerance = 1E-5) # increased tolerance value for win/lin deviations
    expect_known_value(featureTable(ffXCMS), testFile("ff-xcms"))
    expect_known_value(featureTable(ffEP), testFile("ff-envipick"))
})

test_that("verify show output", {
    expect_known_show(ffOpenMS, testFile("ff-show-openms", text = TRUE))
    expect_known_show(ffXCMS, testFile("ff-show-xcms", text = TRUE))
    expect_known_show(ffEP, testFile("ff-show-envipick", text = TRUE))
})

test_that("verify empty object can be generated", {
    expect_length(findFeatures(anaInfo, "openms", thr = 1E9, logPath = NULL), 0)
    expect_length(suppressWarnings(findFeatures(anaInfo, "xcms", snthresh = 1E9)), 0)
    expect_length(findFeatures(epAnaInfo, "envipick", minint = 1E8, maxint = 1E9), 0)
})