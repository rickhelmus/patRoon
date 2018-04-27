context("features")

ffOpenMS <- findFeatures(getTestAnaInfo(), "openms", logPath = NULL)
ffXCMS <- findFeatures(getTestAnaInfo(), "xcms")

# generate mzXML files for enviPick
exDataFiles <- list.files(patRoonData::exampleDataPath(), "\\.mzML$", full.names = TRUE)
convertMSFiles(getTestAnaInfo(), "mzML", "mzXML", getWorkPath())
ffEP <- findFeatures(getTestAnaInfo(getWorkPath()), "envipick")

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
