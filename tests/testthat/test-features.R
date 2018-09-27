context("features")

anaInfo <- getTestAnaInfo()[1:3, ]
anaInfoOne <- getTestAnaInfo()[4, ]

ffOpenMS <- findFeatures(anaInfo, "openms", logPath = NULL)
ffXCMS <- findFeatures(anaInfoOne, "xcms")

# generate mzXML files for enviPick
exDataFiles <- list.files(patRoonData::exampleDataPath(), "\\.mzML$", full.names = TRUE)
convertMSFiles(anaInfoOne, "mzML", "mzXML", getWorkPath())
epAnaInfo <- getTestAnaInfo(getWorkPath())
ffEP <- findFeatures(epAnaInfo, "envipick")
ffEmpty <- findFeatures(anaInfoOne, "openms", thr = 1E9, logPath = NULL)

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
    expect_length(ffEmpty, 0)
    expect_length(suppressWarnings(findFeatures(anaInfoOne, "xcms", snthresh = 1E9)), 0)
    expect_length(findFeatures(epAnaInfo, "envipick", minint = 1E8, maxint = 1E9), 0)
})

test_that("basic subsetting", {
    expect_length(ffOpenMS["nope"], 0)
    expect_equivalent(analyses(ffOpenMS[1:2]), anaInfo$analysis[1:2])
    expect_equivalent(analyses(ffOpenMS[anaInfo$analysis[2:3]]), anaInfo$analysis[2:3])
    expect_equivalent(analyses(ffOpenMS[c(FALSE, TRUE, FALSE)]), anaInfo$analysis[2])
    expect_equal(length(ffOpenMS[FALSE]), 0)
    expect_length(ffEmpty[1:5], 0)
    
    expect_equivalent(ffOpenMS[[2]], featureTable(ffOpenMS)[[2]])
    expect_equivalent(ffOpenMS[[analyses(ffOpenMS)[2]]], featureTable(ffOpenMS)[[2]])
    expect_equivalent(callDollar(ffOpenMS, analyses(ffOpenMS)[2]), ffOpenMS[[2]])
})
