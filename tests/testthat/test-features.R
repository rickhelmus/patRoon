context("features")

anaInfo <- getTestAnaInfo()[1:3, ]
anaInfoOne <- getTestAnaInfo()[4, ]

ffOpenMS <- findFeatures(anaInfo, "openms", logPath = NULL)
ffXCMS <- findFeatures(anaInfoOne, "xcms")

# generate mzXML files for enviPick
exDataFiles <- list.files(patRoonData::exampleDataPath(), "\\.mzML$", full.names = TRUE)
# convertMSFiles(file.path(anaInfoOne$path, paste0(anaInfoOne$analysis, ".mzML")),
#                outPath = getWorkPath(), to = "mzXML", algorithm = "openms",
#                logPath = NULL)
epAnaInfo <- makeMZXMLs(anaInfoOne)
ffEP <- findFeatures(epAnaInfo, "envipick")

ffEmpty <- findFeatures(anaInfoOne, "openms", noiseThrInt = 1E9, logPath = NULL)

if (doDATests())
{
    anaInfoDA <- getDAAnaInfo()[1, ]
    ffDA <- findFeatures(getDAAnaInfo()[1, ], "bruker")

    # NOTE: use 2nd analysis here so first can be re-used for MS peaklists/formulas...
    ffDAEmpty <- findFeatures(getDAAnaInfo()[2, ], "bruker", endRange = 0.01, doFMF = "force")
}

test_that("verify feature finder output", {
    # Don't store ID column: not reproducible
    expect_known_value(sapply(featureTable(ffOpenMS), function(fts) fts[, -"ID"], simplify = FALSE),
                       testFile("ff-openms"), tolerance = 1E-5) # increased tolerance value for win/lin deviations
    expect_known_value(featureTable(ffXCMS), testFile("ff-xcms"))
    expect_known_value(featureTable(ffEP), testFile("ff-envipick"))

    skip_if_not(doDATests())
    expect_known_value(featureTable(ffDA), testFile("ff-DA"))
})

test_that("verify show output", {
    expect_known_show(ffOpenMS, testFile("ff-show-openms", text = TRUE))
    expect_known_show(ffXCMS, testFile("ff-show-xcms", text = TRUE))
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
    expect_equivalent(analyses(ffOpenMS[c(FALSE, TRUE, FALSE)]), anaInfo$analysis[2])
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
