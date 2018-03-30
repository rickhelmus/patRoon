context("test-features.R")

ffOpenMS <- findFeatures(testAnaInfo, "openms", logPath = NULL)
ffXCMS <- findFeatures(testAnaInfo, "xcms")
ffEP <- findFeatures(testAnaInfo, "envipick")

test_that("verify feature finder output", {
    expect_known_value(ffOpenMS, testFile("ff-openms"))
    expect_known_value(ffXCMS, testFile("ff-xcms"))
    expect_known_value(ffEP, testFile("ff-envipick"))
})

test_that("verify show output", {
    expect_known_output(show(ffOpenMS), testFile("ff-show-openms"))
    expect_known_output(show(ffXCMS), testFile("ff-show-xcms"))
    expect_known_output(show(ffEP), testFile("ff-show-envipick"))
})
