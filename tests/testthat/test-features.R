context("test-features.R")

test_that("feature finders work", {
    expect_known_value(findFeatures(testAnaInfo, "openms", logPath = NULL), testFile("ff-openms"))
    expect_known_value(findFeatures(testAnaInfo, "xcms"), testFile("ff-xcms"))
    expect_known_value(findFeatures(testAnaInfo, "envipick"), testFile("ff-envipick"))
})
