context("features")

ffOpenMS <- findFeatures(testAnaInfo, "openms", logPath = NULL)
ffXCMS <- findFeatures(testAnaInfo, "xcms")

if (runWithEnviPick)
    ffEP <- findFeatures(testAnaInfoMzXML, "envipick")

test_that("verify feature finder output", {
    expect_known_value(featureTable(ffOpenMS), testFile("ff-openms"))
    expect_known_value(featureTable(ffXCMS), testFile("ff-xcms"))
    skip_if_not(runWithEnviPick)
    expect_known_value(featureTable(ffEP), testFile("ff-envipick"))
})

test_that("verify show output", {
    expect_known_output(show(ffOpenMS), testFile("ff-show-openms", text = TRUE))
    expect_known_output(show(ffXCMS), testFile("ff-show-xcms", text = TRUE))
    skip_if_not(runWithEnviPick)
    expect_known_output(show(ffEP), testFile("ff-show-envipick", text = TRUE))
})
