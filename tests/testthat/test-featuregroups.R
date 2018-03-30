context("feature groups")

fList <- findFeatures(testAnaInfo, "openms", logPath = NULL)

fgOpenMS <- groupFeatures(fList, "openms")
fgXCMS <- groupFeatures(fList, "xcms")

test_that("verify feature grouping output", {
    expect_known_value(groups(fgOpenMS), testFile("fg-openms"))
    expect_known_value(groups(fgXCMS), testFile("fg-xcms"))
})

test_that("verify show output", {
    expect_known_output(show(fgOpenMS), testFile("fg-show-openms", text = TRUE))
    expect_known_output(show(fgXCMS), testFile("fg-show-xcms", text = TRUE))
})
