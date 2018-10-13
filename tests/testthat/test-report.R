context("reporting")

fGroups <- getTestFGroups()[, 1:25]

test_that("repeated reportMD() calls work", {
    for (i in seq_len(15))
        expect_reportMD(makeReportMD(fGroups, reportPlots = "none"))
})
