context("reporting")

fGroups <- getTestFGroups()[, 1:25]

test_that("repeated reportHTML() calls work", {
    for (i in seq_len(15))
        expect_reportHTML(makeReportHTML(fGroups, reportPlots = "none"))
})
