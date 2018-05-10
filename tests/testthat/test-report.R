context("reporting")

fGroups <- getTestFGroups()[, 1:25]

test_that("repeated reportMD() calls work", {
    for (i in seq_len(15))
        expect_file(reportMD(fGroups, getWorkPath(), reportChord = FALSE, reportFGroups = FALSE), getWorkPath("report.html"))
})
