# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

context("reporting")

fGroups <- getTestFGroups()[1:2, 1:25]
makeReportHTML(fGroups)
makeReportHTML(fGroups, path = getWorkPath("reportNoChroms"), overrideSettings = list(features = list(chromatograms = list(large = FALSE, small = FALSE, features = FALSE))))

test_that("basic reportHTML() usage", {
    expect_lt(file.size(getWorkPath("reportNoChroms/report.html")), file.size(getWorkPath("report.html")))
    checkmate::expect_file_exists({ genReportSettingsFile(getWorkPath("test_settings.yml")); getWorkPath("test_settings.yml")})
})

plotPath <- getWorkPath("reportCleanup/report_files/plots")
mkdirp(plotPath)
cat("dummy", file = file.path(plotPath, "unused_new.svg"))
cat("dummy", file =  file.path(plotPath, "unused_old.svg"))
Sys.setFileTime(file.path(plotPath, "unused_old.svg"), as.POSIXct("1970-01-01"))
makeReportHTML(fGroups, path = getWorkPath("reportCleanup"))
test_that("removal of old plot files works", {
    checkmate::expect_file_exists(file.path(plotPath, "unused_new.svg"))
    expect_false(file.exists(file.path(plotPath, "unused_old.svg")))
})

test_that("repeated report() calls work", {
    for (i in seq_len(15))
        expect_reportHTML(makeReportHTML(fGroups))
})
