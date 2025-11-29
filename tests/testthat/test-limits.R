# SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

test_that("limits tests", {
    tmpLimitsFile <- withr::local_tempfile(fileext = ".yml")
    genLimitsFile(tmpLimitsFile)
    checkmate::expect_file_exists(tmpLimitsFile)
    expect_equal(makeFileHash(tmpLimitsFile), makeFileHash(system.file("misc", "limits.yml", package = "patRoon")))

    # test customization    
    lims <- readYAML(tmpLimitsFile)
    lims$retention$narrow <- lims$retention$narrow * 2
    tmpLimitsFileMod <- withr::local_tempfile(fileext = ".yml")
    writeYAML(lims, tmpLimitsFileMod)
    withOpt(path.limits = tmpLimitsFileMod, {
        expect_equal(defaultLim("retention", "narrow"), readYAML(tmpLimitsFile)$retention$narrow * 2)
    })
    
    # change again to test caching
    lims$retention$narrow <- lims$retention$narrow * 2
    writeYAML(lims, tmpLimitsFileMod)
    withOpt(path.limits = tmpLimitsFileMod, {
        expect_equal(defaultLim("retention", "narrow"), readYAML(tmpLimitsFile)$retention$narrow * 4) # 2*2
    })
    
    tmpLimitsFileIMSAgilent <- withr::local_tempfile(fileext = ".yml")
    genLimitsFile(tmpLimitsFileIMSAgilent, "agilent")
    withOpt(path.limits = tmpLimitsFileIMSAgilent, {
        expect_equal(defaultLim("mobility", "narrow"), readYAML(tmpLimitsFile)$mobility_agilent$narrow)
    })
})
