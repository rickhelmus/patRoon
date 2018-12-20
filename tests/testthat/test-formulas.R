context("formulas")

fGroups <- getTestFGroups(getTestAnaInfo()[4:5, ])
# convert to screening results to simplify things a bit
fGroups <- groupFeaturesScreening(fGroups, screenTargets(fGroups, patRoonData::targets))

fGroupsEmpty <- getEmptyTestFGroups()
plists <- generateMSPeakLists(fGroups, "mzr")
plistsEmpty <- getEmptyPLists()
plistsEmptyMS <- filter(plists, absMSIntThr = 1E9)
plistsEmptyMSMS <- filter(plists, absMSMSIntThr = 1E9)

doSIRIUS <- !is.null(getOption("patRoon.path.SIRIUS")) && nzchar(getOption("patRoon.path.SIRIUS"))

formsGF <- generateFormulas(fGroups, "genform", plists)
formsGFEmpty <- generateFormulas(fGroupsEmpty, "genform", plistsEmpty)
formsGFEmptyPL <- generateFormulas(fGroups, "genform", plistsEmpty)
formsGFEmptyPLMS <- generateFormulas(fGroups, "genform", plistsEmptyMS)
formsGFWithMSMS <- filter(formsGF, minExplainedMSMSPeaks = 1)

if (doSIRIUS)
{
    formsSIR <- generateFormulas(fGroups, "sirius", plists, logPath = NULL, calculateFeatures = FALSE)
    formsSIREmpty <- generateFormulas(fGroupsEmpty, "sirius", plistsEmpty, logPath = NULL)
    formsSIREmptyPL <- generateFormulas(fGroups, "sirius", plistsEmpty, logPath = NULL)
    formsSIREmptyPLMS <- generateFormulas(fGroups, "sirius", plistsEmptyMS, logPath = NULL)
}

test_that("verify formula generation", {
    expect_known_value(formsGF, testFile("formulas-gf"))
    expect_length(formsGFEmpty, 0)
    expect_length(formsGFEmptyPL, 0)
    expect_length(formsGFEmptyPLMS, 0)
    expect_equal(generateFormulas(fGroups, "genform", plistsEmptyMSMS),
                 generateFormulas(fGroups, "genform", plists, MSMode = "ms"))

    expect_gt(length(generateFormulas(fGroups, "genform", plists, featThreshold = 0)),
              length(formsGF))

    skip_if_not(doSIRIUS)
    expect_known_value(formsSIR, testFile("formulas-sir"))
    expect_length(formsSIREmpty, 0)
    expect_length(formsSIREmptyPL, 0)
    expect_length(formsSIREmptyPLMS, 0)
})

test_that("verify formula show output", {
    expect_known_show(formsGF, testFile("formulas-gf", text = TRUE))
    skip_if_not(doSIRIUS)
    expect_known_show(formsSIR, testFile("formulas-sir", text = TRUE))
})

test_that("basic subsetting", {
    expect_length(formsGF["nope"], 0)
    expect_equivalent(groupNames(formsGF[1:2]), groupNames(formsGF)[1:2])
    expect_equivalent(groupNames(formsGF[groupNames(formsGF)[2:3]]), groupNames(formsGF)[2:3])
    expect_equivalent(groupNames(formsGF[c(FALSE, TRUE)]), groupNames(formsGF)[c(FALSE, TRUE)])
    expect_equal(length(formsGF[FALSE]), 0)
    expect_length(formsGFEmpty[1:5], 0)

    expect_equivalent(formsGF[[2, 15]], featureFormulas(formsGF)[[2]][[groupNames(formsGF)[15]]])
    expect_equivalent(formsGF[[analyses(formsGF)[2], groupNames(formsGF)[15]]],
                      featureFormulas(formsGF)[[2]][[groupNames(formsGF)[15]]])

    expect_equivalent(formsGF[[15]], formulaTable(formsGF)[[15]])
    expect_equivalent(formsGF[[groupNames(formsGF)[15]]], formulaTable(formsGF)[[15]])
    expect_equivalent(callDollar(formsGF, groupNames(formsGF)[4]), formsGF[[4]])
})


test_that("filtering works", {
    expect_equal(length(filter(formsGF, topMost = 1)), length(groupNames(formsGF)))
})

test_that("makeTable() works", {
    expect_equal(uniqueN(makeTable(formsGF, maxFormulas = 1),
                         by = c("group", "formula")),
                 length(groupNames(formsGF)))
    # maxFragFormulas = 1: amount of unique (MS/MS) formulas should be same as
    # amount of unique fragments
    expect_equal(uniqueN(makeTable(formsGFWithMSMS, maxFragFormulas = 1),
                         by = c("group", "formula", "byMSMS", "frag_formula")),
                 uniqueN(makeTable(formsGFWithMSMS), by = c("group", "formula")))
    expect_equal(uniqueN(makeTable(formsGF, maxFormulas = 1, maxFragFormulas = 1)),
                 length(groupNames(formsGF)))

    checkmate::expect_names(names(makeTable(formsGF, countElements = c("C", "H"))),
                            must.include = c("C", "H"))
    checkmate::expect_names(names(makeTable(formsGF, countFragElements = c("C", "H"))),
                            must.include = c("frag_C", "frag_H"))


})

if (doSIRIUS)
    fCons <- consensus(formsGF, formsSIR)

test_that("consensus works", {
    expect_length(consensus(formsGF, formsGFEmpty), length(formsGF))

    skip_if_not(doSIRIUS)
    expect_known_value(fCons, testFile("formulas-cons"))
    expect_known_show(fCons, testFile("formulas-cons", text = TRUE))
    expect_lt(length(consensus(formsGF, formsSIR, formThreshold = 1)), length(fCons))
    expect_length(consensus(formsGFEmpty,
                            generateFormulas(fGroupsEmpty, "sirius", plistsEmpty)), 0)

})

test_that("reporting works", {
    expect_error(reportCSV(fGroups, getWorkPath(), formulas = formsGF), NA)
    for (grp in groupNames(formsGF))
        checkmate::expect_file_exists(getWorkPath("formulas", sprintf("%s-%s.csv", class(fGroups), grp)))


    expect_error(reportPDF(fGroups, getWorkPath(), reportFGroups = FALSE, formulas = formsGFWithMSMS,
                           MSPeakLists = plists), NA)
    for (grp in groupNames(formsGFWithMSMS))
        checkmate::expect_file_exists(getWorkPath("formulas", sprintf("%s-%s.pdf", class(fGroups), grp)))

    expect_reportMD(makeReportMD(fGroups, reportPlots = "formulas",
                                 formulas = formsGFWithMSMS, MSPeakLists = plists))
})

plotPrec <- formsGFWithMSMS[[1]][["formula"]][1]
test_that("plotting works", {
    expect_doppel("form-spec", function() plotSpec(formsGFWithMSMS, plotPrec, groupNames(formsGFWithMSMS)[1],
                                                   MSPeakLists = plists))

    # ggplot2 versions don't really work with vdiffr at the moment :(
    # expect_doppel("spec-gg", plotSpec(formsGFWithMSMS, fTable[byMSMS == TRUE, formula][1],
    #                                                 fTable[byMSMS == TRUE, group][1], plists,
    #                                                 useGGPlot2 = TRUE))
    expect_plot(print(plotSpec(formsGFWithMSMS, plotPrec, groupNames(formsGFWithMSMS)[1], MSPeakLists = plists,
                               useGGPlot2 = TRUE)))
})
