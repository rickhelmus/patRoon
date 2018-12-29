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
formsGFWithMSMS <- filter(formsGF, minExplainedFragPeaks = 1)
formsGFOC <- generateFormulas(fGroups, "genform", plists, oc = TRUE)
formsGFMS <- generateFormulas(fGroups, "genform", plists, MSMode = "ms")

if (doSIRIUS)
{
    formsSIR <- generateFormulas(fGroups, "sirius", plists, logPath = NULL, calculateFeatures = FALSE)
    formsSIREmpty <- generateFormulas(fGroupsEmpty, "sirius", plistsEmpty, logPath = NULL)
    formsSIREmptyPL <- generateFormulas(fGroups, "sirius", plistsEmpty, logPath = NULL)
    formsSIREmptyPLMS <- generateFormulas(fGroups, "sirius", plistsEmptyMS, logPath = NULL)
}

if (doDATests())
{
    # HACK: use first standard as its compounds are not touched by MS peak lists
    fgDA <- groupFeatures(findFeatures(getDAAnaInfo()[1, ], "bruker"), "openms")
    formsDA <- generateFormulas(fgDA, "bruker")
}

test_that("verify formula generation", {
    expect_known_value(formsGF, testFile("formulas-gf"))
    expect_length(formsGFEmpty, 0)
    expect_length(formsGFEmptyPL, 0)
    expect_length(formsGFEmptyPLMS, 0)
    expect_equal(generateFormulas(fGroups, "genform", plistsEmptyMSMS), formsGFMS)

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

# extra separate block: can't have >1 skip statements...
test_that("verify DA formula generation", {
    skip_if_not(doDATests())
    expect_known_value(formsDA, testFile("formulas-DA"))
    expect_known_show(formsDA, testFile("formulas-DA", text = TRUE))
})

test_that("basic subsetting", {
    expect_length(formsGF["nope"], 0)
    expect_equivalent(groupNames(formsGF[1:2]), groupNames(formsGF)[1:2])
    expect_equivalent(groupNames(formsGF[groupNames(formsGF)[2:3]]), groupNames(formsGF)[2:3])
    expect_equivalent(groupNames(formsGF[c(FALSE, TRUE)]), groupNames(formsGF)[c(FALSE, TRUE)])
    expect_equal(length(formsGF[FALSE]), 0)
    expect_length(formsGFEmpty[1:5], 0)

    expect_equivalent(formsGF[[2, 15]], formulaTable(formsGF, TRUE)[[2]][[groupNames(formsGF)[15]]])
    expect_equivalent(formsGF[[analyses(formsGF)[2], groupNames(formsGF)[15]]],
                      formulaTable(formsGF, TRUE)[[2]][[groupNames(formsGF)[15]]])

    expect_equivalent(formsGF[[15]], formulaTable(formsGF)[[15]])
    expect_equivalent(formsGF[[groupNames(formsGF)[15]]], formulaTable(formsGF)[[15]])
    expect_equivalent(callDollar(formsGF, groupNames(formsGF)[4]), formsGF[[4]])
})

test_that("filtering works", {
    expect_true(all(as.data.table(filter(formsGF, minExplainedFragPeaks = 1))$byMSMS))
    expect_lt(length(filter(formsGF, minExplainedFragPeaks = 2)),
              length(filter(formsGF, minExplainedFragPeaks = 1)))
    expect_equal(length(filter(formsGFMS, minExplainedFragPeaks = 1)), 0)

    expect_length(filter(formsGF, elements = "C1-100"), length(formsGFOC)) # all should contain carbon
    expect_length(filter(formsGF, elements = c("Na1-100", "C1-100")), length(formsGFOC)) # no sodium, but carbon should be there
    expect_length(filter(formsGF, elements = c("H1-100", "C1-100")), length(formsGF)) # presence of both shouldn't affect results
    expect_length(filter(formsGF, elements = "Na1-100"), 0) # no sodium
    expect_length(filter(formsGF, elements = "Na0-100"), length(formsGF)) # no sodium, but optional

    # same for fragments
    expect_lte(length(filter(formsGFOC, fragElements = "C1-100")),
               length(filter(formsGFOC, minExplainedFragPeaks = 1)))
    expect_lt(length(filter(formsGFWithMSMS, fragElements = "C")), length(formsGFWithMSMS)) # fragments may contain only single carbon
    expect_length(filter(formsGFWithMSMS, fragElements = "Na1-100"), 0)
    expect_length(filter(formsGFWithMSMS, fragElements = "Na0-100"), length(formsGFWithMSMS))
    expect_length(filter(formsGFMS, fragElements = "C0-100"), 0) # no MS/MS
    
    expect_length(filter(formsGFWithMSMS, lossElements = "C0-100"), length(formsGFWithMSMS))
    expect_length(filter(formsGFWithMSMS, lossElements = "Na0-100"), length(formsGFWithMSMS))
    expect_gt(length(filter(formsGFWithMSMS, lossElements = "C1-100")), 0) # NL might be empty, at least some should contain carbon though!
    expect_length(filter(formsGFWithMSMS, lossElements = "Na1-100"), 0) # no sodium
    expect_length(filter(formsGFMS, lossElements = "C0-100"), 0) # no MS/MS
    
    expect_equal(length(filter(formsGF, topMost = 1)), length(groupNames(formsGF)))
    expect_range(length(filter(formsGF, topMost = 2)),
                 c(length(groupNames(formsGF)), length(groupNames(formsGF)) * 2))
    
    expect_equivalent(filter(formsGF, scoreLimits = list(isoScore = c(-Inf, Inf))), formsGF)
    expect_lt(length(filter(formsGF, scoreLimits = list(isoScore = c(0.1, Inf)))), length(formsGF))
    expect_length(filter(formsGF, scoreLimits = list(MSMSScore = c(0, Inf))), length(formsGFWithMSMS)) # should filter away MS only formulas
    
    expect_lt(length(filter(formsGF, OM = TRUE)), length(formsGF))
})

OMTab <- as.data.table(formsGF, OM = TRUE)
test_that("as.data.table() works", {
    expect_setequal(as.data.table(formsGF)$group, groupNames(formsGF))
    expect_setequal(as.data.table(formsGF, average = TRUE)$group, groupNames(formsGF))
    expect_setequal(as.data.table(formsGF, maxFormulas = 1)$group, groupNames(formsGF))
    expect_setequal(as.data.table(formsGF, maxFragFormulas = 1)$group, groupNames(formsGF))
    expect_setequal(as.data.table(formsGF, maxFormulas = 1, maxFragFormulas = 1)$group, groupNames(formsGF))

    expect_equal(uniqueN(as.data.table(formsGF, maxFormulas = 1),
                         by = c("group", "formula")),
                 length(groupNames(formsGF)))
    # maxFragFormulas = 1: amount of unique (MS/MS) formulas should be same as
    # amount of unique fragments
    expect_equal(uniqueN(as.data.table(formsGFWithMSMS, maxFragFormulas = 1),
                         by = c("group", "formula", "byMSMS", "frag_formula")),
                 uniqueN(as.data.table(formsGFWithMSMS), by = c("group", "formula")))
    expect_equal(uniqueN(as.data.table(formsGF, maxFormulas = 1, maxFragFormulas = 1), by = "group"),
                 length(groupNames(formsGF)))
    expect_equal(uniqueN(as.data.table(formsGF, average = TRUE), by = "group"),
                 length(groupNames(formsGF)))
    expect_equal(as.data.table(formsGFMS, maxFragFormulas = 1), as.data.table(formsGFMS))
    
    checkmate::expect_names(names(as.data.table(formsGF, countElements = c("C", "H"))),
                            must.include = c("C", "H"))
    checkmate::expect_names(names(as.data.table(formsGF, countFragElements = c("C", "H"))),
                            must.include = c("frag_C", "frag_H"))
    expect_false(any(names(as.data.table(formsGFMS, countFragElements = c("C", "H"))) %in% c("frag_C", "frag_H")))
    
    checkmate::qexpectr(OMTab[, c(unlist(strsplit("CHNOPS", "")), paste0(unlist(strsplit("HNOPS", "")), "C"),
                                  "DBE_AI", "AI")], "N+")
    checkmate::expect_character(OMTab[["classification"]], min.chars = 1, any.missing = FALSE, len = nrow(OMTab))
})

if (doSIRIUS)
    fCons <- consensus(formsGF, formsSIR)

test_that("consensus works", {
    expect_length(consensus(formsGF, formsGFEmpty), length(formsGF))

    skip_if_not(doSIRIUS)
    expect_known_value(fCons, testFile("formulas-cons"))
    expect_known_show(fCons, testFile("formulas-cons", text = TRUE))
    expect_setequal(groupNames(consensus(formsGF, formsSIR)), union(groupNames(formsGF), groupNames(formsSIR)))
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
