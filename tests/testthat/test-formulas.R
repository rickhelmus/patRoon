context("formulas")

fGroups <- getFormFGroups()
fGroupsEmpty <- getEmptyTestFGroups()
plists <- generateMSPeakLists(fGroups, "mzr")
plistsEmpty <- plists[FALSE, reAverage = TRUE]
plistsEmptyMS <- removeMSPlists(plists, "MS")
plistsEmptyMSMS <- removeMSPlists(plists, "MSMS")

doSIRIUS <- TRUE # !is.null(getOption("patRoon.path.SIRIUS")) && nzchar(getOption("patRoon.path.SIRIUS"))

formsGF <- doGenForms(fGroups, plists, "genform")
formsGFEmpty <- doGenForms(fGroupsEmpty, plistsEmpty, "genform")
formsGFEmptyPL <- doGenForms(fGroups, plistsEmpty, "genform")
formsGFEmptyPLMS <- doGenForms(fGroups, plistsEmptyMS, "genform")
formsGFWithMSMS <- filter(formsGF, minExplainedPeaks = 1)
formsGFOC <- doGenForms(fGroups, plists, "genform", oc = TRUE)
formsGFMS <- doGenForms(fGroups, plists, "genform", MSMode = "ms")

if (doSIRIUS)
{
    formsSIR <- doGenForms(fGroups, plists, "sirius", calculateFeatures = FALSE)
    formsSIREmpty <- doGenForms(fGroupsEmpty, plistsEmpty, "sirius")
    formsSIREmptyPL <- doGenForms(fGroups, plistsEmpty, "sirius")
    formsSIREmptyPLMS <- doGenForms(fGroups, plistsEmptyMS, "sirius")
    
    if (FALSE)
        updateSIRIUSFormFPsProj(fGroups, plists)
    
    formsSIRFPs <- doGenFormsSIRFPs(fGroups, plists)
}

if (doDATests())
{
    # HACK: use first standard as its compounds are not touched by MS peak lists
    fgDA <- getTestFGroupsDA(getDAAnaInfo("std1"))
    plistsDA <- generateMSPeakLists(fgDA, "brukerfmf")
    formsDA <- doGenForms(fgDA, plistsDA, "bruker")
}

test_that("verify formula generation", {
    expect_known_value(formsGF, testFile("formulas-gf"))
    expect_length(formsGFEmpty, 0)
    expect_length(formsGFEmptyPL, 0)
    expect_length(formsGFEmptyPLMS, 0)
    expect_equal(doGenForms(fGroups, plistsEmptyMSMS, "genform"), formsGFMS)

    expect_gt(length(doGenForms(fGroups, plists, "genform", featThresholdAnn = 0)),
              length(formsGF))

    skip_if_not(doSIRIUS)
    expect_known_value(formsSIR, testFile("formulas-sir"))
    expect_length(formsSIREmpty, 0)
    expect_length(formsSIREmptyPL, 0)
    expect_length(formsSIREmptyPLMS, 0)
    expect_known_value(formsSIRFPs, testFile("formulas-sir-fps"))
})

test_that("verify formula show output", {
    expect_known_show(formsGF, testFile("formulas-gf", text = TRUE))
    skip_if_not(doSIRIUS)
    expect_known_show(formsSIR, testFile("formulas-sir", text = TRUE))
    expect_known_show(formsSIRFPs, testFile("formulas-sir-fps", text = TRUE))
})

# extra separate block: can't have >1 skip statements...
test_that("verify DA formula generation", {
    skip_if_not(doDATests())
    expect_known_value(formsDA, testFile("formulas-DA"))
    expect_known_show(formsDA, testFile("formulas-DA", text = TRUE))
})

test_that("verify fingerprints", {
    skip_if(!doSIRIUS || testWithSets())
    expect_gt(length(formsSIRFPs), 0)
    testSIRFPSubset(formsSIRFPs)
})

test_that("basic subsetting", {
    expect_length(formsGF["nope"], 0)
    expect_equivalent(groupNames(formsGF[1:2]), groupNames(formsGF)[1:2])
    expect_equivalent(groupNames(formsGF[groupNames(formsGF)[2:3]]), groupNames(formsGF)[2:3])
    expect_equivalent(groupNames(formsGF[c(FALSE, TRUE)]), groupNames(formsGF)[c(FALSE, TRUE)])
    expect_equal(length(formsGF[FALSE]), 0)
    expect_length(formsGFEmpty[1:5], 0)

    expect_equivalent(formsGF[[2, 5]], annotations(formsGF, TRUE)[[2]][[groupNames(formsGF)[5]]])
    expect_equivalent(formsGF[[analyses(formsGF)[2], groupNames(formsGF)[5]]],
                      annotations(formsGF, TRUE)[[2]][[groupNames(formsGF)[5]]])

    expect_equivalent(formsGF[[5]], annotations(formsGF)[[5]])
    expect_equivalent(formsGF[[groupNames(formsGF)[5]]], annotations(formsGF)[[5]])
    expect_equivalent(callDollar(formsGF, groupNames(formsGF)[4]), formsGF[[4]])
})

if (!testWithSets())
{
    formsGFMST1 <- filter(formsGFMS, topMost = 1)
    formsGFMST1N <- filter(formsGFMS, topMost = 1, negate = TRUE)
}

test_that("delete and filter", {
    checkmate::expect_names(groupNames(delete(formsGF, i = 1)), disjunct.from = groupNames(formsGF)[1])
    checkmate::expect_names(groupNames(delete(formsGF, i = groupNames(formsGF)[1])), disjunct.from = groupNames(formsGF)[1])
    expect_length(delete(formsGF, i = groupNames(formsGF)), 0)
    expect_false(delete(formsGF, j = 1)[[1]]$UID[1] == formsGF[[1]]$UID[1])
    expect_length(delete(formsGF, j = 1), length(formsGF) - length(groupNames(formsGF)))
    expect_false(delete(formsGF, j = function(...) 1)[[1]]$UID[1] == formsGF[[1]]$UID[1])
    expect_length(delete(formsGF, j = function(...) 1), length(formsGF) - length(groupNames(formsGF)))
    expect_length(delete(formsGF, j = function(...) TRUE), 0)
    expect_equal(delete(formsGF, i = character()), formsGF)
    expect_equal(delete(formsGF, j = integer()), formsGF)
    expect_length(delete(formsGF), 0)
    expect_length(delete(formsGFEmpty), 0)
    
    expect_true(all(as.data.table(filter(formsGF, minExplainedPeaks = 1))$explainedPeaks >= 1))
    expect_true(all(as.data.table(filter(formsGF, minExplainedPeaks = 1, negate = TRUE))$explainedPeaks == 0))
    expect_lt(length(filter(formsGF, minExplainedPeaks = 2)),
              length(filter(formsGF, minExplainedPeaks = 1)))
    expect_gt(length(filter(formsGF, minExplainedPeaks = 2, negate = TRUE)),
              length(filter(formsGF, minExplainedPeaks = 1, negate = TRUE)))
    expect_length(filter(formsGFMS, minExplainedPeaks = 1), 0)
    expect_length(filter(formsGFMS, minExplainedPeaks = 1, negate = TRUE), length(formsGFMS))

    expect_true(all(grepl("^C", as.data.table(filter(formsGF, elements = "C1-100"))$neutral_formula))) # all should contain carbon
    expect_length(filter(formsGFOC, elements = "C1-100", negate = TRUE), 0)
    expect_true(all(grepl("^C", as.data.table(filter(formsGF, elements = c("Na1-100", "C1-100")))$neutral_formula))) # no sodium, but carbon should be there
    expect_length(filter(formsGFOC, elements = c("Na1-100", "C1-100"), negate = TRUE), length(formsGFOC)) # no sodium
    expect_length(filter(formsGF, elements = c("H0-100", "C1-100")), length(formsGF)) # presence of both shouldn't affect results
    expect_length(filter(formsGF, elements = c("H90-100", "C0-100"), negate = TRUE), length(formsGF)) # all formulae shouldn't have 90-100 Hs
    expect_length(filter(formsGF, elements = "Na1-100"), 0) # no sodium
    expect_length(filter(formsGF, elements = "Na1-100", negate = TRUE), length(formsGF))
    expect_length(filter(formsGF, elements = "Na0-100"), length(formsGF)) # no sodium, but optional
    expect_length(filter(formsGF, elements = "Na0-100", negate = TRUE), 0)

    # same for fragments
    expect_lte(length(filter(formsGFOC, fragElements = "C1-100")),
               length(filter(formsGFOC, minExplainedPeaks = 1)))
    expect_lte(length(filter(formsGFOC, fragElements = "C1-100", negate = TRUE)),
               length(filter(formsGFOC, minExplainedPeaks = 1)))
    expect_lt(length(filter(formsGFWithMSMS, fragElements = "C")),
              length(formsGFWithMSMS)) # >=1 fragments may contain only single carbon
    expect_lt(length(filter(formsGFWithMSMS, fragElements = "C", negate = TRUE)),
                     length(formsGFWithMSMS)) # >=1 fragments may not contain only single carbon
    expect_length(filter(formsGFWithMSMS, fragElements = "Na1-100"), 0)
    expect_length(filter(formsGFWithMSMS, fragElements = "Na1-100", negate = TRUE), length(formsGFWithMSMS))
    expect_length(filter(formsGFWithMSMS, fragElements = "Na0-100"), length(formsGFWithMSMS))
    expect_length(filter(formsGFWithMSMS, fragElements = "Na0-100", negate = TRUE), 0)
    expect_length(filter(formsGFMS, fragElements = "C0-100"), 0) # no MS/MS
    expect_length(filter(formsGFMS, fragElements = "C0-100", negate = TRUE), 0)

    expect_length(filter(formsGFWithMSMS, lossElements = "C0-100"), length(formsGFWithMSMS))
    expect_length(filter(formsGFWithMSMS, lossElements = "C0-100", negate = TRUE), 0)
    expect_length(filter(formsGFWithMSMS, lossElements = "Na0-100"), length(formsGFWithMSMS))
    expect_length(filter(formsGFWithMSMS, lossElements = "Na0-100", negate = TRUE), 0)
    expect_gt(length(filter(formsGFWithMSMS, lossElements = "C1-100")), 0) # NL might be empty, at least some should contain carbon though!
    expect_gt(length(filter(formsGFWithMSMS, lossElements = "C1-100", negate = TRUE)), 0)
    expect_length(filter(formsGFWithMSMS, lossElements = "Na1-100"), 0) # no sodium
    expect_length(filter(formsGFWithMSMS, lossElements = "Na1-100", negate = TRUE), length(formsGFWithMSMS))
    expect_length(filter(formsGFMS, lossElements = "C0-100"), 0) # no MS/MS
    expect_length(filter(formsGFMS, lossElements = "C0-100", negate = TRUE), 0)

    expect_equal(length(groupNames(filter(formsGF, topMost = 1))), length(groupNames(formsGF)))
    expect_gte(length(filter(formsGF, topMost = 1, negate = TRUE)), length(groupNames(formsGF)))
    expect_range(length(filter(formsGF, topMost = 2)),
                 c(length(groupNames(formsGF)), length(groupNames(formsGF)) * 2))
    expect_range(length(filter(formsGF, topMost = 2, negate = TRUE)),
                 c(length(groupNames(formsGF)), length(groupNames(formsGF)) * 2))

    expect_equivalent(filter(formsGF, scoreLimits = list(isoScore = c(-Inf, Inf))), formsGF)
    expect_length(filter(formsGF, scoreLimits = list(isoScore = c(-Inf, Inf)), negate = TRUE), 0)
    expect_lt(length(filter(formsGF, scoreLimits = list(isoScore = c(0.7, Inf)))), length(formsGF))
    expect_lt(length(filter(formsGF, scoreLimits = list(isoScore = c(0.7, Inf)), negate = TRUE)), length(formsGF))

    expect_lt(length(filter(formsGF, OM = TRUE)), length(formsGF))
    expect_lt(length(filter(formsGF, OM = TRUE, negate = TRUE)), length(formsGF))
    expect_length(formsGF, sum(length(filter(formsGF, OM = TRUE)),
                               length(filter(formsGF, OM = TRUE, negate = TRUE))))
    
    skip_if(testWithSets())
    
    # in case of ties between pos/neg the isoScore is sometimes not the highest --> skip test with sets for now
    expect_true(all(sapply(annotations(formsGFMST1[groupNames(formsGFMST1N)]), function(a) max(a$isoScore)) >=
                        sapply(annotations(formsGFMST1N), function(a) max(a$isoScore))))
})

test_that("as.data.table() works", {
    testFeatAnnADT(formsGF)

    normScName <- if (testWithSets()) "isoScore-positive" else "isoScore"
    expect_range(na.omit(as.data.table(formsGF, normalizeScores = "max")[[normScName]]), c(0, 1))
    expect_setequal(as.data.table(formsGF, average = TRUE)$group, groupNames(formsGF))
    expect_equal(uniqueN(as.data.table(formsGF, average = TRUE), by = "group"),
                 length(groupNames(formsGF)))
    expect_true(all(is.na(unlist(as.data.table(formsGFMS,
                                               countFragElements = c("C", "H"))[, c("frag_C", "frag_H"), with = FALSE]))))
})

if (doSIRIUS)
    fCons <- doFormCons(formsGF, formsSIR, MSPeakLists = plists)

test_that("consensus works", {
    expect_length(doFormCons(formsGF, formsGFEmpty, MSPeakLists = plists), length(formsGF))

    skip_if_not(doSIRIUS)
    expect_known_value(fCons, testFile("formulas-cons"))
    expect_known_show(fCons, testFile("formulas-cons", text = TRUE))
    expect_setequal(groupNames(doFormCons(formsGF, formsSIR, MSPeakLists = plists)),
                    union(groupNames(formsGF), groupNames(formsSIR)))
    expect_lt(length(doFormCons(formsGF, formsSIR, MSPeakLists = plists, relMinAbundance = 1)), length(fCons))
    expect_length(doFormCons(formsGFEmpty, formsSIREmpty, MSPeakLists = plists), 0)

    expect_equal(sum(lengths(list(doFormCons(formsGF, formsSIR, MSPeakLists = plists, uniqueFrom = 1),
                                  doFormCons(formsGF, formsSIR, MSPeakLists = plists, uniqueFrom = 2),
                                  doFormCons(formsGF, formsSIR, MSPeakLists = plists, relMinAbundance = 1)))),
                 length(fCons))
    expect_equal(sum(lengths(list(doFormCons(formsGF, formsSIR, MSPeakLists = plists, uniqueFrom = 1:2, uniqueOuter = TRUE),
                                  doFormCons(formsGF, formsSIR, MSPeakLists = plists, relMinAbundance = 1)))),
                 length(fCons))
    expect_length(doFormCons(formsGF, formsSIR, MSPeakLists = plists, uniqueFrom = 1:2), length(fCons))
    expect_lt(length(doFormCons(formsGF, formsSIR, MSPeakLists = plists, uniqueFrom = 1:2, uniqueOuter = TRUE)),
              length(fCons))
    expect_length(doFormCons(formsGFEmpty, formsSIREmpty, MSPeakLists = plists, uniqueFrom = 1), 0)
    expect_length(doFormCons(formsGFEmpty, formsSIREmpty, MSPeakLists = plists, uniqueFrom = 1, uniqueOuter = TRUE), 0)
})

anPLGroup <- screenInfo(fGroups)[name == "1H-benzotriazole"]$group
anPL <- annotatedPeakList(formsGF, index = 1, groupName = anPLGroup, MSPeakLists = plists)
anPLOnly <- annotatedPeakList(formsGF, index = 1, groupName = anPLGroup, MSPeakLists = plists, onlyAnnotated = TRUE)

if (doSIRIUS)
    anPLCons <- annotatedPeakList(fCons, index = 1, groupName = anPLGroup, MSPeakLists = plists, onlyAnnotated = TRUE)

test_that("annotation works", {
    skip_if_not(doSIRIUS)

    expect_lt(nrow(anPLOnly), nrow(anPL))
    expect_true(any(is.na(anPL$ion_formula)))
    expect_false(any(is.na(anPLOnly$ion_formula)))
    expect_true(all(formsGF[[anPLGroup]]$fragInfo[[1]]$ion_formula %in% anPLOnly$ion_formula))
    
    skip_if(!doSIRIUS)
    expect_true(any(grepl("genform", anPLCons$mergedBy)))
    expect_true(any(grepl("sirius", anPLCons$mergedBy)))
})

reportGroups <- groupNames(formsGFWithMSMS)[1:5]
test_that("reporting works", {
    expect_error(reportCSV(fGroups, getWorkPath(), formulas = formsGF), NA)
    for (grp in groupNames(formsGF))
        checkmate::expect_file_exists(getWorkPath("formulas", sprintf("%s-%s.csv", class(fGroups), grp)))

    expect_error(reportPDF(fGroups[, reportGroups], getWorkPath(), reportFGroups = FALSE, formulas = formsGFWithMSMS,
                           MSPeakLists = plists), NA)
    for (grp in reportGroups)
        checkmate::expect_file_exists(getWorkPath("formulas", sprintf("%s-%s.pdf", class(fGroups), grp)))

    expect_reportHTML(makeReportHTML(fGroups[, reportGroups], formulas = formsGFWithMSMS, MSPeakLists = plists))
})

test_that("reporting empty objects works", {
    expect_error(reportCSV(fGroups, getWorkPath(), formulas = formsGFEmpty), NA)
    expect_error(reportCSV(fGroupsEmpty, getWorkPath(), formulas = formsGF), NA)
    expect_error(reportPDF(fGroups, getWorkPath(), reportFGroups = FALSE, formulas = formsGFEmpty,
                           MSPeakLists = plistsEmpty), NA)
    expect_error(reportPDF(fGroupsEmpty, getWorkPath(), reportFGroups = FALSE, formulas = formsGF,
                           MSPeakLists = plists), NA)
    expect_reportHTML(makeReportHTML(fGroups, formulas = formsGFEmpty, MSPeakLists = plistsEmpty))
    expect_error(makeReportHTML(fGroupsEmpty, formulas = formsGF, MSPeakLists = plists), NA)
})

anPLGroup2 <- screenInfo(fGroups)[name == "DEET"]$group
test_that("plotting works", {
    expect_doppel("form-spec", function() plotSpectrum(formsGFWithMSMS, index = 1, anPLGroup, MSPeakLists = plists))
    expect_doppel("form-spec_sim", function() plotSpectrum(formsGFWithMSMS, index = c(1, 1), c(anPLGroup, anPLGroup2),
                                                           MSPeakLists = plists))

    expect_doppel("form-scores", function() plotScores(formsGFWithMSMS, index = 1, anPLGroup))

    skip_if_not(doSIRIUS)
    expect_doppel("venn", function() plotVenn(formsGF, formsSIR))
    expect_error(plotVenn(formsGFEmpty, formsSIREmpty))
    expect_equal(expect_plot(plotVenn(formsGF, formsSIR))$areas[2], length(formsSIR))
    expect_equal(expect_plot(plotVenn(formsGF, formsSIREmpty))$areas[1], length(formsGF))
    expect_equal(expect_plot(plotVenn(formsGFEmpty, formsSIR))$areas[2], length(formsSIR))
    expect_equal(expect_plot(plotVenn(formsGF, formsSIREmpty))$intersectionCounts, 0)

    expect_ggplot(plotUpSet(formsGF, formsSIR))
    expect_error(plotUpSet(formsGFEmpty, formsSIREmpty))
    expect_error(plotUpSet(formsGF, formsSIREmpty))

    expect_equal(expect_plot(plotVenn(formsGF, formsSIR))$intersectionCounts,
                 length(doFormCons(formsGF, formsSIR, MSPeakLists = plists, relMinAbundance = 1)))
})

if (testWithSets())
{
    fgOneEmptySet <- makeOneEmptySetFGroups(fGroups)
    formsGFOneEmptySet <- doGenForms(fgOneEmptySet, plists, "genform")
    formsGFAvgSpecCols <- doGenForms(fGroups, plists, "genform", setAvgSpecificScores = TRUE)
}

test_that("sets functionality", {
    skip_if_not(testWithSets())
    
    expect_equal(formsGF, formsGF[, sets = sets(formsGF)])
    expect_length(formsGF[, sets = character()], 0)
    expect_equal(sets(filter(formsGF, sets = "positive", negate = TRUE)), "negative")
    expect_setequal(groupNames(formsGF), unique(unlist(lapply(setObjects(formsGF), groupNames))))
    expect_setequal(groupNames(unset(formsGF, "positive")), groupNames(setObjects(formsGF)[[1]]))
    expect_setequal(groupNames(unset(formsGFOneEmptySet, "positive")), groupNames(setObjects(formsGFOneEmptySet)[[1]]))
    expect_length(unset(formsGFOneEmptySet, "negative"), 0)
    
    expect_length(doGenForms(fgOneEmptySet, plists, "genform", setThreshold = 1), 0)
    expect_length(doGenForms(fgOneEmptySet, plists, "genform", setThresholdAnn = 1), length(formsGFOneEmptySet))
    
    # setAvgSpecificScores=FALSE (default)
    checkmate::expect_names(names(as.data.table(formsGF)),
                            must.include = paste0("isoScore-", sets(formsGFAvgSpecCols)[1]))
    checkmate::expect_names(names(as.data.table(formsGFAvgSpecCols)), must.include = "isoScore")
    
    expect_doppel("form-spec-set", function() plotSpectrum(formsGFWithMSMS, index = 1, anPLGroup, MSPeakLists = plists,
                                                           perSet = FALSE))
    expect_doppel("form-spec-set-perset", function() plotSpectrum(formsGFWithMSMS, index = 1, anPLGroup,
                                                                  MSPeakLists = plists, perSet = TRUE, mirror = FALSE))
    expect_doppel("form-spec-set-mirror", function() plotSpectrum(formsGFWithMSMS, index = 1, anPLGroup,
                                                                  MSPeakLists = plists, perSet = TRUE, mirror = TRUE))
    
    skip_if_not(doSIRIUS)
    
    expect_gt(length(setObjects(formsSIRFPs)[[1]]@fingerprints), 0)
    expect_gt(length(setObjects(formsSIRFPs)[[2]]@fingerprints), 0)
    testSIRFPSubset(setObjects(formsSIRFPs)[[1]])
})
