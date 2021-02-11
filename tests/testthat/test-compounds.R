context("compounds")

fGroups <- getCompFGroups()

plists <- generateMSPeakLists(fGroups, "mzr")
plistsEmpty <- plists[FALSE, reAverage = TRUE]
plistsEmptyMS <- removeMSPlists(plists, "MS")
fGroupsEmpty <- getEmptyTestFGroups()

doMetFrag <- !is.null(getOption("patRoon.path.MetFragCL")) && nzchar(getOption("patRoon.path.MetFragCL"))
doSIRIUS <- !is.null(getOption("patRoon.path.SIRIUS")) && nzchar(getOption("patRoon.path.SIRIUS"))

if (doMetFrag)
{
    compsMF <- callMF(fGroups, plists)
    ct <- compoundTable(compsMF)
    compsMFEmpty <- callMF(fGroupsEmpty, plists = plistsEmpty)
    compsMFEmptyPL <- callMF(fGroups, plists = plistsEmpty)
}

if (doSIRIUS)
{
    compsSIR <- doGenComps(fGroups, plists, "sirius")
    compsSIREmpty <- doGenComps(fGroupsEmpty, plistsEmpty, "sirius")
    compsSIREmptyPL <- doGenComps(fGroups, plistsEmpty, "sirius")
}

test_that("verify MetFragCL compound generation", {
    skip_if_not(doMetFrag)
    expect_known_value(compsMF, testFile("compounds-mf"))
    expect_known_show(compsMF, testFile("compounds-mf", text = TRUE))
    expect_length(compsMF, 5) # should be one compound per feature group
    # make sure that all suspect names correspond to identified compounds
    expect_true(all(sapply(names(ct), function(grp) nrow(ct[[grp]]) == 1 && ct[[grp]]$identifier == screenInfo(fGroups)[group == grp]$name)))
    expect_length(compsMFEmpty, 0)
    expect_length(compsMFEmptyPL, 0)
})

test_that("verify SIRIUS compound generation", {
    skip_if_not(doSIRIUS)
    expect_known_value(compsSIR, testFile("compounds-sir"))
    expect_known_show(compsSIR, testFile("compounds-sir", text = TRUE))
    expect_length(compsSIREmpty, 0)
    expect_length(compsSIREmptyPL, 0)
    expect_length(doGenComps(fGroups, plistsEmptyMS, "sirius"), 0)
})

hasCompounds <- doMetFrag || doSIRIUS

if (doMetFrag)
{
    # include some isomers to test filtering... (sirius should already have multiple compounds for feature groups)
    compsMFIso <- callMF(fGroups, plists, db = file.path(getTestDataPath(), "test-mf-db-isomers.csv"))
}

# continue with one or another...
comps <- if (doMetFrag) compsMFIso else if (doSIRIUS) compsSIR
compsEmpty <- if (doMetFrag) compsMFEmptyPL else if (doSIRIUS) compsSIREmptyPL
compsExplained <- filter(comps, minExplainedPeaks = 1)

test_that("filtering works", {
    skip_if_not(hasCompounds)

    expect_lte(length(filter(comps, topMost = 1)), length(fGroups))
    expect_lte(length(filter(comps, topMost = 5)), 5 * length(fGroups))
    expect_lte(length(filter(comps, topMost = 1, negate = TRUE)), length(fGroups))
    expect_lte(length(filter(comps, topMost = 5, negate = TRUE)), 5 * length(fGroups))
    expect_true(all(as.data.table(filter(comps, topMost = 1))$score >
                        as.data.table(filter(comps, topMost = 1, negate = TRUE))$score))

    expect_lte(length(filter(comps, minExplainedPeaks = 2)), length(comps))
    expect_length(filter(comps, minExplainedPeaks = 1E6), 0)
    expect_length(filter(compsEmpty, minExplainedPeaks = 2, topMost = 1), 0)
    expect_equivalent(filter(comps, scoreLimits = list(fragScore = c(-Inf, Inf))), comps)
    expect_lte(length(filter(comps, minExplainedPeaks = 2, negate = TRUE)), length(comps))
    expect_equivalent(filter(comps, minExplainedPeaks = 1E6, negate = TRUE), comps)
    expect_length(filter(compsEmpty, minExplainedPeaks = 2, topMost = 1, negate = TRUE), 0)
    expect_length(filter(comps, scoreLimits = list(fragScore = c(-Inf, Inf)), negate = TRUE), 0)

    expect_length(filter(comps, elements = "C1-100"), length(comps)) # all should contain carbon
    expect_length(filter(comps, elements = c("Na1-100", "C1-100")), length(comps)) # no sodium, but carbon should be there
    expect_length(filter(comps, elements = c("H1-100", "C1-100")), length(comps)) # presence of both shouldn't affect results
    expect_length(filter(comps, elements = "Na1-100"), 0) # no sodium
    expect_length(filter(comps, elements = "Na0-100"), length(comps)) # no sodium, but optional
    expect_length(filter(comps, elements = "C1-100", negate = TRUE), 0)
    expect_length(filter(comps, elements = c("Na1-100", "C1-100"), negate = TRUE), length(comps))
    expect_length(filter(comps, elements = c("H1-100", "C1-100"), negate = TRUE), 0)
    expect_length(filter(comps, elements = "Na1-100", negate = TRUE), length(comps))
    expect_length(filter(comps, elements = "Na0-100", negate = TRUE), 0)

    expect_lte(length(filter(compsExplained, fragElements = "C1-100")), length(compsExplained))
    expect_length(filter(compsExplained, fragElements = "C"), 0) # fragments may contain only single carbon
    expect_length(filter(compsExplained, fragElements = "Na1-100"), 0)
    expect_length(filter(compsExplained, fragElements = "Na0-100"), length(compsExplained))
    expect_lte(length(filter(compsExplained, fragElements = "C1-100", negate = TRUE)), length(compsExplained))
    expect_length(filter(compsExplained, fragElements = "C", negate = TRUE), length(compsExplained))
    expect_length(filter(compsExplained, fragElements = "Na1-100", negate = TRUE), length(compsExplained))
    expect_length(filter(compsExplained, fragElements = "Na0-100", negate = TRUE), 0)

    expect_length(filter(compsExplained, lossElements = "C0-100"), length(compsExplained))
    expect_length(filter(compsExplained, lossElements = "Na0-100"), length(compsExplained))
    expect_gt(length(filter(compsExplained, lossElements = "C1-100")), 0) # NL might be empty, at least some should contain carbon though!
    expect_length(filter(compsExplained, lossElements = "Na1-100"), 0) # no sodium
    expect_length(filter(compsExplained, lossElements = "C0-100", negate = TRUE), 0)
    expect_length(filter(compsExplained, lossElements = "Na0-100", negate = TRUE), 0)
    expect_gt(length(filter(compsExplained, lossElements = "C1-100", negate = TRUE)), 0)
    expect_length(filter(compsExplained, lossElements = "Na1-100", negate = TRUE), length(compsExplained))

    skip_if_not(doMetFrag)
    expect_lt(length(filter(compsMFIso, minScore = 0.75)), length(compsMFIso))
    expect_lt(length(filter(compsMFIso, minScore = 0.75, negate = TRUE)), length(compsMFIso))
    expect_lt(length(filter(compsMFIso, minFragScore = 200)), length(compsMFIso))
    expect_lt(length(filter(compsMFIso, minFragScore = 200, negate = TRUE)), length(compsMFIso))
    expect_lt(length(filter(compsMFIso, scoreLimits = list(fragScore = c(200, Inf)), negate = TRUE)),
              length(compsMFIso))

    skip_if_not(doSIRIUS)
    expect_lt(length(filter(compsSIR, minScore = -200)), length(compsSIR))
    expect_lt(length(filter(compsSIR, minScore = -200, negate = TRUE)), length(compsSIR))
    expect_lt(length(filter(compsSIR, scoreLimits = list(score = c(-200, Inf)))), length(compsSIR))
    expect_lt(length(filter(compsSIR, scoreLimits = list(score = c(-200, Inf)), negate = TRUE)),
              length(compsSIR))
})

test_that("basic subsetting", {
    skip_if_not(hasCompounds)

    expect_length(comps["nope"], 0)
    expect_equivalent(groupNames(comps[1:2]), groupNames(comps)[1:2])
    expect_equivalent(groupNames(comps[groupNames(comps)[2:3]]), groupNames(comps)[2:3])
    expect_equivalent(groupNames(comps[c(FALSE, TRUE)]), groupNames(comps)[c(FALSE, TRUE)])
    expect_equal(length(comps[FALSE]), 0)
    expect_length(compsEmpty[1:5], 0)

    expect_equivalent(comps[[1]], compoundTable(comps)[[1]])
    expect_equivalent(comps[[groupNames(comps)[1]]], compoundTable(comps)[[1]])
    expect_equivalent(callDollar(comps, groupNames(comps)[1]), comps[[1]])
})

test_that("basic usage", {
    expect_equal(nrow(as.data.table(comps)), length(comps))
    checkmate::expect_names(names(as.data.table(comps, fGroups = fGroups)),
                            must.include = c("ret", "group_mz"))
    checkmate::expect_names(names(as.data.table(comps, fragments = TRUE)),
                            must.include = c("frag_formula", "frag_mz"))
    expect_gt(nrow(as.data.table(comps, fragments = TRUE)), length(comps))
    expect_range(as.data.table(comps, normalizeScores = "max")$fragScore, c(0, 1))
})

if (doMetFrag)
{
    forms <- doGenForms(fGroups, "genform", plists)
    compsMFIsoF <- addFormulaScoring(compsMFIso, forms)
}

test_that("formula scoring works", {
    skip_if_not(doMetFrag)
    expect_lt(length(filter(compsMFIsoF, minFormulaScore = 3)), length(compsMFIsoF))
    expect_error(addFormulaScoring(compsMFEmptyPL, forms), NA)
})

# on a clean system, i.e. where ~/.jnati/repo/jnniinchi is not yet initialized, starting multiple
# MetFrag processes in parallel (i.e. when maxProcAmount>1) may result in errors. This should now
# be fixed by setting a small delay between starting up processes (delayBetweenProc arg of executeMultiProcess())
jnatiTestDir <- file.path(tempdir(), "jnati-test")
test_that("MetFrag uninitialized jniinchi workaround", {
    skip_if_not(doMetFrag)

    # temporarily change jnati workdir so it can be safely wiped
    withr::with_envvar(c(JAVA_TOOL_OPTIONS = sprintf("-Djnati.dir=%s", jnatiTestDir)), {
        withr::with_options(c(patRoon.cache.mode = "none"), {
            for (n in seq_len(5))
            {
                info = sprintf("iter: %d", n)
                unlink(jnatiTestDir, recursive = TRUE)
                expect_warning(compsJNI <- callMF(fGroups, plists), NA, info = info)
                expect_equal(compsJNI, compsMF, info = info)
            }
        })
    })
})

test_that("MetFrag can timeout", {
    skip_if_not(doMetFrag)
    withr::with_options(c(patRoon.cache.mode = "none"), {
        # call with unreasonably short timeout...
        expect_warning(compsTO <- callMF(fGroups, plists, to = 1))
        expect_lt(length(compsTO), length(compsMF))
    })
})

if (doMetFrag && doSIRIUS)
    compsCons <- consensus(compsMF, compsSIR)

test_that("consensus works", {
    skip_if_not(hasCompounds)
    expect_length(consensus(comps, compsEmpty), length(comps))

    skip_if_not(doMetFrag && doSIRIUS)
    expect_known_value(compsCons, testFile("compounds-cons"))
    expect_known_show(compsCons, testFile("compounds-cons", text = TRUE))
    expect_lt(length(consensus(compsMF, compsSIR, relMinAbundance = 1)), length(compsCons))
    expect_length(consensus(compsMFEmptyPL, compsSIREmptyPL), 0)

    expect_equal(length(consensus(compsMF, compsSIR, uniqueFrom = 1)) +
                 length(consensus(compsMF, compsSIR, uniqueFrom = 2)) +
                 length(consensus(compsMF, compsSIR, relMinAbundance = 1)), length(compsCons))
    expect_equal(length(consensus(compsMF, compsSIR, uniqueFrom = 1:2, uniqueOuter = TRUE)) +
                 length(consensus(compsMF, compsSIR, relMinAbundance = 1)), length(compsCons))
    expect_length(consensus(compsMF, compsSIR, uniqueFrom = 1:2), length(compsCons))
    expect_lt(length(consensus(compsMF, compsSIR, uniqueFrom = 1:2, uniqueOuter = TRUE)), length(compsCons))
    expect_length(consensus(compsMFEmptyPL, compsSIREmptyPL, uniqueFrom = 1), 0)
    expect_length(consensus(compsMFEmptyPL, compsSIREmptyPL, uniqueFrom = 1, uniqueOuter = TRUE), 0)
})

if (doMetFrag && doSIRIUS)
{
    anPL <- annotatedPeakList(compsCons, index = 1, groupName = groupNames(compsCons)[1],
                              MSPeakLists = plists, formulas = forms)
    anPLOnly <- annotatedPeakList(compsCons, index = 1, groupName = groupNames(compsCons)[1],
                                  MSPeakLists = plists, formulas = forms, onlyAnnotated = TRUE)
}

test_that("annotation works", {
    skip_if_not(doMetFrag && doSIRIUS)

    expect_lt(nrow(anPLOnly), nrow(anPL))
    expect_true(any(is.na(anPL$formula)))
    expect_false(any(is.na(anPLOnly$formula)))
    expect_true(all(compsCons[[1]]$fragInfo[[1]]$formula %in% anPLOnly$formula))
    expect_true(any(grepl("metfrag", anPLOnly$mergedBy)))
    expect_true(any(grepl("sirius", anPLOnly$mergedBy)))
    expect_true(any(grepl("genform", anPLOnly$mergedBy)))
})

test_that("reporting works", {
    skip_if_not(hasCompounds)

    expect_error(reportCSV(fGroups, getWorkPath(), compounds = comps), NA)
    for (grp in names(compoundTable(comps)))
        checkmate::expect_file_exists(getWorkPath("compounds", sprintf("%s-%s.csv", class(fGroups), grp)))

    expect_error(reportPDF(fGroups, getWorkPath(), reportFGroups = FALSE, compounds = comps,
                           MSPeakLists = plists), NA)
    for (grp in names(compoundTable(comps)))
        checkmate::expect_file_exists(getWorkPath("compounds", sprintf("%s-%s.pdf", class(fGroups), grp)))

    expect_reportHTML(makeReportHTML(fGroups, reportPlots = "none",
                                     compounds = comps, MSPeakLists = plists))
})

test_that("reporting empty objects works", {
    skip_if_not(hasCompounds)
    expect_error(reportCSV(fGroups, getWorkPath(), compounds = compsEmpty), NA)
    expect_error(reportCSV(fGroupsEmpty, getWorkPath(), compounds = comps), NA)
    expect_error(reportPDF(fGroups, getWorkPath(), reportFGroups = FALSE, compounds = compsEmpty,
                           MSPeakLists = plistsEmpty), NA)
    expect_error(reportPDF(fGroupsEmpty, getWorkPath(), reportFGroups = FALSE, compounds = comps,
                           MSPeakLists = plists), NA)
    expect_reportHTML(makeReportHTML(fGroups, reportPlots = "none",
                                     compounds = compsEmpty, MSPeakLists = plistsEmpty))
    expect_error(makeReportHTML(fGroupsEmpty, reportPlots = "none",
                                compounds = comps, MSPeakLists = plists), NA)
})

test_that("plotting works", {
    skip_if_not(doMetFrag)

    # plotting structure seems to be difficult to do reproducible between systems, so disable for vdiffr now...
    expect_doppel("compound-spec", function() plotSpectrum(compsMFIso, 1, names(compoundTable(compsMFIso))[2], plists, plotStruct = FALSE))
    expect_plot(plotSpectrum(compsMFIso, 1, names(compoundTable(compsMFIso))[2], plists, plotStruct = TRUE))
    # expect_doppel("spec-gg", plotSpectrum(compsMFIso, 1, names(compoundTable(compsMFIso))[1], plists, useGGPlot2 = TRUE))
    expect_ggplot(plotSpectrum(compsMFIso, 1, names(compoundTable(compsMFIso))[2], plists, useGGPlot2 = TRUE))

    # plotStructure gives an empty plot??
    # expect_doppel("struct", function() plotStructure(compsMFIso, 1, names(compoundTable(compsMFIso))[1]))
    expect_plot(plotStructure(compsMFIso, 1, names(compoundTable(compsMFIso))[2]))
    expect_doppel("scores", function() plotScores(compsMFIso, 1, names(compoundTable(compsMFIso))[2]))

    skip_if_not(doSIRIUS)
    expect_doppel("venn", function() plotVenn(compsMF, compsSIR))
    expect_error(plotVenn(compsMFEmpty, compsSIREmpty))
    expect_equal(expect_plot(plotVenn(compsMF, compsSIR))$areas[2], length(compsSIR))
    expect_equal(expect_plot(plotVenn(compsMF, compsSIREmpty))$areas[1], length(compsMF))
    expect_equal(expect_plot(plotVenn(compsMFEmpty, compsSIR))$areas[2], length(compsSIR))
    expect_equal(expect_plot(plotVenn(compsMF, compsSIR))$intersectionCounts,
                 length(consensus(compsMF, compsSIR, relMinAbundance = 1)))
    expect_equal(expect_plot(plotVenn(compsMF, compsSIREmpty))$intersectionCounts, 0)

    expect_ggplot(plotUpSet(compsMF, compsSIR))
    expect_error(plotUpSet(compsMFEmpty, compsSIREmpty))
    expect_error(plotUpSet(compsMF, compsSIREmpty))
})
