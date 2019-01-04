context("compounds")

fGroups <- getTestFGroups(getTestAnaInfo()[4, ])

# convert to screening results to simplify things a bit
fGroups <- groupFeaturesScreening(fGroups, screenTargets(fGroups, patRoonData::targets))

mfTestDBPath <- file.path(getTestDataPath(), "test-mf-db.csv")
mfTestDB <- fread(mfTestDBPath)

# just focus on 5 targets, these are named exactly the same as in the MetFrag test DB
fGroupsSub <- fGroups[, mfTestDB$Name]

plists <- generateMSPeakLists(fGroupsSub, "mzr")
plistsEmpty <- getEmptyPLists()
plistsEmptyMS <- removeMSPlists(plists, "MS")
fGroupsEmpty <- getEmptyTestFGroups()

doMetFrag <- !is.null(getOption("patRoon.path.metFragCL")) && nzchar(getOption("patRoon.path.metFragCL"))
doSIRIUS <- !is.null(getOption("patRoon.path.SIRIUS")) && nzchar(getOption("patRoon.path.SIRIUS"))

callMF <- function(fGroups, plists, db = mfTestDBPath, to = 300)
{
    generateCompounds(fGroups, plists, "metfrag", logPath = NULL,
                      adduct = 1, isPositive = TRUE, timeout = to,
                      database = "LocalCSV", scoreTypes = "FragmenterScore",
                      extraOpts = list(LocalDatabasePath = db))
}

if (doMetFrag)
{
    compsMF <- callMF(fGroupsSub, plists)
    ct <- compoundTable(compsMF)
    compsMFEmpty <- callMF(fGroupsEmpty, plists = plistsEmpty)
    compsMFEmptyPL <- callMF(fGroupsSub, plists = plistsEmpty)
}

if (doSIRIUS)
{
    compsSIR <- generateCompounds(fGroupsSub, plists, "sirius", logPath = NULL)
    compsSIREmpty <- generateCompounds(fGroupsEmpty, plistsEmpty, "sirius", logPath = NULL)
    compsSIREmptyPL <- generateCompounds(fGroupsSub, plistsEmpty, "sirius", logPath = NULL)
}

test_that("verify MetFragCL compound generation", {
    skip_if_not(doMetFrag)
    expect_known_value(compsMF, testFile("compounds-mf"))
    expect_known_show(compsMF, testFile("compounds-mf", text = TRUE))
    expect_length(compsMF, 5) # should be one compound per feature group
    # make sure that all feature group names (=targets) correspond to identified compounds
    expect_true(all(sapply(names(ct), function(grp) nrow(ct[[grp]]) == 1 && ct[[grp]]$identifier == grp)))
    expect_length(compsMFEmpty, 0)
    expect_length(compsMFEmptyPL, 0)
})

test_that("verify SIRIUS compound generation", {
    skip_if_not(doSIRIUS)
    expect_known_value(compsSIR, testFile("compounds-sir"))
    expect_known_show(compsSIR, testFile("compounds-sir", text = TRUE))
    expect_length(compsSIREmpty, 0)
    expect_length(compsSIREmptyPL, 0)
    expect_length(generateCompounds(fGroupsSub, plistsEmptyMS, "sirius", logPath = NULL), 0)
})

hasCompounds <- doMetFrag || doSIRIUS

if (doMetFrag)
{
    # include some isomers to test filtering... (sirius should already have multiple compounds for feature groups)
    compsMFIso <- callMF(fGroupsSub, plists, db = file.path(getTestDataPath(), "test-mf-db-isomers.csv"))
}

# continue with one or another...
comps <- if (doMetFrag) compsMFIso else if (doSIRIUS) compsSIR
compsEmpty <- if (doMetFrag) compsMFEmptyPL else if (doSIRIUS) compsSIREmptyPL

test_that("filtering works", {
    skip_if_not(hasCompounds)

    expect_lte(length(filter(comps, topMost = 1)), length(fGroupsSub))
    expect_lte(length(filter(comps, topMost = 5)), 5 * length(fGroupsSub))
    expect_lte(length(filter(comps, minExplainedPeaks = 2)), length(comps))
    expect_length(filter(comps, minExplainedPeaks = 1E6), 0)
    expect_length(filter(compsEmpty, minExplainedPeaks = 2, topMost = 1), 0)

    skip_if_not(doMetFrag)
    expect_lt(length(filter(compsMFIso, minScore = 2)), length(compsMFIso))
    expect_lt(length(filter(compsMFIso, minFragScore = 200)), length(compsMFIso))

    skip_if_not(doSIRIUS)
    expect_lt(length(filter(compsSIR, minScore = -200)), length(compsSIR))
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
})

if (doMetFrag)
{
    forms <- generateFormulas(fGroupsSub, "genform", plists)
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
                expect_warning(compsJNI <- callMF(fGroupsSub, plists), NA, info = info)
                expect_equal(compsJNI, compsMF, info = info)
            }
        })
    })
})

test_that("MetFrag can timeout", {
    skip_if_not(doMetFrag)
    withr::with_options(c(patRoon.cache.mode = "none"), {
        # call with unreasonably short timeout...
        expect_warning(compsTO <- callMF(fGroupsSub, plists, to = 1))
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
    expect_lt(length(consensus(compsMF, compsSIR, compThreshold = 1)), length(compsCons))
    expect_length(consensus(compsMFEmptyPL, compsSIREmptyPL), 0)
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

    expect_reportMD(makeReportMD(fGroups, reportPlots = "none",
                                 compounds = comps, MSPeakLists = plists))
})

test_that("reporting empty objects works", {
    skip_if_not(hasCompounds)
    expect_error(reportCSV(fGroups, getWorkPath(), compounds = compsEmpty), NA)
    expect_error(reportPDF(fGroups, getWorkPath(), reportFGroups = FALSE, compounds = compsEmpty,
                           MSPeakLists = plistsEmpty), NA)
    expect_reportMD(makeReportMD(fGroups, reportPlots = "none",
                                 compounds = compsEmpty, MSPeakLists = plistsEmpty))
})

test_that("plotting works", {
    skip_if_not(doMetFrag)

    # plotting structure seems to be difficult to do reproducible between systems, so disable for vdiffr now...
    expect_doppel("compound-spec", function() plotSpec(compsMFIso, 1, names(compoundTable(compsMFIso))[1], plists, plotStruct = FALSE))
    expect_plot(plotSpec(compsMFIso, 1, names(compoundTable(compsMFIso))[1], plists, plotStruct = TRUE))
    # expect_doppel("spec-gg", plotSpec(compsMFIso, 1, names(compoundTable(compsMFIso))[1], plists, useGGPlot2 = TRUE))
    expect_plot(print(plotSpec(compsMFIso, 1, names(compoundTable(compsMFIso))[1], plists, useGGPlot2 = TRUE)))

    # plotStructure gives an empty plot??
    # expect_doppel("struct", function() plotStructure(compsMFIso, 1, names(compoundTable(compsMFIso))[1]))
    expect_plot(plotStructure(compsMFIso, 1, names(compoundTable(compsMFIso))[1]))
    expect_doppel("scores", function() plotScores(compsMFIso, 1, names(compoundTable(compsMFIso))[1]))
})
