context("compounds")

fGroups <- getTestFGroups(getTestAnaInfo()[4, ])

# convert to screening results to simplify things a bit
fGroups <- groupFeaturesScreening(fGroups, screenTargets(fGroups, patRoonData::targets))

mfTestDBPath <- file.path(getTestDataPath(), "test-mf-db.csv")
mfTestDB <- fread(mfTestDBPath)

# just focus on 5 targets, these are named exactly the same as in the MetFrag test DB
fGroupsSub <- fGroups[, mfTestDB$Name]

plists <- generateMSPeakLists(fGroupsSub, "mzr")

doMetFrag <- !is.null(getOption("patRoon.path.metFragCL")) && nzchar(getOption("patRoon.path.metFragCL"))
doSIRIUS <- !is.null(getOption("patRoon.path.SIRIUS")) && nzchar(getOption("patRoon.path.SIRIUS"))

if (doMetFrag)
{
    compsMF <- generateCompounds(fGroupsSub, plists, "metfrag", logPath = NULL,
                                 adduct = 1, isPositive = TRUE, database = "LocalCSV",
                                 scoreTypes = "FragmenterScore",
                                 extraOpts = list(LocalDatabasePath = mfTestDBPath))
    ct <- compoundTable(compsMF)
}

if (doSIRIUS)
    compsSIR <- generateCompounds(fGroupsSub, plists, "sirius", logPath = NULL)

test_that("verify MetFragCL compound generation", {
    skip_if_not(doMetFrag)
    expect_known_value(compsMF, testFile("compounds-mf"))
    expect_known_show(compsMF, testFile("compounds-mf", text = TRUE))
    expect_length(compsMF, 5) # should be one compound per feature group
    # make sure that all feature group names (=targets) correspond to identified compounds
    expect_true(all(sapply(names(ct), function(grp) nrow(ct[[grp]]) == 1 && ct[[grp]]$identifier == grp)))
})

test_that("verify SIRIUS compound generation", {
    skip_if_not(doSIRIUS)
    expect_known_value(compsSIR, testFile("compounds-sir"))
    expect_known_show(compsSIR, testFile("compounds-sir", text = TRUE))
})

hasCompounds <- doMetFrag || doSIRIUS

if (doMetFrag)
{
    # include some isomers to test filtering... (sirius should already have multiple compounds for feature groups)
    compsMFIso <- generateCompounds(fGroupsSub, plists, "metfrag", logPath = NULL,
                                    adduct = 1, isPositive = TRUE, database = "LocalCSV",
                                    scoreTypes = "FragmenterScore",
                                 extraOpts = list(LocalDatabasePath = file.path(getTestDataPath(), "test-mf-db-isomers.csv")))
}

# continue with one or another...
comps <- if (doMetFrag) compsMFIso else if (doSIRIUS) compsSIR

test_that("filtering works", {
    skip_if_not(hasCompounds)
    
    expect_length(filter(comps, topMost = 1), length(fGroupsSub))
    expect_lte(length(filter(comps, topMost = 5)), 5 * length(fGroupsSub))
    expect_lt(length(filter(comps, minExplainedPeaks = 2)), length(comps))
    
    skip_if_not(doMetFrag)
    expect_lt(length(filter(compsMFIso, minScore = 2)), length(compsMFIso))
    expect_lt(length(filter(compsMFIso, minFragScore = 200)), length(compsMFIso))
    
    skip_if_not(doSIRIUS)
    expect_lt(length(filter(compsSIR, minScore = -200)), length(compsSIR))
})

if (doMetFrag)
{
    forms <- consensus(generateFormulas(fGroupsSub, "genform", plists), fGroups = fGroupsSub)
    compsMFIsoF <- addFormulaScoring(compsMFIso, forms)
}

test_that("formula scoring works", {
    skip_if_not(doMetFrag)
    expect_lt(length(filter(compsMFIsoF, minFormulaScore = 3)), length(compsMFIsoF))
})

if (doMetFrag && doSIRIUS)
    compsCons <- consensus(compsMF, compsSIR)

test_that("consensus works", {
    skip_if_not(doMetFrag && doSIRIUS)
    expect_known_value(compsCons, testFile("compounds-cons"))
    expect_known_show(compsCons, testFile("compounds-cons", text = TRUE))
    expect_lt(length(consensus(compsMF, compsSIR, compThreshold = 1)), length(compsCons))
})

test_that("feature group filtering", {
    skip_if_not(hasCompounds)
    expect_named(filter(fGroups, compounds = comps), names(compoundTable(comps)))
})
