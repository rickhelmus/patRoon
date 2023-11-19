context("MS peak lists")

fGroups <- getTestFGroups(getTestAnaInfoAnn())[, 1:100]
plists <- generateMSPeakLists(fGroups, "mzr")
plistsMSMS <- filter(plists, withMSMS = TRUE)

if (doDATests())
{
    # NOTE: use different analyses than first: this call will clearout all DA
    # compounds, forcing other calls (from features/formulas tests) to re-run
    # FMF
    fgDA <- getTestFGroupsDA(getDAAnaInfo("std2"))[, 1:25]

    # NOTE: set bgsubtr to FALSE: subtraction might remove pecursor peaks of
    # (often wrong) low intensity features and result in warnings during
    # averaging
    plistsDA <- generateMSPeakLists(fgDA, "bruker", save = FALSE, bgsubtr = FALSE)
    plistsDAEmpty <- generateMSPeakLists(fgDA[, "nope"], "bruker", save = FALSE, bgsubtr = FALSE)

    fgDA2 <- getTestFGroupsDA(getDAAnaInfo("std1"))[, 1:25]
    plistsDAFMF <- generateMSPeakLists(fgDA2, "brukerfmf")
}

# remove Ion Mobility data as it inconsistently is present or not
plistsNoIM <- plists
plistsNoIM@metadata <- lapply(plistsNoIM@metadata, function(mda) lapply(mda, function(mdf) lapply(mdf, function(mds) mds[, setdiff(names(mds), "ionMobilityDriftTime"), with = FALSE])))
if (testWithSets())
    plistsNoIM@analysisInfo <- data.table() # remove as it is system dependent

test_that("verify generation of MS peak lists", {
    expect_known_value(plistsNoIM, testFile("plists-mzr"))

    skip_if_not(doDATests())
    expect_known_value(plistsDA, testFile("plists-DA"))
    expect_known_value(plistsDAFMF, testFile("plists-DAFMF"))
})

test_that("verify show output", {
    expect_known_show(plists, testFile("plists-mzr", text = TRUE))

    skip_if_not(doDATests())
    expect_known_show(plistsDA, testFile("plists-DA", text = TRUE))
    expect_known_show(plistsDAFMF, testFile("plists-DAFMF", text = TRUE))
})

checkIntLimit <- function(plists, relative, doMin, doMSMS, plistsOrig = NULL)
{
    plists <- removePrecursors(plists)

    lim <- if (doMin) min else max

    intPLLimit <- function(pl, plorig = NULL)
    {
        spec <- specorig <- NULL
        if (!doMSMS && !is.null(pl[["MS"]]))
        {
            spec <- pl[["MS"]]
            if (!is.null(plorig))
                specorig <- plorig[["MS"]]
        }
        else if (doMSMS && !is.null(pl$MSMS))
        {
            spec <- pl[["MSMS"]]
            if (!is.null(plorig))
                specorig <- plorig[["MSMS"]]
        }
        else
            return(NA)

        if (is.null(spec) || nrow(spec) == 0)
            return(NA)

        ret <- lim(spec$intensity)
        if (relative)
        {
            if (!is.null(spec[["set"]]))
            {
                spec <- copy(spec)
                spec[, intrel := intensity / max(specorig[set == .BY]$intensity), by = "set"]
                ret <- lim(spec$intrel)
            }
            else
                ret <- ret / max(specorig$intensity)
        }
        return(ret)
    }

    ftLim <- lim(sapply(names(plists@peakLists), function(ana) lim(sapply(names(plists@peakLists[[ana]]), function(grp)
    {
        intPLLimit(plists[[ana, grp]], if (!is.null(plistsOrig)) plistsOrig[[ana, grp]] else NULL)
    }), na.rm = TRUE)), na.rm = TRUE)

    fgLim <- lim(sapply(names(plists@averagedPeakLists), function(grp)
    {
        intPLLimit(plists[[grp]], if (!is.null(plistsOrig)) plistsOrig[[grp]] else NULL)
    }), na.rm = TRUE)


    return(lim(ftLim, fgLim))
}

checkPeaksLimit <- function(plists, doMin, doMSMS)
{
    lim <- if (doMin) min else max
    pl <- peakLists(removePrecursors(plists))
    return(lim(sapply(pl, function(ana) lim(sapply(ana, function(grp)
    {
        if (!doMSMS && !is.null(grp[["MS"]]))
            return(nrow(grp[["MS"]]))
        else if (doMSMS && !is.null(grp[["MSMS"]]))
            return(nrow(grp[["MSMS"]]))
        else
            return(NA)
    }), na.rm = TRUE)), na.rm = TRUE))
}

isoTestBy <- if (testWithSets()) c("group", "set") else "group"

test_that("delete and filter", {
    checkmate::expect_names(groupNames(delete(plists, i = 1)), disjunct.from = groupNames(plists)[1])
    checkmate::expect_names(groupNames(delete(plists, i = groupNames(plists)[1])), disjunct.from = groupNames(plists)[1])
    checkmate::expect_names(analyses(delete(plists, k = 1)), disjunct.from = analyses(plists)[1])
    checkmate::expect_names(analyses(delete(plists, k = analyses(plists)[1])), disjunct.from = analyses(plists)[1])
    expect_length(delete(plists, i = groupNames(plists)), 0)
    expect_length(delete(plists, k = analyses(plists)), 0)
    expect_equal(delete(plists, j = 1)[[1]]$MS$ID[1], plists[[1]]$MS$ID[2])
    expect_equal(delete(plists, j = function(...) 1)[[1]]$MS$ID[1], plists[[1]]$MS$ID[2])
    expect_length(delete(plists, j = function(...) TRUE), 0)
    expect_equal(delete(plists, i = character()), plists)
    expect_equal(delete(plists, j = integer()), plists)
    expect_equal(delete(plists, k = character()), plists)
    expect_length(delete(plists), 0)
    
    expect_gte(checkIntLimit(filter(plists, absMSIntThr = 2500), FALSE, TRUE, FALSE), 2500)
    expect_gte(checkIntLimit(filter(plists, absMSMSIntThr = 2500), FALSE, TRUE, TRUE), 2500)
    expect_lte(checkIntLimit(filter(plists, absMSIntThr = 2500, negate = TRUE), FALSE, FALSE, FALSE), 2500)
    expect_lte(checkIntLimit(filter(plists, absMSMSIntThr = 2500, negate = TRUE), FALSE, FALSE, TRUE), 2500)

    expect_gte(checkIntLimit(filter(plists, relMSIntThr = 0.2), TRUE, TRUE, FALSE, plists), 0.2)
    expect_gte(checkIntLimit(filter(plists, relMSMSIntThr = 0.2), TRUE, TRUE, TRUE, plists), 0.2)
    expect_lte(checkIntLimit(filter(plists, relMSIntThr = 0.2, negate = TRUE), TRUE, FALSE, FALSE, plists), 0.2)
    expect_lte(checkIntLimit(filter(plists, relMSMSIntThr = 0.2, negate = TRUE), TRUE, FALSE, TRUE, plists), 0.2)

    expect_lte(checkPeaksLimit(filter(plists, topMSPeaks = 5), FALSE, FALSE), 5)
    expect_lte(checkPeaksLimit(filter(plists, topMSMSPeaks = 5), FALSE, TRUE), 5)
    expect_lte(checkPeaksLimit(filter(plists, topMSPeaks = 5, negate = TRUE), FALSE, FALSE), 5)
    expect_lte(checkPeaksLimit(filter(plists, topMSMSPeaks = 5, negate = TRUE), FALSE, TRUE), 5)

    expect_true(all(sapply(averagedPeakLists(plistsMSMS), function(pl) !is.null(pl[["MSMS"]]))))
    expect_true(all(sapply(averagedPeakLists(filter(plists, withMSMS = TRUE, negate = TRUE)),
                           function(pl) is.null(pl[["MSMS"]]))))

    expect_lte(length(filter(plists, topMSMSPeaks = 5, retainPrecursorMSMS = FALSE)),
               length(filter(plists, topMSMSPeaks = 5)))

    # isotopes for MS peak lists should not exceed M+5 (default)
    expect_lte(max(as.data.table(filter(
        plists, isolatePrec = getDefIsolatePrecParams()))[type == "MS", diff(range(mz)), by = isoTestBy][["V1"]]), 5)
    # half with z=2
    expect_lte(max(as.data.table(filter(
        plists, isolatePrec = getDefIsolatePrecParams(z=2)))[type == "MS", diff(range(mz)), by = isoTestBy][["V1"]]), 2.5)

    # UNDONE: deisotope?
})

plistsEmpty <- removePrecursors(filter(plists, absMSIntThr = 1E9, absMSMSIntThr = 1E9))
plistsEmptyMS <- removePrecursors(filter(plists, absMSIntThr = 1E9))

test_that("empty object", {
    expect_length(plistsEmpty, 0)
    expect_length(delete(plistsEmpty), 0)
    expect_length(filter(plistsEmpty, relMSIntThr = 0.2, relMSMSIntThr = 0.2, topMSPeaks = 10,
                         topMSMSPeaks = 10), 0)
    expect_length(generateMSPeakLists(getEmptyTestFGroups(), "mzr"), 0)
    expect_lt(length(plistsEmptyMS), length(plists))
    expect_gt(length(plistsEmptyMS), 0)

    skip_if_not(doDATests())
    expect_length(plistsDAEmpty, 0)
})

test_that("basic functionality", {
    expect_length(plists["nope"], 0)
    expect_equivalent(analyses(plists[1:2]), analyses(fGroups)[1:2])
    expect_equivalent(analyses(plists[analyses(fGroups)[1:2]]), analyses(fGroups)[1:2])
    expect_equivalent(analyses(plists[c(FALSE, TRUE)]), analyses(fGroups)[c(FALSE, TRUE)])
    expect_equivalent(groupNames(plists[, 1:2]), groupNames(plists)[1:2])
    expect_equivalent(groupNames(plists[, groupNames(plists)[2:3]]), groupNames(plists)[2:3])
    expect_equivalent(groupNames(plists[, c(FALSE, TRUE)]), groupNames(plists)[c(FALSE, TRUE)])
    expect_equal(length(plists[FALSE, reAverage = TRUE]), 0)
    expect_length(plistsEmpty[1:5], 0)

    expect_equivalent(plists[[2, 15]], peakLists(plists)[[2]][[groupNames(plists)[15]]])
    expect_equivalent(plists[[analyses(plists)[2], groupNames(plists)[15]]], peakLists(plists)[[2]][[groupNames(plists)[15]]])

    expect_equivalent(plists[[20]], averagedPeakLists(plists)[[groupNames(plists)[20]]])
    expect_equivalent(plists[[groupNames(plists)[20]]], averagedPeakLists(plists)[[groupNames(plists)[20]]])

    expect_equal(nrow(as.data.table(plists, averaged = TRUE)),
                 sum(unlist(recursiveApplyDT(averagedPeakLists(plists), nrow))))
    expect_equal(nrow(as.data.table(plists, averaged = FALSE)), length(plists))
    checkmate::expect_names(names(as.data.table(plists, fGroups = fGroups)),
                            must.include = c("ret", "group_mz"))
})

testSpecSim <- function(obj, groupName1, groupName2, ..., expectNA = FALSE, expectOne = FALSE)
{
    sim <- spectrumSimilarity(obj, groupName1 = groupName1, groupName2 = groupName2, NAToZero = FALSE, ...)
    simNoNA <- spectrumSimilarity(obj, groupName1 = groupName1, groupName2 = groupName2, NAToZero = TRUE, ...)
    
    expect_range(simNoNA, c(0, 1))
    
    expect_length(sim, if (is.null(groupName2)) length(groupName1)^2 else length(groupName1) * length(groupName2))
    expect_true(!expectNA || any(is.na(sim)))
    expect_true(!expectOne || all(numEQ(sim, 1)))
    
    if (length(groupName1) == 1 && (is.null(groupName2) || length(groupName2) == 1))
    {
        expect_false(checkmate::testMatrix(sim)) # should be dropped
        checkmate::expect_matrix(spectrumSimilarity(obj, groupName1 = groupName1, groupName2 = groupName2, NAToZero = FALSE,
                                                    drop = FALSE, ...), nrows = 1, ncols = 1, row.names = "unique",
                                 col.names = "unique")
    }
    else
        checkmate::expect_matrix(sim, nrows = length(groupName1),
                                 ncols = length(if (is.null(groupName2)) groupName1 else groupName2),
                                 row.names = "unique", col.names = "unique")
}

simFG1 <- groupNames(plists)[2]; simFG2 <- groupNames(plists)[3]
test_that("spectral similarity", {
    testSpecSim(plists, simFG1, simFG2, analysis1 = analyses(plists)[1],
                analysis2 = analyses(plists)[2], MSLevel = 1)
    testSpecSim(plists, simFG1, simFG2, analysis1 = analyses(plists)[1],
                analysis2 = NULL, MSLevel = 1)
    testSpecSim(plists, simFG1, simFG2, analysis1 = NULL,
                analysis2 = analyses(plists)[1], MSLevel = 1)
    testSpecSim(plists, simFG1, simFG2, MSLevel = 1)

    testSpecSim(plists, simFG1, simFG1, analysis1 = analyses(plists)[1],
                analysis2 = analyses(plists)[1], MSLevel = 1, expectOne = TRUE)
    testSpecSim(plists, simFG1, simFG1, MSLevel = 1, expectOne = TRUE)
    testSpecSim(plists, simFG1, NULL, MSLevel = 1, expectOne = TRUE)
    
    testSpecSim(plists, groupNames(plists), groupNames(plists),
                analysis1 = rep(analyses(plists)[1], length(groupNames(plists))),
                analysis2 = rep(analyses(plists)[1], length(groupNames(plists))), MSLevel = 1, expectNA = testWithSets())
    testSpecSim(plists, groupNames(plists), groupNames(plists), MSLevel = 1, expectNA = testWithSets())
    testSpecSim(plists, groupNames(plists), groupNames(plists), MSLevel = 2, expectNA = TRUE)
    testSpecSim(plists, groupNames(plists)[1:5], groupNames(plists)[6:10], MSLevel = 1)
    
    testSpecSim(plists, groupNames(plists), NULL, MSLevel = 1,
                analysis1 = rep(analyses(plists)[1], length(groupNames(plists))), expectNA = testWithSets())
    testSpecSim(plists, groupNames(plists), NULL, MSLevel = 1, expectNA = testWithSets())
    testSpecSim(plists, groupNames(plists), NULL, MSLevel = 2, expectNA = TRUE)
    
    testSpecSim(plistsMSMS, groupNames(plistsMSMS)[1], groupNames(plistsMSMS)[2], MSLevel = 2)
    testSpecSim(plistsMSMS, groupNames(plistsMSMS)[1], groupNames(plistsMSMS)[1], MSLevel = 2, expectOne = TRUE)
    testSpecSim(plistsMSMS, groupNames(plistsMSMS), groupNames(plistsMSMS), MSLevel = 2)
    testSpecSim(plistsMSMS, groupNames(plistsMSMS), NULL, MSLevel = 2)
    
    expect_null(spectrumSimilarity(plistsEmpty, groupNames(plistsEmpty)[1], groupNames(plistsEmpty)[1]))
    checkmate::expect_scalar_na(spectrumSimilarity(plistsEmptyMS, groupNames(plistsEmptyMS)[1], groupNames(plistsEmptyMS)[1]))
})

makeTestPL <- function()
{
    spec <- data.table(mz = seq(50, 300, 50),
                       intensity = seq(50E4, 300E4, length.out = 6),
                       precursor = c(rep(FALSE, 5), TRUE),
                       ID = seq_len(3))
    
    ret <- MSPeakLists(algorithm = "test")
    ret@peakLists[["ana"]][["spec1"]] <- list(MS = copy(spec))
    ret@peakLists[["ana"]][["spec2"]] <- list(MS = copy(spec))
    ret@averagedPeakLists[["spec1"]] <- list(MS = copy(spec))
    ret@averagedPeakLists[["spec2"]] <- list(MS = copy(spec))
    return(ret)
}

testPL <- makeTestPL()
testPLMethod <- makeTestPL()
testPLMethod@averagedPeakLists[["spec2"]]$MS[, intensity := intensity * runif(.N, 0.1)]
testPLPrec <- makeTestPL()
testPLPrec@averagedPeakLists[["spec2"]]$MS <- testPLPrec@averagedPeakLists[["spec2"]]$MS[precursor == FALSE]
testPLInt <- makeTestPL()
testPLInt@averagedPeakLists[["spec2"]]$MS <- testPLInt@averagedPeakLists[["spec2"]]$MS[(intensity/max(intensity)) >= 0.5]
# for testing shifts: make sure precursor is changed and that there is no overlap when shift is applied by adding .1
# (UNDONE?)
testPLShiftPrec <- makeTestPL()
testPLShiftPrec@averagedPeakLists[["spec2"]]$MS[, mz := mz + 50.1]
testPLShiftBoth <- makeTestPL()
testPLShiftBoth@averagedPeakLists[["spec2"]]$MS[4:6, mz := mz + 50.1]
doSimTestPL <- function(obj, ...) spectrumSimilarity(obj, "spec1", "spec2", ...)

test_that("spectral similarity params", {
    expect_lt(doSimTestPL(testPLMethod, specSimParams = getDefSpecSimParams(method = "cosine", relMinIntensity = 0)), 1)
    expect_equal(doSimTestPL(testPLMethod, specSimParams = getDefSpecSimParams(method = "cosine", relMinIntensity = 0,
                                                                               mzWeight = 1, intWeight = 0)), 1)
    expect_equal(doSimTestPL(testPLMethod, specSimParams = getDefSpecSimParams(method = "jaccard", relMinIntensity = 0)), 1)
    expect_equal(doSimTestPL(testPLPrec, specSimParams = getDefSpecSimParams(removePrecursor = TRUE)), 1)
    expect_equal(doSimTestPL(testPLInt, specSimParams = getDefSpecSimParams(relMinIntensity = 0.5)), 1)
    expect_equal(doSimTestPL(testPL, specSimParams = getDefSpecSimParams(minPeaks = nrow(testPL[[1]]$MS))), 1)
    expect_true(is.na(doSimTestPL(testPL, specSimParams = getDefSpecSimParams(minPeaks = 100))))
    expect_lt(doSimTestPL(testPLShiftPrec, specSimParams = getDefSpecSimParams(shift = "none")), 1)
    expect_equal(doSimTestPL(testPLShiftPrec, specSimParams = getDefSpecSimParams(shift = "precursor")), 1)
    expect_equal(doSimTestPL(testPLShiftPrec, specSimParams = getDefSpecSimParams(shift = "both")), 1)
    expect_lt(doSimTestPL(testPLShiftBoth, specSimParams = getDefSpecSimParams(shift = "none")), 1)
    expect_equal(doSimTestPL(testPLShiftBoth, specSimParams = getDefSpecSimParams(shift = "both")), 1)
})

test_that("plotting works", {
    expect_doppel("mspl-spec-ms", function() plotSpectrum(plists, groupName = groupNames(plists)[70],
                                                          analysis = analyses(plists)[1], MSLevel = 1))
    expect_doppel("mspl-spec-msms", function() plotSpectrum(plistsMSMS, groupName = groupNames(plistsMSMS)[2],
                                                            analysis = analyses(plistsMSMS)[1], MSLevel = 2))
    expect_doppel("mspl-spec-avg-ms", function() plotSpectrum(plists, groupName = groupNames(plists)[32],
                                                              MSLevel = 1))
    expect_doppel("mspl-spec-avg-msms", function() plotSpectrum(plistsMSMS, groupName = groupNames(plistsMSMS)[14],
                                                                MSLevel = 2))
    
    expect_doppel("mspl-spec_sim-ms", function() plotSpectrum(plists,
                                                              groupName = c(groupNames(plistsMSMS)[2], groupNames(plistsMSMS)[3]),
                                                              analysis = c(analyses(plistsMSMS)[1], analyses(plistsMSMS)[1]),
                                                              MSLevel = 1))
    expect_doppel("mspl-spec_sim-msms", function() plotSpectrum(plistsMSMS,
                                                                groupName = c(groupNames(plistsMSMS)[2], groupNames(plistsMSMS)[3]),
                                                                analysis = c(analyses(plistsMSMS)[1], analyses(plistsMSMS)[1]),
                                                                MSLevel = 2))
    expect_doppel("mspl-spec_sim-avg-ms", function() plotSpectrum(plists,
                                                                  groupName = c(groupNames(plists)[32], groupName = groupNames(plists)[33]),
                                                                  MSLevel = 1))
    expect_doppel("mspl-spec_sim-avg-msms", function() plotSpectrum(plistsMSMS,
                                                                    groupName = c(groupNames(plistsMSMS)[14], groupNames(plistsMSMS)[15]),
                                                                    MSLevel = 2))
})

if (testWithSets())
{
    fgOneEmptySet <- makeOneEmptySetFGroups(fGroups)
    plistsOneEmptySet <- generateMSPeakLists(fgOneEmptySet, "mzr")
    plotSetSpecFG <- groupNames(setObjects(plistsMSMS)[[2]])[12]
}

test_that("sets functionality", {
    skip_if_not(testWithSets())
    
    expect_equal(analysisInfo(plists[, sets = "positive"], TRUE)[, 1:4], getTestAnaInfoPos(getTestAnaInfoAnn()),
                 check.attributes = FALSE)
    expect_equal(plists, plists[, sets = sets(plists)])
    expect_length(plists[, sets = character()], 0)
    expect_equal(sets(filter(plists, sets = "positive", negate = TRUE)), "negative")
    expect_setequal(groupNames(plists), unique(unlist(lapply(setObjects(plists), groupNames))))
    expect_setequal(groupNames(unset(plists, "positive")), groupNames(setObjects(plists)[[1]]))
    expect_setequal(groupNames(unset(plistsOneEmptySet, "positive")), groupNames(setObjects(plistsOneEmptySet)[[1]]))
    expect_length(unset(plistsOneEmptySet, "negative"), 0)
    
    expect_doppel("mspl-spec-set", function() plotSpectrum(plistsMSMS, groupName = plotSetSpecFG,
                                                           MSLevel = 2, perSet = FALSE))
    expect_doppel("mspl-spec-set-perset", function() plotSpectrum(plistsMSMS, groupName = plotSetSpecFG,
                                                                  MSLevel = 2, perSet = TRUE, mirror = FALSE))
    expect_doppel("mspl-spec-set-mirror", function() plotSpectrum(plistsMSMS, groupName = plotSetSpecFG,
                                                                  MSLevel = 2, perSet = TRUE, mirror = TRUE))
})
