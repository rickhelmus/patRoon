# SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

context("MS peak lists")

fGroups <- getTestFGroups(getTestAnaInfoAnn())[, 1:100]
ovFGroup <- names(overlap(fGroups, which = c("positive", "negative"), aggregate = "set"))[1]
plists <- generateMSPeakLists(fGroups)
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
    plistsDA <- generateMSPeakListsDA(fgDA, save = FALSE, bgsubtr = FALSE)
    plistsDAEmpty <- generateMSPeakListsDA(fgDA[, "nope"], save = FALSE, bgsubtr = FALSE)

    fgDA2 <- getTestFGroupsDA(getDAAnaInfo("std1"))[, 1:25]
    plistsDAFMF <- generateMSPeakListsDAFMF(fgDA2)
}

# remove some inconsistent metadata
plistsNoIM <- plists
rmCols <- c("electronBeamEnergy")
clearMD <- function(md) lapply(md, function(mda) lapply(mda, function(mdf) lapply(mdf, function(mds) mds[, setdiff(names(mds), rmCols), with = FALSE])))
plistsNoIM@metadata <- clearMD(plistsNoIM@metadata)
plistsNoIM@setObjects <- lapply(setObjects(plistsNoIM), function(so) { so@metadata <- clearMD(so@metadata); so })
plistsNoIM@analysisInfo <- data.table() # remove as it is system dependent

test_that("verify generation of MS peak lists", {
    expect_known_value(plistsNoIM, testFile("plists"))

    skip_if_not(doDATests())
    expect_known_value(plistsDA, testFile("plists-DA"))
    expect_known_value(plistsDAFMF, testFile("plists-DAFMF"))
})

test_that("verify show output", {
    expect_known_show(plists, testFile("plists", text = TRUE))

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

testMSPLAbundance <- function(obj, avg, rel, thr)
{
    # NOTE: exclude data with relative abundance==1: these may be averaged from less than available spectra
    colRel <- if (!avg) "abundance_rel" else "fgroup_abundance_rel"
    tab <- as.data.table(obj, averaged = avg)[get(colRel) != 1]
    col <- sprintf("%sabundance_%s", if (!avg) "" else "fgroup_", if (rel) "rel" else "abs")
    expect_gte(min(tab[[col]]), thr, label = col)
}

getMSPLProp <- function(obj, avg, col)
{
    tab <- as.data.table(obj, averaged = avg)
    return(tab[precursor == FALSE][[col]]) # ignore precursors: these are usually not filtered away
}

testMaxMZOverPrec <- function(obj, thr, neg)
{
    tab <- as.data.table(obj)
    tab[, precursor_mz := if (any(precursor)) mz[precursor == TRUE] else NA_real_, by = c("group", "set", "type")]
    tab <- tab[!is.na(precursor_mz) & precursor == FALSE] # ignore precursors: these are usually not filtered away
    if (!neg)
        expect_true(all((tab$mz - tab$precursor_mz) <= thr))
    else
        expect_true(all((tab$mz - tab$precursor_mz) > thr))
}

test_that("avg params", {
    expect_gte(min(as.data.table(generateMSPeakLists(fGroups, avgFeatParams = getDefAvgPListParams(minIntensityPost = 2500)), averaged = FALSE)$intensity), 2500)
    expect_gte(min(as.data.table(generateMSPeakLists(fGroups, avgFGroupParams = getDefAvgPListParams(minIntensityPost = 2500)))$intensity), 2500)
    testMSPLAbundance(generateMSPeakLists(fGroups, avgFeatParams = getDefAvgPListParams(minAbundanceRel = 0.5)), FALSE, TRUE, 0.5)
    testMSPLAbundance(generateMSPeakLists(fGroups, avgFeatParams = getDefAvgPListParams(minAbundanceAbs = 10)), FALSE, FALSE, 10)
    # UNDONE: these tests are not so useful as there are only two analyses to average... so the abundance will also be between one (in case of feature being in only one analyses) or two
    # testMSPLAbundance(generateMSPeakLists(fGroups, avgFGroupParams = getDefAvgPListParams(minAbundanceRel = 0.5)), TRUE, TRUE, 0.5)
    # testMSPLAbundance(generateMSPeakLists(fGroups, avgFGroupParams = getDefAvgPListParams(minAbundanceAbs = 2)), TRUE, FALSE, 2)
})

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
    expect_equal(delete(plists, k = 1, reAverage = FALSE)[[ovFGroup]], plists[[ovFGroup]])
    expect_false(isTRUE(all.equal(delete(plists, k = 1, reAverage = TRUE)[[ovFGroup]], plists[[ovFGroup]])))
    expect_length(delete(plists), 0)

    # we only test MSLevel selection here: the mechanism is the same for all filters
    expect_gte(checkIntLimit(filter(plists, MSLevel = 1, absMinIntensity = 2500), FALSE, TRUE, FALSE), 2500)
    expect_gte(checkIntLimit(filter(plists, MSLevel = 2, absMinIntensity = 2500), FALSE, TRUE, TRUE), 2500)
    expect_lte(checkIntLimit(filter(plists, MSLevel = 1, absMinIntensity = 2500, negate = TRUE), FALSE, FALSE, FALSE), 2500)
    expect_lte(checkIntLimit(filter(plists, MSLevel = 2, absMinIntensity = 2500, negate = TRUE), FALSE, FALSE, TRUE), 2500)

    expect_gte(checkIntLimit(filter(plists, MSLevel = 1, relMinIntensity = 0.2), TRUE, TRUE, FALSE, plists), 0.2)
    expect_gte(checkIntLimit(filter(plists, MSLevel = 2, relMinIntensity = 0.2), TRUE, TRUE, TRUE, plists), 0.2)
    expect_lte(checkIntLimit(filter(plists, MSLevel = 1, relMinIntensity = 0.2, negate = TRUE), TRUE, FALSE, FALSE, plists), 0.2)
    expect_lte(checkIntLimit(filter(plists, MSLevel = 2, relMinIntensity = 0.2, negate = TRUE), TRUE, FALSE, TRUE, plists), 0.2)
    
    expect_lte(checkPeaksLimit(filter(plists, MSLevel = 1, topMostPeaks = 5), FALSE, FALSE), 5)
    expect_lte(checkPeaksLimit(filter(plists, MSLevel = 2, topMostPeaks = 5), FALSE, TRUE), 5)
    expect_lte(checkPeaksLimit(filter(plists, MSLevel = 1, topMostPeaks = 5, negate = TRUE), FALSE, FALSE), 5)
    expect_lte(checkPeaksLimit(filter(plists, MSLevel = 2, topMostPeaks = 5, negate = TRUE), FALSE, TRUE), 5)

    expect_true(all(sapply(averagedPeakLists(plistsMSMS), function(pl) !is.null(pl[["MSMS"]]))))
    expect_true(all(sapply(averagedPeakLists(filter(plists, withMSMS = TRUE, negate = TRUE)),
                           function(pl) is.null(pl[["MSMS"]]))))

    expect_min_gte(getMSPLProp(filter(plists, minAbundanceFeatRel = 0.5), "abundance_rel", avg = FALSE), 0.5)
    expect_min_gte(getMSPLProp(filter(plists, minAbundanceFeatAbs = 10), "abundance_abs", avg = FALSE), 10)
    expect_max_lt(getMSPLProp(filter(plists, minAbundanceFeatRel = 0.5, negate = TRUE), "abundance_rel", avg = FALSE), 0.5)
    expect_max_lt(getMSPLProp(filter(plists, minAbundanceFeatAbs = 10, negate = TRUE), "abundance_abs", avg = FALSE), 10)
    
    expect_min_gte(getMSPLProp(filter(plists, minAbundanceFeatRel = 0.5), "feat_abundance_rel", avg = TRUE), 0.5)
    expect_min_gte(getMSPLProp(filter(plists, minAbundanceFeatAbs = 10), "feat_abundance_abs", avg = TRUE), 10)
    expect_max_lt(getMSPLProp(filter(plists, minAbundanceFeatRel = 0.5, negate = TRUE), "feat_abundance_rel", avg = TRUE), 0.5)
    expect_max_lt(getMSPLProp(filter(plists, minAbundanceFeatAbs = 10, negate = TRUE), "feat_abundance_abs", avg = TRUE), 10)
    
    expect_min_gte(getMSPLProp(filter(plists, minAbundanceFGroupRel = 0.5), "fgroup_abundance_rel", avg = TRUE), 0.5)
    expect_min_gte(getMSPLProp(filter(plists, minAbundanceFGroupAbs = 2), "fgroup_abundance_abs", avg = TRUE), 2)
    expect_max_lt(getMSPLProp(filter(plists, minAbundanceFGroupRel = 1, negate = TRUE), "fgroup_abundance_rel", avg = TRUE), 1)
    expect_max_lt(getMSPLProp(filter(plists, minAbundanceFGroupAbs = 2, negate = TRUE), "fgroup_abundance_abs", avg = TRUE), 2)

    testMaxMZOverPrec(filter(plists, maxMZOverPrec = 0.5), 0.5, FALSE)
    testMaxMZOverPrec(filter(plists, maxMZOverPrec = 0.5, negate = TRUE), 0.5, TRUE)
    
    rmMZ <- plists[[1]][[1]]$mz[1]
    expect_true(all(abs(getMSPLProp(filter(plists, removeMZs = list(rmMZ, rmMZ), mzWindow = 0.005), TRUE, "mz") - rmMZ) > 0.005))
    expect_true(all(abs(getMSPLProp(filter(plists, removeMZs = list(rmMZ, rmMZ), mzWindow = 0.005, negate = TRUE), TRUE, "mz") - rmMZ) <= 0.005))
    
    expect_lte(length(filter(plists, MSLevel = 2, topMostPeaks = 5, retainPrecursor = FALSE)),
               length(filter(plists, MSLevel = 2, topMostPeaks = 5)))

    # isotopes for MS peak lists should not exceed M+5 (default)
    expect_lte(max(as.data.table(filter(
        plists, MSLevel = 1, isolatePrec = getDefIsolatePrecParams()))[type == "MS", diff(range(mz)), by = c("group", "set")][["V1"]]), 5)
    # half with z=2
    expect_lte(max(as.data.table(filter(
        plists, MSLevel = 1, isolatePrec = getDefIsolatePrecParams(z=2)))[type == "MS", diff(range(mz)), by = c("group", "set")][["V1"]]), 2.5)

    # UNDONE: deisotope?
})

plistsEmpty <- removePrecursors(filter(plists, absMinIntensity = 1E9))
plistsEmptyMS <- removePrecursors(filter(plists, MSLevel = 1, absMinIntensity = 1E9))

test_that("empty object", {
    expect_length(plistsEmpty, 0)
    expect_length(delete(plistsEmpty), 0)
    expect_length(filter(plistsEmpty, relMinIntensity = 0.2, topMostPeaks = 10), 0)
    expect_length(generateMSPeakLists(getEmptyTestFGroups()), 0)
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
    expect_equivalent(plists[1:2, reAverage = FALSE][[ovFGroup]], plists[[ovFGroup]])
    expect_false(isTRUE(all.equal(plists[1:2, reAverage = TRUE][[ovFGroup]], plists[[ovFGroup]])))
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
                analysis2 = rep(analyses(plists)[1], length(groupNames(plists))), MSLevel = 1, expectNA = TRUE)
    testSpecSim(plists, groupNames(plists), groupNames(plists), MSLevel = 1, expectNA = TRUE)
    testSpecSim(plists, groupNames(plists), groupNames(plists), MSLevel = 2, expectNA = TRUE)
    testSpecSim(plists, groupNames(plists)[1:5], groupNames(plists)[6:10], MSLevel = 1)
    
    testSpecSim(plists, groupNames(plists), NULL, MSLevel = 1,
                analysis1 = rep(analyses(plists)[1], length(groupNames(plists))), expectNA = TRUE)
    testSpecSim(plists, groupNames(plists), NULL, MSLevel = 1, expectNA = TRUE)
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

fgOneEmptySet <- makeOneEmptySetFGroups(fGroups)
plistsOneEmptySet <- generateMSPeakLists(fgOneEmptySet)
plotSetSpecFG <- groupNames(setObjects(plistsMSMS)[[2]])[12]

test_that("sets functionality", {
    expect_equal(analysisInfo(plists[, sets = "positive"])[, -"set"], getTestAnaInfoPos(getTestAnaInfoAnn()),
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

test_that("IMS tests", {
    fGroupsIMS <- getTestFGroupsIMS()[, 1:50]
    fGroupsIMS <- doAssignMobs(fGroupsIMS)
    plistsIMS <- generateMSPeakLists(fGroupsIMS)
    specSims <- spectrumSimilarityMobility(plistsIMS, fGroupsIMS, doFGroups = TRUE)
    checkmate::expect_data_table(specSims)
    checkmate::expect_names(names(specSims), identical.to = c("group", "ims_parent_group", "similarity"))
    expect_range(specSims$similarity, c(0, 1))
    expect_setequal(c(unique(specSims$group), unique(specSims$ims_parent_group)), groupNames(plistsIMS))
    
    specSimsFeats <- spectrumSimilarityMobility(plistsIMS, fGroupsIMS, doFGroups = FALSE)
    checkmate::expect_data_table(specSimsFeats)
    checkmate::expect_names(names(specSimsFeats), identical.to = c("group", "ims_parent_group", "analysis", "similarity"))
    expect_range(specSimsFeats$similarity, c(0, 1))
    expect_true(anyNA(specSimsFeats$similarity)) # since this is sets data the IMS features will only be present in one set
    expect_setequal(c(unique(specSimsFeats$group), unique(specSimsFeats$ims_parent_group)), groupNames(plistsIMS))
    
    expect_null(spectrumSimilarityMobility(delete(plistsIMS), fGroupsIMS))
    expect_warning(spectrumSimilarityMobility(plistsIMS, delete(fGroupsIMS)), "No relevant")
    expect_warning(spectrumSimilarityMobility(plistsIMS, delete(fGroupsIMS), warn = FALSE), NA)
    expect_error(spectrumSimilarityMobility(plistsIMS, fGroupsIMS[, IMS = FALSE]), "No mobility")
    expect_warning(spectrumSimilarityMobility(plistsIMS, fGroupsIMS[, IMS = TRUE]), "No relevant")
})

doGetBG <- function(ai = getTestAnaInfoNS()[getTestAnaInfoNS()$replicate == "solvent-pos", ], minBPIntensity = 1E5, ...)
{
    getBGMSMSPeaks(ai, minBPIntensity = minBPIntensity, ...)
}
bg <- doGetBG()
test_that("BG MSMS", {
    expect_known_value(bg, testFile("bg-msms"))
    checkmate::expect_data_table(bg, any.missing = FALSE)
    checkmate::expect_names(names(bg), must.include = c("mz", "intensity", "abundance_rel_ana", "abundance_abs_ana", 
                                                        "abundance_rel_spec", "abundance_abs_spec"))
    expect_equal(bg, doGetBG(ai = getTestAnaInfoNS(), replicates = "solvent-pos"))
    expect_lt(sum(doGetBG(minBPIntensity = 3E5)$abundance_abs_spec), sum(bg$abundance_abs_spec))
    expect_min_gte(doGetBG(avgSpectraParams = getDefAvgPListParams(minAbundanceRel = 0.5, topMost = 25))$abundance_rel_spec, 0.5)
    expect_min_gte(doGetBG(avgSpectraParams = getDefAvgPListParams(minAbundanceAbs = 50, topMost = 25))$abundance_abs_spec, 50)
    expect_min_gte(doGetBG(avgAnalysesParams = getDefAvgPListParams(minAbundanceRel = 0.5, topMost = 25))$abundance_rel_ana, 0.5)
    expect_min_gte(doGetBG(avgAnalysesParams = getDefAvgPListParams(minAbundanceAbs = 2, topMost = 25))$abundance_abs_ana, 2)
})
