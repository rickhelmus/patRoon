context("MS peak lists")

fGroups <- getTestFGroupsAnn()[, 1:100]
plists <- generateMSPeakLists(fGroups, "mzr")
plistsMSMS <- filter(plists, withMSMS = TRUE)

if (doDATests())
{
    # NOTE: use different analyses than first: this call will clearout all DA
    # compounds, forcing other calls (from features/formulas tests) to re-run
    # FMF
    fgDA <- groupFeatures(findFeatures(getDAAnaInfo()[2, ], "bruker"), "openms")

    # NOTE: set bgsubtr to FALSE: subtraction might remove pecursor peaks of
    # (often wrong) low intensity features and result in warnings during
    # averaging
    plistsDA <- generateMSPeakLists(fgDA, "bruker", save = FALSE, bgsubtr = FALSE)
    plistsDAEmpty <- generateMSPeakLists(fgDA["nope"], "bruker", save = FALSE, bgsubtr = FALSE)

    fgDA2 <- groupFeatures(findFeatures(getDAAnaInfo()[1, ], "bruker"), "openms")
    plistsDAFMF <- generateMSPeakLists(fgDA2, "brukerfmf")
}

# remove Ion Mobility data as it inconsistently is present or not
plistsNoIM <- plists
plistsNoIM@metadata <- lapply(plistsNoIM@metadata, function(mda) lapply(mda, function(mdf) lapply(mdf, function(mds) mds[, setdiff(names(mds), "ionMobilityDriftTime"), with = FALSE])))

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
        if (!doMSMS && !is.null(pl[["MS"]]))
        {
            ints <- pl[["MS"]]$intensity
            if (!is.null(plorig))
                intsorig <- plorig[["MS"]]$intensity
        }
        else if (doMSMS && !is.null(pl$MSMS))
        {
            ints <- pl[["MSMS"]]$intensity
            if (!is.null(plorig))
                intsorig <- plorig[["MSMS"]]$intensity
        }
        else
            return(NA)

        if (length(ints) == 0)
            return(NA)

        ret <- lim(ints)
        if (relative)
            ret <- ret / max(intsorig)

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

test_that("filtering", {
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
        plists, isolatePrec = getDefIsolatePrecParams()))[type == "MS", diff(range(mz)), by = "group"][[2]]), 5)
    # half with z=2
    expect_lte(max(as.data.table(filter(
        plists, isolatePrec = getDefIsolatePrecParams(z=2)))[type == "MS", diff(range(mz)), by = "group"][[2]]), 2.5)

    # UNDONE: deisotope?
})

plistsEmpty <- removePrecursors(filter(plists, absMSIntThr = 1E9, absMSMSIntThr = 1E9))
plistsEmptyMS <- removePrecursors(filter(plists, absMSIntThr = 1E9))

test_that("empty object", {
    expect_length(plistsEmpty, 0)
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
    expect_equivalent(analyses(plists[c(FALSE, TRUE)]), analyses(fGroups)[2])
    expect_equivalent(groupNames(plists[, 1:2]), groupNames(plists)[1:2])
    expect_equivalent(groupNames(plists[, groupNames(plists)[2:3]]), groupNames(plists)[2:3])
    expect_equivalent(groupNames(plists[, c(FALSE, TRUE)]), groupNames(plists)[c(FALSE, TRUE)])
    expect_equal(length(plists[FALSE]), 0)
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

test_that("plotting works", {
    expect_doppel("mspl-spec-ms", function() plotSpectrum(plists, groupName = groupNames(plists)[70],
                                                          analysis = analyses(plists)[1], MSLevel = 1))
    expect_doppel("mspl-spec-msms", function() plotSpectrum(plistsMSMS, groupName = groupNames(plistsMSMS)[2],
                                                            analysis = analyses(plistsMSMS)[1], MSLevel = 2))
    expect_doppel("mspl-spec-avg-ms", function() plotSpectrum(plists, groupName = groupNames(plists)[70],
                                                              MSLevel = 1))
    expect_doppel("mspl-spec-avg-msms", function() plotSpectrum(plistsMSMS, groupName = groupNames(plistsMSMS)[2],
                                                                MSLevel = 2))

    expect_ggplot(plotSpectrum(plists, groupName = groupNames(plists)[70],
                               analysis = analyses(plists)[1], MSLevel = 1,
                               useGGPlot2 = TRUE))
    expect_ggplot(plotSpectrum(plistsMSMS, groupName = groupNames(plistsMSMS)[2],
                               analysis = analyses(plistsMSMS)[1], MSLevel = 2,
                               useGGPlot2 = TRUE))
    expect_ggplot(plotSpectrum(plists, groupName = groupNames(plists)[70], MSLevel = 1,
                               useGGPlot2 = TRUE))
    expect_ggplot(plotSpectrum(plistsMSMS, groupName = groupNames(plistsMSMS)[2], MSLevel = 2,
                               useGGPlot2 = TRUE))
})
