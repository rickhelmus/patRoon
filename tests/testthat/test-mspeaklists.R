context("MS peak lists")

fGroups <- getTestFGroups()[4:5, 1:100]
plists <- generateMSPeakLists(fGroups, "mzr")

if (doDATests())
{
    # NOTE: use different analyses than first: this call will clearout all DA
    # compounds, forcing other calls (from features/formulas tests) to re-run
    # FMF
    fgDA <- groupFeatures(findFeatures(getDAAnaInfo()[2, ], "bruker"), "openms")
    plistsDA <- generateMSPeakLists(fgDA, "bruker", save = FALSE)
    plistsDAEmpty <- generateMSPeakLists(fgDA["nope"], "bruker", save = FALSE)
    
    fgDA2 <- groupFeatures(findFeatures(getDAAnaInfo()[1, ], "bruker"), "openms")
    plistsDAFMF <- generateMSPeakLists(fgDA2, "brukerfmf")
}

test_that("verify generation of MS peak lists", {
    expect_known_value(plists, testFile("plists-mzr"))
    
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

checkMinInt <- function(plists, relative, doMSMS)
{
    minIntPL <- function(pl)
    {
        if (!doMSMS)
            ints <- pl$MS$intensity
        else if (!is.null(pl$MSMS))
            ints <- pl$MSMS$intensity
        else
            return(NA)
        
        if (length(ints) == 0)
            return(NA)
        
        ret <- min(ints)
        if (relative)
            ret <- ret / max(ints)
        return(ret)
    }
    
    ftMin <- min(sapply(peakLists(plists), function(ana) min(sapply(ana, minIntPL), na.rm = TRUE)), na.rm = TRUE)
    fgMin <- min(sapply(averagedPeakLists(plists), minIntPL), na.rm = TRUE)
    return(min(ftMin, fgMin))
}

checkMaxPeaks <- function(plists, doMSMS)
{
    pl <- peakLists(plists)
    return(max(sapply(pl, function(ana) max(sapply(ana, function(grp)
    {
        if (!doMSMS)
            return(nrow(grp$MS))
        else if (!is.null(grp$MSMS))
            return(nrow(grp$MSMS))
        else
            return(NA)
    }), na.rm = TRUE)), na.rm = TRUE))
}

test_that("filtering", {
    expect_gte(checkMinInt(filter(plists, absMSIntThr = 500), FALSE, FALSE), 500)
    expect_gte(checkMinInt(filter(plists, absMSMSIntThr = 500), FALSE, TRUE), 500)
    expect_gte(checkMinInt(filter(plists, relMSIntThr = 0.2), TRUE, FALSE), 0.2)
    expect_gte(checkMinInt(filter(plists, relMSMSIntThr = 0.2), TRUE, TRUE), 0.2)
    expect_lte(checkMaxPeaks(filter(plists, topMSPeaks = 10), FALSE), 10)
    expect_lte(checkMaxPeaks(filter(plists, topMSMSPeaks = 10), TRUE), 10)
    # UNDONE: deisotope?
})

plistsEmpty <- filter(plists, absMSIntThr = 1E9, absMSMSIntThr = 1E9)
plistsEmptyMS <- filter(plists, absMSIntThr = 1E9)

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
})

test_that("plotting works", {
    expect_doppel("mspl-spec-ms", function() plotSpec(plists, groupName = groupNames(plists)[70],
                                                      analysis = analyses(plists)[1], MSLevel = 1))
    expect_doppel("mspl-spec-msms", function() plotSpec(plists, groupName = groupNames(plists)[70],
                                                        analysis = analyses(plists)[1], MSLevel = 2))
    expect_doppel("mspl-spec-avg-ms", function() plotSpec(plists, groupName = groupNames(plists)[70],
                                                          MSLevel = 1))
    expect_doppel("mspl-spec-avg-msms", function() plotSpec(plists, groupName = groupNames(plists)[70],
                                                            MSLevel = 2))

    expect_plot(print(plotSpec(plists, groupName = groupNames(plists)[70],
                               analysis = analyses(plists)[1], MSLevel = 1)))
    expect_plot(print(plotSpec(plists, groupName = groupNames(plists)[70],
                               analysis = analyses(plists)[1], MSLevel = 2)))
    expect_plot(print(plotSpec(plists, groupName = groupNames(plists)[70], MSLevel = 1)))
    expect_plot(print(plotSpec(plists, groupName = groupNames(plists)[70], MSLevel = 2)))
})
