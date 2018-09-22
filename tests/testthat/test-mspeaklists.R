context("MS peak lists")

fGroups <- getTestFGroups()[4:5, 1:100]
plists <- generateMSPeakLists(fGroups, "mzr")

test_that("verify generation of MS peak lists", {
    # expect_known_hash(plists, "d4b026d4e3") # use a hash here because of the large resulting file
    expect_known_value(plists, testFile("plists-mzr"))
})

test_that("verify show output", {
    expect_known_show(plists, testFile("plists-mzr", text = TRUE))
})

checkMinInt <- function(plists, relative, doMSMS)
{
    pl <- peakLists(plists)
    return(min(sapply(pl, function(ana) min(sapply(ana, function(grp)
    {
        if (!doMSMS)
            ints <- grp$MS$intensity
        else if (!is.null(grp$MSMS))
            ints <- grp$MSMS$intensity
        else
            return(NA)

        if (length(ints) == 0)
            return(NA)
        
        ret <- min(ints)
        if (relative)
            ret <- ret / max(ints)
        return(ret)
    }), na.rm = TRUE)), na.rm = TRUE))
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
})

test_that("feature group filtering", {
    expect_equal(filterBy(plists, fGroups), fGroups)
    expect_lt(length(filterBy(plists, fGroups, onlyMSMS = TRUE)), length(fGroups))
    expect_length(filterBy(plistsEmpty, fGroups), 0)
})

test_that("basic subsetting", {
    expect_length(plists["nope"], 0)
    expect_equivalent(analyses(plists[1:2]), analyses(fGroups)[1:2])
    expect_equivalent(analyses(plists[analyses(fGroups)[1:2]]), analyses(fGroups)[1:2])
    expect_equivalent(analyses(plists[c(FALSE, TRUE)]), analyses(fGroups)[2])
    expect_equivalent(groupNames(plists[, 1:2]), groupNames(plists)[1:2])
    expect_equivalent(groupNames(plists[, groupNames(plists)[2:3]]), groupNames(plists)[2:3])
    expect_equivalent(groupNames(plists[, c(FALSE, TRUE)]), groupNames(plists)[c(FALSE, TRUE)])
    expect_equal(length(plists[FALSE]), 0)
    expect_length(plistsEmpty[1:5], 0)
})
