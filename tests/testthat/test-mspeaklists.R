context("MS peak lists")

fGroups <- getTestFGroups()
plists <- generateMSPeakLists(fGroups, "mzr")

test_that("verify generation of MS peak lists", {
    expect_known_hash(plists, "d4b026d4e3") # use a hash here because of the large resulting file
})

test_that("verify show output", {
    expect_known_output(show(plists), testFile("plists-mzr", text = TRUE))
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