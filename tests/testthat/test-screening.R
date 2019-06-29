context("screening")

fGroups <- getTestFGroups(getTestAnaInfo())
scr <- screenSuspects(fGroups, patRoonData::targets)
scrF <- screenSuspects(getFeatures(fGroups), patRoonData::targets)

test_that("target screening is OK", {
    expect_equal(nrow(scr), nrow(patRoonData::targets))
    expect_length(unique(scrF$name), nrow(patRoonData::targets))
    expect_known_value(scr, testFile("screening"))
    expect_known_value(scrF[, -"feature"], testFile("screening-feat"))
    expect_length(groupFeaturesScreening(fGroups, scr), nrow(patRoonData::targets))

    # check suspects without retention
    expect_gte(nrow(screenSuspects(fGroups, patRoonData::targets[, -3])), nrow(scr))
    expect_gte(nrow(screenSuspects(getFeatures(fGroups), patRoonData::targets[, -3])), nrow(scrF))
})
