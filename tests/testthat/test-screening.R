context("screening")

fGroups <- getTestFGroups(getTestAnaInfo())
scr <- screenTargets(fGroups, patRoonData::targets)
scrF <- screenTargets(getFeatures(fGroups), patRoonData::targets)

test_that("target screening is OK", {
    expect_equal(nrow(scr), nrow(patRoonData::targets))
    expect_length(unique(scrF$name), nrow(patRoonData::targets))
    expect_known_value(scr, testFile("screening"))
    expect_known_value(scrF[, -"feature"], testFile("screening-feat"))
    expect_length(groupFeaturesScreening(fGroups, scr), nrow(patRoonData::targets))
})
