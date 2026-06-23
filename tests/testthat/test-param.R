# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

# UNDONE: just some dummy tests now

test_that("all param classes can be constructed without errors", {
    # Feature params
    expect_error(FeaturesOpenMSParam(), NA)
    expect_error(FeaturesXCMS3Param(), NA)
    expect_error(FeaturesEnviPickParam(), NA)
    expect_error(FeaturesKPIC2Param(), NA)
    expect_error(FeaturesPiekParam(), NA)
    expect_error(FeaturesSAFDParam(), NA)

    # Component params
    expect_error(ComponentsOpenMSParam(), NA)
    expect_error(ComponentsIntClustParam(), NA)
    expect_error(ComponentsRAMClustRParam(), NA)
    expect_error(ComponentsNontargetParam(), NA)
    expect_error(ComponentsSpecClustParam(), NA)
    expect_error(ComponentsCliqueMSParam(), NA)
    expect_error(ComponentsCAMERAParam(), NA)
    expect_error(ComponentsTPsParam(), NA)

    # Feature group params
    expect_error(FeatureGroupsOpenMSParam(), NA)
    expect_error(FeatureGroupsXCMS3Param(), NA)
    expect_error(FeatureGroupsKPIC2Param(), NA)
    expect_error(FeatureGroupsGreedyParam(), NA)

    # Other
    expect_error(MSPeakListsParam(), NA)
})
