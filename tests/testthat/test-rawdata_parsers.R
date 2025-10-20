# SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

context("raw data parsing")

anaInfo <- getTestAnaInfo()

anaInfoOne <- getTestAnaInfo()[4, ]

fList <- getTestFeatures(noiseThrInt = 1E5)

fgOpenMS <- groupFeatures(fList, "openms")

test_that("TICs", {
    checkmate::expect_data_table(getTICs(anaInfoOne))
    expect_equal(unique(getTICs(anaInfoOne, MSLevel = 2)[["MSLevel"]]), 2)
    expect_equal(colnames(getTICs(anaInfoOne)), c("analysis", "replicate", "ret", "MSLevel", "intensity"))
    expect_equal(nrow(getTICs(anaInfoOne)), 759)
    expect_lt(nrow(getTICs(anaInfoOne, retentionRange = c(10, 200))), 759)
    checkmate::expect_data_table(getTICs(fList, MSLevel = 2))
    checkmate::expect_data_table(getTICs(fgOpenMS))
})

test_that("BPCs", {
    checkmate::expect_data_table(getBPCs(anaInfoOne))
    expect_equal(colnames(getBPCs(anaInfoOne)), c("analysis", "replicate", "ret", "MSLevel", "intensity"))
    expect_equal(nrow(getBPCs(anaInfoOne)), 759)
    checkmate::expect_data_table(getBPCs(fList))
    checkmate::expect_data_table(getBPCs(fgOpenMS))
})

test_that("test plot TICs and BPCs", {
    expect_doppel("raw-tic", function() plotTICs(anaInfoOne))
    expect_doppel("raw-tic-col", function() plotTICs(fList, groupBy = "replicate", MSLevel = 1))
    expect_doppel("raw-ms2", function() plotTICs(fgOpenMS[2:4, ], groupBy = "replicate", MSLevel = 2))
    expect_doppel("raw-bpc", function() plotBPCs(anaInfo[2:3, ]))
    expect_doppel("raw-tic-feat", function() plotTICs(fList[1, ]))
    expect_doppel("raw-tic-fg", function() plotTICs(fgOpenMS[1, ]))
})

test_that("EICs", {
    anaInfoIMSOne <- as.data.table(getTestAnaInfoIMS()[4, ])
    susps <- as.data.table(patRoonDataIMS::suspectsPos)[name %in% c("Sulfamethoxazole", "Benzotriazole", "Metoprolol")]
    
    EICInfoList <- setNames(list(
        data.table::data.table(mzmin = susps$mz-0.01, mzmax = susps$mz+0.01, retmin = 200, retmax = 300, mobmin = 0.5, mobmax = 0.9)
    ), anaInfoIMSOne$analysis)
    EICInfoListOne <- EICInfoList
    EICInfoListOne[[1]] <- EICInfoListOne[[1]][1]
    
    # Test adjacency and intensity filters
    # NOTE: the thresholds were based on manual EIC inspection...
    eicsAdjP <- doGetEICs(anaInfoIMSOne, EICInfoListOne, gapFactor = 3, minEICAdjIntensity = 5E3,
                          minEICAdjPoints = 6, minEICAdjTime = 0, mode = "simple")
    eicsAdjT <- doGetEICs(anaInfoIMSOne, EICInfoListOne, gapFactor = 3, minEICAdjIntensity = 5E3,
                          minEICAdjPoints = 0, minEICAdjTime = 4, mode = "simple")
    eicsI <- doGetEICs(anaInfoIMSOne, EICInfoListOne, gapFactor = 3, minEICIntensity = 3E4, mode = "simple")
    eicsNoF <- doGetEICs(anaInfoIMSOne, EICInfoListOne, gapFactor = 3, mode = "simple")
    
    expect_length(pruneList(eicsAdjP[[1]], checkZeroRows = TRUE), 0)
    expect_length(pruneList(eicsAdjT[[1]], checkZeroRows = TRUE), 0)
    expect_length(pruneList(eicsI[[1]], checkZeroRows = TRUE), 0)
    expect_gt(length(pruneList(eicsNoF[[1]])), 0)
    
    # Test smoothing and frame summing affects mz or mobility results
    eicsSmoothed <- doGetEICs(anaInfoIMSOne, EICInfoListOne, gapFactor = 3, mode = "full", sumFramesMZ = 7,
                              smoothWindowMZ = 3, smoothExtMZ = 0.02, sumFramesMob = 7, smoothWindowMob = 7,
                              smoothExtMob = 0.1)
    eicsUnsmoothed <- doGetEICs(anaInfoIMSOne, EICInfoList[1], gapFactor = 3, mode = "full")

    expect_true(any(eicsSmoothed[[1]][[1]][, "mz"] != eicsUnsmoothed[[1]][[1]][, "mz"]))
    expect_true(any(eicsSmoothed[[1]][[1]][, "mzBP"] != eicsUnsmoothed[[1]][[1]][, "mzBP"]))
    expect_true(any(eicsSmoothed[[1]][[1]][, "mobility"] != eicsUnsmoothed[[1]][[1]][, "mobility"]))
    expect_true(any(eicsSmoothed[[1]][[1]][, "mobilityBP"] != eicsUnsmoothed[[1]][[1]][, "mobilityBP"]))
    
    # Test topMost filter
    eicsTop1 <- doGetEICs(anaInfoIMSOne, EICInfoList, gapFactor = 3, topMost = 1, mode = "full")
    expect_length(pruneList(eicsTop1[[1]], checkZeroRows = TRUE), 1)
    
    # Test pad=TRUE
    eicsPadded <- doGetEICs(anaInfoIMSOne, EICInfoListOne, gapFactor = 3, pad = TRUE, mode = "simple")
    eicsUnpadded <- doGetEICs(anaInfoIMSOne, EICInfoListOne, gapFactor = 3, pad = FALSE, mode = "simple")
    expect_true(nrow(eicsPadded[[1]][[1]]) >= nrow(eicsUnpadded[[1]][[1]]))
    expect_false(any(eicsUnpadded[[1]][[1]][, "intensity"] == 0))
    expect_true(any(eicsPadded[[1]][[1]][, "intensity"] == 0))
})

test_that("EIMs", {
    anaInfoIMSOne <- as.data.table(getTestAnaInfoIMS()[4, ])
    susps <- as.data.table(patRoonDataIMS::suspectsPos)[name == "Sulfamethoxazole"]
    
    EIMInfoList <- setNames(list(
        data.table::data.table(mzmin = susps$mz-0.01, mzmax = susps$mz+0.01, retmin = susps$rt-6, retmax = susps$rt+6,
                               mobmin = 0, mobmax = 2)
    ), anaInfoIMSOne$analysis)
    
    # Test compress=TRUE
    eimsCompressed <- doGetEIMs(anaInfoIMSOne, EIMInfoList, minIntensity = 0, smooth = "none", smLength = 0,
                                sgOrder = 0, compress = TRUE)
    eimsUncompressed <- doGetEIMs(anaInfoIMSOne, EIMInfoList, minIntensity = 0, smooth = "none", smLength = 0,
                                  sgOrder = 0, compress = FALSE)
    expect_lte(nrow(eimsCompressed[[1]][[1]]), nrow(eimsUncompressed[[1]][[1]]))
    expect_gt(sum(eimsUncompressed[[1]][[1]][, "intensity"] == 0), sum(eimsCompressed[[1]][[1]][, "intensity"] == 0))
})
