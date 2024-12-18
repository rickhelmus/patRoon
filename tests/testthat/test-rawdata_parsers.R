# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
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
    expect_equal(colnames(getTICs(anaInfoOne)), c("analysis", "group", "ret", "MSLevel", "intensity"))
    expect_equal(nrow(getTICs(anaInfoOne)), 759)
    expect_lt(nrow(getTICs(anaInfoOne, retentionRange = c(10, 200))), 759)
    checkmate::expect_data_table(getTICs(fList, MSLevel = 2))
    checkmate::expect_data_table(getTICs(fgOpenMS))
})

test_that("BPCs", {
    checkmate::expect_data_table(getBPCs(anaInfoOne))
    expect_equal(colnames(getBPCs(anaInfoOne)), c("analysis", "group", "ret", "MSLevel", "mz", "intensity"))
    expect_equal(nrow(getBPCs(anaInfoOne)), 759)
    checkmate::expect_data_table(getBPCs(fList))
    checkmate::expect_data_table(getBPCs(fgOpenMS))
})

test_that("test plot TICs and BPCs", {
    expect_doppel("raw-tic", function() plotTICs(anaInfoOne))
    expect_doppel("raw-tic-col", function() plotTICs(fList, colourBy = "rGroups", MSLevel = 1))
    expect_doppel("raw-ms2", function() plotTICs(fgOpenMS[2:4, ], colourBy = "rGroups", MSLevel = 2))
    expect_doppel("raw-bpc", function() plotBPCs(anaInfo[2:3, ]))
    expect_doppel("raw-tic-feat", function() plotTICs(fList[1, ]))
    expect_doppel("raw-tic-fg", function() plotTICs(fgOpenMS[1, ]))
})
