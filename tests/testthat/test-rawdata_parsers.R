context("raw data parsing")

anaInfo <- getTestAnaInfo()

anaInfoOne <- getTestAnaInfo()[4, ]

fList <- getTestFeatures(noiseThrInt = 1E5)

fgOpenMS <- groupFeatures(fList, "openms")

test_that("TICs", {
    expect_s3_class(getTICs(anaInfoOne), "data.table")
    expect_equal(unique(getTICs(anaInfoOne, MSLevel = 2)[["MSLevel"]]), 2)
    expect_equal(colnames(getTICs(anaInfoOne)), c("analysis", "group", "ret", "MSLevel", "intensity"))
    expect_equal(nrow(getTICs(anaInfoOne)), 759)
    expect_lt(nrow(getTICs(anaInfoOne, retentionRange = c(10, 200))), 759)
    expect_s3_class(getTICs(fList, MSLevel = 2), "data.table")
    expect_s3_class(getTICs(fgOpenMS), "data.table")
})

test_that("BPCs", {
    expect_s3_class(getBPCs(anaInfoOne), "data.table")
    expect_equal(colnames(getBPCs(anaInfoOne)), c("analysis", "group", "ret", "MSLevel", "mz", "intensity"))
    expect_equal(nrow(getBPCs(anaInfoOne)), 759)
    expect_s3_class(getBPCs(fList), "data.table")
    expect_s3_class(getBPCs(fgOpenMS), "data.table")
})

test_that("test plot TICs and BPCs", {
    expect_plot(plotTICs(anaInfoOne))
    expect_plot(plotTICs(fList, colourBy = "rGroups", MSLevel = 1))
    expect_plot(plotTICs(fgOpenMS[2:4, ], colourBy = "rGroups", MSLevel = 2))
    expect_plot(plotBPCs(anaInfo[2:3, ]))
    expect_plot(plotTICs(fList[1, ]))
    expect_plot(plotTICs(fgOpenMS[1, ]))
})
