context("components")

# set localMZRange=0 to keep isotopes
fGroups <- getTestFGroups(getTestAnaInfo()[4:5, ], localMZRange = 0)
# reduced set for CAMERA/RAMClustR; for nontarget we keep all to get less common homologues
fGroupsSimple <- fGroups[, 1:50]


# fix seed for reproducible clustering, suppress warnings about <5 samples
withr::with_seed(20, suppressWarnings(compsRC <- generateComponents(fGroupsSimple, "ramclustr", ionization = "positive")))
# UNDONE: getting unknown NaN warnings here...
suppressWarnings(compsCAM <- generateComponents(fGroupsSimple, "camera", ionization = "positive"))
compsNT <- generateComponents(fGroups, "nontarget", ionization = "positive")
compsInt <- generateComponents(fGroupsSimple, "intclust", average = FALSE) # no averaging: only one rep group
fGroupsEmpty <- getEmptyTestFGroups()
compsEmpty <- components()

test_that("components generation works", {
    # For RC/CAM: don't store their internal objects as they contain irreproducible file names
    expect_known_value(list(componentTable(compsRC), componentInfo(compsRC)), testFile("components-rc"))
    expect_known_value(list(componentTable(compsCAM), componentInfo(compsCAM)), testFile("components-cam"))
    expect_known_value(compsInt, testFile("components-int"))
    expect_length(compsEmpty, 0)
    expect_length(generateComponents(fGroupsEmpty, "ramclustr", ionization = "positive"), 0)
    expect_length(generateComponents(fGroupsEmpty, "camera", ionization = "positive"), 0)
    expect_length(generateComponents(fGroupsEmpty, "intclust"), 0)
    skip_if(length(compsNT) == 0)
    expect_known_value(compsNT, testFile("components-nt"))
    expect_length(generateComponents(fGroupsEmpty, "nontarget", ionization = "positive"), 0)
})

test_that("verify components show", {
    expect_known_show(compsRC, testFile("components-rc", text = TRUE))
    expect_known_show(compsCAM, testFile("components-cam", text = TRUE))
    expect_known_show(compsInt, testFile("components-int", text = TRUE))
    skip_if(length(compsNT) == 0)
    expect_known_show(compsNT, testFile("components-nt", text = TRUE))
})

test_that("basic subsetting", {
    expect_length(compsRC["nope"], 0)
    expect_equivalent(names(compsRC[1:2]), names(compsRC)[1:2])
    expect_equivalent(names(compsRC[names(compsRC)[1:2]]), names(compsRC)[1:2])
    expect_equivalent(names(compsRC[c(FALSE, TRUE)]), names(compsRC)[c(FALSE, TRUE)])
    expect_equivalent(groupNames(compsRC[, 1:2]), groupNames(compsRC)[1:2])
    expect_equivalent(groupNames(compsRC[, groupNames(compsRC)[2:3]]), groupNames(compsRC)[2:3])
    expect_equivalent(groupNames(compsRC[, c(FALSE, TRUE)]), groupNames(compsRC)[c(FALSE, TRUE)])
    expect_equal(length(compsRC[FALSE]), 0)
    expect_length(compsEmpty[1:5], 0)

    expect_equivalent(compsRC[[1, 1]], componentTable(compsRC)[[1]][group == groupNames(compsRC)[1]])
    expect_equivalent(compsRC[[names(compsRC)[1], groupNames(compsRC)[1]]], componentTable(compsRC)[[1]][group == groupNames(compsRC)[1]])
    expect_equivalent(compsRC[[2]], componentTable(compsRC)[[2]])
    expect_equivalent(callDollar(compsRC, names(compsRC)[3]), compsRC[[3]])
})

test_that("filtering works", {
    expect_length(filter(compsRC, size = c(0, 100)), length(compsRC))
    expect_length(filter(compsRC, size = c(50, 100)), 0)

    expect_length(filter(compsEmpty, size = c(0, 100)), 0)

    # shouldn't filter if related data is not there
    expect_equal(groupNames(filter(compsInt, adducts = TRUE)), groupNames(compsInt))
    expect_equal(groupNames(filter(compsInt, isotopes = TRUE)), groupNames(compsInt))

    expect_lt(length(groupNames(filter(compsRC, adducts = TRUE))), length(groupNames(compsRC)))
    expect_lt(length(groupNames(filter(compsRC, adducts = FALSE))), length(groupNames(compsRC)))
    expect_lt(length(groupNames(filter(compsRC, adducts = "[M+H]+"))), length(groupNames(compsRC)))
    expect_true(all(sapply(componentTable(filter(compsRC, adducts = "[M+H]+")),
                           function(cmp) all(cmp$adduct_ion == "[M+H]+"))))

    expect_lt(length(groupNames(filter(compsRC, isotopes = TRUE))), length(groupNames(compsRC)))
    expect_lt(length(groupNames(filter(compsRC, isotopes = FALSE))), length(groupNames(compsRC)))
    expect_lt(length(groupNames(filter(compsRC, isotopes = 0:1))), length(groupNames(compsRC)))
    expect_true(all(sapply(componentTable(filter(compsRC, isotopes = 0)),
                           function(cmp) all(cmp$isonr == 0))))

    skip_if(length(compsNT) == 0)
    expect_length(filter(compsNT, rtIncrement = c(0, 1000)), length(compsNT))
    expect_length(filter(compsNT, rtIncrement = c(100, 1000)), 0)
    expect_length(filter(compsNT, mzIncrement = c(0, 1000)), length(compsNT))
    expect_length(filter(compsNT, mzIncrement = c(1000, 10000)), 0)
})

test_that("basic usage works", {
    expect_equal(length(unique(as.data.table(compsCAM)$name)), length(compsCAM))

    expect_equivalent(findFGroup(compsCAM, compsCAM[[3]]$group[3]), 3)
    expect_length(findFGroup(compsCAM, "none"), 0)
    expect_length(findFGroup(compsEmpty, "1"), 0)
})

test_that("consensus works", {
    expect_length(consensus(compsRC, compsCAM), length(compsRC) + length(compsCAM))
    expect_error(consensus(compsRC, compsEmpty), "non-empty")
    expect_error(consensus(compsEmpty, components(componentInfo = data.table(),
                                                  algorithm = "empty2")),
                 "non-empty")
})

test_that("intensity clustered components", {
    expect_equivalent(length(treeCut(compsInt, k = 5)), 5)
    expect_equivalent(treeCutDynamic(compsInt), compsInt)
})

test_that("reporting works", {
    expect_file(reportCSV(fGroupsSimple, getWorkPath(), components = compsRC),
                getWorkPath("components.csv"))

    expect_file(reportPDF(fGroupsSimple, getWorkPath(), reportFGroups = FALSE, components = compsRC),
                getWorkPath("components.pdf"))

    # UNDONE: setting reportPlots="none" in combination with plotting components
    # crashes pandoc2, using "eics" for now...
    expect_reportMD(makeReportMD(fGroupsSimple, reportPlots = "eics", components = compsRC))
})

test_that("reporting empty object works", {
    expect_error(reportCSV(fGroupsSimple, getWorkPath(), components = compsEmpty), NA)
    expect_error(reportPDF(fGroupsSimple, getWorkPath(), reportFGroups = FALSE, components = compsEmpty), NA)
    expect_reportMD(makeReportMD(fGroupsSimple, reportPlots = "none", components = compsEmpty))
})

test_that("plotting works", {
    # specs don't work because of legend ... :(
    # expect_doppel("compon-spec", function() plotSpec(compsNT, 1))
    # expect_doppel("compon-spec-mark", function() plotSpec(compsNT, 1, markFGroup = names(fGroups)[1]))
    expect_plot(plotSpec(compsRC, 1))
    expect_plot(plotSpec(compsRC, 1, markFGroup = names(fGroupsSimple)[1]))
    expect_plot(print(plotSpec(compsRC, 1, useGGPlot2 = TRUE)))
    expect_doppel("eic-component", function() plotEIC(compsRC, 1, fGroupsSimple))

    expect_plot(plot(compsInt))
    expect_doppel("component-ic-int", function() plotInt(compsInt, index = 1))
    expect_doppel("component-ic-sil", function() plotSilhouettes(compsInt, 2:6))
    expect_doppel("component-ic-heat", function() plotHeatMap(compsInt, interactive = FALSE))
})
