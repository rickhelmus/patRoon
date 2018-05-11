context("components")

fGroups <- getTestFGroups(getTestAnaInfo()[4:5, ])

# fix seed for reproducible clustering
withr::with_seed(20, compsRC <- generateComponents(fGroups, "ramclustr", ionization = "positive"))
# UNDONE: getting unknown NaN warnings here...
suppressWarnings(compsCAM <- generateComponents(fGroups, "camera", ionization = "positive"))
compsNT <- generateComponents(fGroups, "nontarget", ionization = "positive")
fGroupsEmpty <- getEmptyTestFGroups()
compsEmpty <- components()

test_that("components generation works", {
    # For RC/CAM: don't store their internal objects as they contain irreproducible file names
    expect_known_value(list(componentTable(compsRC), componentInfo(compsRC)), testFile("components-rc"))
    expect_known_value(list(componentTable(compsCAM), componentInfo(compsCAM)), testFile("components-cam"))
    expect_length(compsEmpty, 0)
    expect_length(generateComponents(fGroupsEmpty, "ramclustr", ionization = "positive"), 0)
    expect_length(generateComponents(fGroupsEmpty, "camera", ionization = "positive"), 0)
    skip_if(length(compsNT) == 0)
    expect_known_value(compsNT, testFile("components-nt"))
    expect_length(generateComponents(fGroupsEmpty, "nontarget", ionization = "positive"), 0)
})

test_that("verify components show", {
    expect_known_show(compsRC, testFile("components-rc", text = TRUE))
    expect_known_show(compsCAM, testFile("components-cam", text = TRUE))
    skip_if(length(compsNT) == 0)
    expect_known_show(compsNT, testFile("components-nt", text = TRUE))
})

test_that("findFGroup works", {
    expect_equivalent(findFGroup(compsCAM, names(fGroups)[1]), 1)
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

test_that("reporting works", {
    expect_file(reportCSV(fGroups, getWorkPath(), components = compsRC),
                getWorkPath("components.csv"))

    expect_file(reportPDF(fGroups, getWorkPath(), reportFGroups = FALSE, components = compsRC),
                getWorkPath("components.pdf"))
    
    expect_file(reportMD(fGroups, getWorkPath(), reportChord = FALSE, reportFGroups = FALSE,
                         components = compsRC),
                getWorkPath("report.html"))
})

test_that("reporting empty object works", {
    expect_error(reportCSV(fGroups, getWorkPath(), components = compsEmpty), NA)
    expect_error(reportPDF(fGroups, getWorkPath(), reportFGroups = FALSE, components = compsEmpty), NA)
    expect_file(reportMD(fGroups, getWorkPath(), reportChord = FALSE, reportFGroups = FALSE,
                         components = compsEmpty),
                getWorkPath("report.html"))
})

test_that("plotting works", {
    # specs don't work because of legend ... :(
    # expect_doppel("compon-spec", function() plotSpec(compsNT, 1))
    # expect_doppel("compon-spec-mark", function() plotSpec(compsNT, 1, markFGroup = names(fGroups)[1]))
    expect_plot(plotSpec(compsRC, 1))
    expect_plot(plotSpec(compsRC, 1, markFGroup = names(fGroups)[1]))
    expect_plot(print(plotSpec(compsRC, 1, useGGPlot2 = TRUE)))
    expect_doppel("eic-component", function() plotEIC(compsRC, 1, fGroups))
})
