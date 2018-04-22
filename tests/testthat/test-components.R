context("components")

fGroups <- getTestFGroups(getTestAnaInfo()[4:5, ])

if (Sys.info()[["sysname"]] == "Windows")
{
    library(nontarget)
    detach("package:nontarget", unload = TRUE)
}

# fix seed for reproducible clustering
withr::with_seed(20, compsRC <- generateComponents(fGroups, "ramclustr", ionization = "positive"))
# UNDONE: getting unknown NaN warnings here...
suppressWarnings(compsCAM <- generateComponents(fGroups, "camera", ionization = "positive"))
compsNT <- generateComponents(fGroups, "nontarget", ionization = "positive")

test_that("components generation works", {
    # For RC/CAM: don't store their internal objects as they contain irreproducible file names
    expect_known_value(list(componentTable(compsRC), componentInfo(compsRC)), testFile("components-rc"))
    expect_known_value(list(componentTable(compsCAM), componentInfo(compsCAM)), testFile("components-cam"))
    expect_known_value(compsNT, testFile("components-nt"))
})

test_that("verify components show", {
    expect_known_show(compsRC, testFile("components-rc", text = TRUE))
    expect_known_show(compsCAM, testFile("components-cam", text = TRUE))
    expect_known_show(compsNT, testFile("components-nt", text = TRUE))
})

test_that("findFGroup works", {
    expect_equivalent(findFGroup(compsCAM, names(fGroups)[1]), 1)
})

test_that("consensus works", {
    expect_equal(length(consensus(compsRC, compsNT)), length(compsRC) + length(compsNT))
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

test_that("plotting works", {
    # specs don't work because of legend ... :(
    # expect_docker("compon-spec", function() plotSpec(compsNT, 1))
    # expect_docker("compon-spec-mark", function() plotSpec(compsNT, 1, markFGroup = names(fGroups)[1]))
    expect_plot(plotSpec(compsNT, 1))
    expect_plot(plotSpec(compsNT, 1, markFGroup = names(fGroups)[1]))
    expect_plot(print(plotSpec(compsNT, 1, useGGPlot2 = TRUE)))
    expect_docker("eic", function() plotEIC(compsNT, 1, fGroups))
})
