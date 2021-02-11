context("compounds clustering")

hasMetfrag <- !is.null(getOption("patRoon.path.MetFragCL")) && nzchar(getOption("patRoon.path.MetFragCL"))

if (hasMetfrag)
{
    fGroups <- getCompFGroups()
    plists <- generateMSPeakLists(fGroups, "mzr")
    compounds <- callMF(fGroups, plists, db = file.path(getTestDataPath(), "test-mf-db-isomers.csv"))
    compsClust <- makeHCluster(compounds)
    firstGroup <- names(clusters(compsClust))[1]

    compsEmpty <- filter(compounds, minFragScore = 1E4)
    compsClustEmpty <- makeHCluster(compsEmpty)
}

test_that("verify compound cluster generation", {
    skip_if_not(hasMetfrag)

    expect_known_value(compsClust, testFile("compounds-clust"))
    expect_known_show(compsClust, testFile("compounds-clust", text = TRUE))
    # should have clusters for same number of feature groups with compounds
    expect_length(clusters(compsClust), length(compoundTable(compounds)))

    expect_length(compsClustEmpty, 0)
    expect_known_show(compsClustEmpty, testFile("compounds-clust-empty", text = TRUE))
})

test_that("basic subsetting", {
    skip_if_not(hasMetfrag)

    expect_length(compsClust["nope"], 0)
    expect_equivalent(groupNames(compsClust[1:2]), groupNames(compsClust)[1:2])
    expect_equivalent(groupNames(compsClust[groupNames(compsClust)[2:3]]), groupNames(compsClust)[2:3])
    expect_equivalent(groupNames(compsClust[c(FALSE, TRUE)]), groupNames(compsClust)[c(FALSE, TRUE)])
    expect_equal(length(compsClust[FALSE]), 0)
    expect_length(compsClustEmpty[1:5], 0)
})

test_that("override cutting clusters work", {
    skip_if_not(hasMetfrag)

    expect_equivalent(lengths(treeCut(compsClust, k = 2, groupName = firstGroup))[1], 2)
    expect_equal(treeCutDynamic(compsClust, groupName = firstGroup), compsClust)
})

test_that("plotting works", {
    skip_if_not(hasMetfrag)

    # these don't play well with vdiffr...
    expect_plot(plot(compsClust, groupName = firstGroup))
    expect_plot(plotStructure(compsClust, groupName = firstGroup, cluster = 1))

    # expect_doppel("compounds-clust-sil", function() plotSilhouettes(compsClust, kSeq = 2:6, groupName = firstGroup))
    # UNDONE: currently all clusters have only two compounds so can't plot silhouette which needs >2
    expect_error(plotSilhouettes(compsClust, kSeq = 2:6, groupName = firstGroup))
})

test_that("reporting works", {
    skip_if_not(hasMetfrag)

    # also include reporting compounds: cluster number is added to candidate info

    expect_error(reportCSV(fGroups, getWorkPath(), compounds = compounds, compsCluster = compsClust), NA)
    for (grp in names(clusters(compsClust)))
        expect_csv_file(getWorkPath("compounds", sprintf("%s-%s.csv", class(fGroups), grp)), "cluster")

    expect_error(reportPDF(fGroups, getWorkPath(), reportFGroups = FALSE, compounds = compounds,
                           MSPeakLists = plists, compsCluster = compsClust), NA)
    for (grp in names(clusters(compsClust)))
        checkmate::expect_file_exists(getWorkPath("compounds", sprintf("%s-%s-clusters.pdf", class(fGroups), grp))) # UNDONE: check col name

    expect_reportHTML(makeReportHTML(fGroups, reportPlots = "none", compounds = compounds, MSPeakLists = plists, compsCluster = compsClust))
})

test_that("reporting empty object works", {
    skip_if_not(hasMetfrag)

    expect_error(reportCSV(fGroups, getWorkPath(), compsCluster = compsClustEmpty), NA)
    expect_error(reportPDF(fGroups, getWorkPath(), compsCluster = compsClustEmpty), NA)
    expect_reportHTML(makeReportHTML(fGroups, reportPlots = "none", compsCluster = compsClustEmpty))
})
