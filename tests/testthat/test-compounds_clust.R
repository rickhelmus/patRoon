context("compounds clustering")

hasMetfrag <- !is.null(getOption("patRoon.path.metFragCL")) && nzchar(getOption("patRoon.path.metFragCL"))

if (hasMetfrag)
{
    fGroups <- getTestFGroups(getTestAnaInfo()[4, ])[, 1:25]
    plists <- generateMSPeakLists(fGroups, "mzr")
    compounds <- generateCompounds(fGroups, plists, "metfrag", logPath = NULL,
                                   adduct = 1, isPositive = TRUE)
    
    compsClust <- makeHCluster(compounds)
    firstGroup <- names(clusters(compsClust))[1]
    
    # for reference: clearout molecules as these can't be stored well
    compsClustRef <- compsClust
    compsClustRef@molecules <- list()
}

test_that("verify compound cluster generation", {
    skip_if_not(hasMetfrag)

    expect_known_value(compsClustRef, testFile("compounds-clust"))
    expect_known_show(compsClust, testFile("compounds-clust", text = TRUE))
    # should have clusters for same number of feature groups with compounds
    expect_length(clusters(compsClust), length(compoundTable(compounds)))
})

test_that("override cutting clusters work", {
    skip_if_not(hasMetfrag)
    
    expect_equivalent(lengths(treeCut(compsClust, k = 5, groupName = firstGroup))[1], 5)
    expect_equal(treeCutDynamic(compsClust, groupName = firstGroup), compsClust)
})

test_that("plotting works", {
    skip_if_not(hasMetfrag)
    
    # these don't play well with vdiffr...
    expect_plot(plot(compsClust, groupName = firstGroup))
    expect_plot(plotStructure(compsClust, groupName = firstGroup, cluster = 1))
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
    
    expect_file(reportMD(fGroups, getWorkPath(), reportChord = FALSE, reportFGroups = FALSE,
                         compounds = compounds, MSPeakLists = plists),
                getWorkPath("report.html"))
})
