context("formulas")

fGroups <- getTestFGroups(getTestAnaInfo()[4:5, ])
plists <- generateMSPeakLists(fGroups, "mzr")

doSIRIUS <- !is.null(getOption("patRoon.path.SIRIUS")) && nzchar(getOption("patRoon.path.SIRIUS"))

formsGF <- generateFormulas(fGroups, "genform", plists)

if (doSIRIUS)
    formsSIR <- generateFormulas(fGroups, "sirius", plists, logPath = NULL)

test_that("verify formula generation", {
    expect_known_value(formsGF, testFile("formulas-gf"))
    skip_if_not(doSIRIUS)
    expect_known_value(formsSIR, testFile("formulas-sir"))
})

test_that("verify show output", {
    expect_known_show(formsGF, testFile("formulas-gf", text = TRUE))
    skip_if_not(doSIRIUS)
    expect_known_show(formsSIR, testFile("formulas-sir", text = TRUE))
})

fCons <- consensus(formsGF, fGroups = fGroups)
if (doSIRIUS)
    fCons2 <- consensus(formsGF, formsSIR, fGroups = fGroups)

fTable <- formulaTable(fCons)

test_that("consensus works", {
    expect_known_value(fCons, testFile("formulas-cons"))
    expect_known_show(fCons, testFile("formulas-cons", text = TRUE))
    expect_error(consensus(formsGF, formsGF, fGroups = fGroups)) # only support different algos at this point
    expect_gt(length(consensus(formsGF, fGroups = fGroups, formAnaThreshold = 0)), length(fCons))
    expect_lt(length(consensus(formsGF, fGroups = fGroups, maxFormulas = 1)), length(fCons))
    expect_lt(length(consensus(formsGF, fGroups = fGroups, maxFragFormulas = 1)), length(fCons))
    expect_lte(length(consensus(formsGF, fGroups = fGroups, maxFormulas = 1, maxFragFormulas = 1)),
               length(unique(fTable$group)) * 2) # * 2: one MS + one MSMS formula max
    expect_gte(min(formulaTable(consensus(formsGF, fGroups = fGroups, minIntensity = 500))$min_intensity), 500)
    expect_lte(min(formulaTable(consensus(formsGF, fGroups = fGroups, maxIntensity = 10000))$min_intensity), 10000)
    checkmate::expect_names(names(formulaTable(consensus(formsGF, fGroups = fGroups, elements = c("C", "H")))),
                            must.include = c("C", "H"))
    checkmate::expect_names(names(formulaTable(consensus(formsGF, fGroups = fGroups, fragElements = c("C", "H")))),
                            must.include = c("frag_C", "frag_H"))
    
    skip_if_not(doSIRIUS)
    expect_known_value(fCons2, testFile("formulas-cons2"))
    expect_known_show(fCons2, testFile("formulas-cons2", text = TRUE))
    expect_lt(length(consensus(formsGF, formsSIR, fGroups = fGroups, formListThreshold = 1)), length(fCons2))
})

test_that("feature group filtering", {
    expect_setequal(names(filter(fGroups, formConsensus = fCons)), fTable$group)
})

test_that("reporting works", {
    expect_file(reportCSV(fGroups, getWorkPath(), formConsensus = fCons),
                getWorkPath("formulas.csv"))
    
    expect_error(reportPDF(fGroups, getWorkPath(), reportFGroups = FALSE, formConsensus = fCons,
                           MSPeakLists = plists), NA)
    for (grp in unique(fTable[byMSMS == TRUE, group]))
        checkmate::expect_file_exists(getWorkPath("formulas", sprintf("%s-%s.pdf", class(fGroups), grp)))
    
    expect_file(reportMD(fGroups, getWorkPath(), reportChord = FALSE, reportFGroups = FALSE,
                         formConsensus = fCons, MSPeakLists = plists),
                getWorkPath("report.html"))
})

test_that("plotting works", {
    expect_docker("form-spec", function() plotSpec(fCons, fTable[byMSMS == TRUE, formula][1],
                                              fTable[byMSMS == TRUE, group][1], plists))
    
    # ggplot2 versions don't really work with vdiffr at the moment :(
    # expect_docker("spec-gg", plotSpec(fCons, fTable[byMSMS == TRUE, formula][1],
    #                                                 fTable[byMSMS == TRUE, group][1], plists,
    #                                                 useGGPlot2 = TRUE))
    expect_plot(print(plotSpec(fCons, fTable[byMSMS == TRUE, formula][1],
                               fTable[byMSMS == TRUE, group][1], plists,
                               useGGPlot2 = TRUE)))
})
