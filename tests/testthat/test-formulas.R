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
    expect_known_output(show(formsGF), testFile("formulas-gf", text = TRUE))
    skip_if_not(doSIRIUS)
    expect_known_output(show(formsSIR), testFile("formulas-sir", text = TRUE))
})

fCons <- consensus(formsGF, fGroups = fGroups)
if (doSIRIUS)
    fCons2 <- consensus(formsGF, formsSIR, fGroups = fGroups)

test_that("consensus works", {
    expect_known_value(fCons, testFile("formulas-cons"))
    expect_known_output(show(fCons), testFile("formulas-cons", text = TRUE))
    expect_error(consensus(formsGF, formsGF, fGroups = fGroups)) # only support different algos at this point
    expect_gt(length(consensus(formsGF, fGroups = fGroups, formAnaThreshold = 0)), length(fCons))
    expect_lt(length(consensus(formsGF, fGroups = fGroups, maxFormulas = 1)), length(fCons))
    expect_lt(length(consensus(formsGF, fGroups = fGroups, maxFragFormulas = 1)), length(fCons))
    expect_lte(length(consensus(formsGF, fGroups = fGroups, maxFormulas = 1, maxFragFormulas = 1)),
               length(unique(formulaTable(fCons)$group)) * 2) # * 2: one MS + one MSMS formula max
    expect_gte(min(formulaTable(consensus(formsGF, fGroups = fGroups, minIntensity = 500))$min_intensity), 500)
    expect_lte(min(formulaTable(consensus(formsGF, fGroups = fGroups, maxIntensity = 10000))$min_intensity), 10000)
    checkmate::expect_names(names(formulaTable(consensus(formsGF, fGroups = fGroups, elements = c("C", "H")))),
                            must.include = c("C", "H"))
    checkmate::expect_names(names(formulaTable(consensus(formsGF, fGroups = fGroups, fragElements = c("C", "H")))),
                            must.include = c("frag_C", "frag_H"))
    
    skip_if_not(doSIRIUS)
    expect_known_value(fCons2, testFile("formulas-cons2"))
    expect_known_output(show(fCons2), testFile("formulas-cons2", text = TRUE))
    expect_lt(length(consensus(formsGF, formsSIR, fGroups = fGroups, formListThreshold = 1)), length(fCons2))
})

test_that("feature group filtering", {
    expect_setequal(names(filter(fGroups, formConsensus = fCons)), formulaTable(fCons)$group)
})
