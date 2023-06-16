context("prediction")

# UNDONE: test empties

# NOTE: fGroups, peak lists, annotations should all be in sync with test-formulas/test-compounds

# use example data from MS2Quant: the calibrants/eluent are completely unapplicable to patRoonData, but at least we can
# test calculations
calib <- fread(system.file("example_data", "quantification_example.csv", package = "MS2Quant"))
calib[, retention_time := retention_time * 60]
calib <- calib[!is.na(conc_M)]
if (testWithSets())
    calib <- list(calib, calib)
eluent <- fread(system.file("example_data", "eluent.csv", package = "MS2Quant"))
eluent[, time := time * 60]

fGroupsComps <- getCompFGroups()
fGroupsEmpty <- getEmptyTestFGroups()

fGroupsComps <- predictRespFactors(fGroupsComps, calib, eluent, organicModifier = "MeOH", pHAq = 4,
                                   calibConcUnit = "M")
fGroupsComps <- predictTox(fGroupsComps)

doSIRIUS <- !is.null(getOption("patRoon.path.SIRIUS")) && nzchar(getOption("patRoon.path.SIRIUS"))
if (doSIRIUS)
{
    fGroupsForms <- getFormFGroups()
    fGroupsForms <- predictRespFactors(fGroupsForms, calib, eluent, organicModifier = "MeOH", pHAq = 4,
                                       calibConcUnit = "M")
    fGroupsForms <- predictTox(fGroupsForms)
    
    plistsForms <- generateMSPeakLists(fGroupsForms, "mzr")
    plistsComps <- generateMSPeakLists(fGroupsComps, "mzr")
    
    formsSIR <- doGenFormsSIRFPs(fGroupsForms, plistsForms)
    formsSIR <- predictRespFactors(formsSIR, fGroupsForms, calib, eluent, organicModifier = "MeOH", pHAq = 4,
                                   calibConcUnit = "M")
    formsSIR <- predictTox(formsSIR)
    formsTab <- as.data.table(formsSIR)
    
    compsSIR <- doGenCompsSIR(fGroupsComps, plistsComps)
    compsSIR <- predictRespFactors(compsSIR, fGroupsComps, calib, eluent, organicModifier = "MeOH", pHAq = 4,
                                   calibConcUnit = "M", type = "both")
    compsSIR <- predictTox(compsSIR, type = "both")
    compsTab <- as.data.table(compsSIR)
}

test_that("Basics for prediction", {
    expect_known_value(screenInfo(fGroupsComps)[, c("RF_SMILES", "LC50_SMILES"), with = FALSE], testFile("pred-scr"))
    checkmate::expect_names(names(screenInfo(fGroupsComps)), must.include = c("RF_SMILES", "LC50_SMILES"))
    checkmate::expect_data_table(screenInfo(fGroupsComps)[, c("RF_SMILES", "LC50_SMILES"), with = FALSE],
                                 types = "numeric")
    checkmate::expect_numeric(screenInfo(filter(fGroupsComps, minRF = 1E6))$RF_SMILES, lower = 1E6)
    checkmate::expect_numeric(screenInfo(filter(fGroupsComps, maxLC50 = 1E5))$LC50_SMILES, upper = 1E5)
    
    skip_if_not(doSIRIUS)

    expect_known_value(formsTab[, grep("^(RF|LC50)_", names(formsTab), value = TRUE), with = FALSE], testFile("pred-forms"))
    expect_known_value(compsTab[, grep("^(RF|LC50)_", names(compsTab), value = TRUE), with = FALSE], testFile("pred-comps"))
    
    expect_length(grep("^(RF|LC50)_", names(formsTab)), if (testWithSets()) 4 else 2)
    expect_length(grep("^(RF|LC50)_", names(compsTab)), if (testWithSets()) 6 else 3)
    expect_length(grep("^(RF|LC50)_", compsSIR@scoreTypes), if (testWithSets()) 6 else 3)
})

if (doSIRIUS)
{
    fGroupsFormsC <- calculateConcs(fGroupsForms, formsSIR)
    fGroupsFormsC <- calculateTox(fGroupsFormsC, formsSIR)
    fGroupsCompsC <- calculateConcs(fGroupsComps, compsSIR)
    fGroupsCompsC <- calculateTox(fGroupsCompsC, compsSIR)
    calcTab <- as.data.table(fGroupsCompsC)
    calcTabNoPref <- as.data.table(fGroupsCompsC, concAggrParams = getDefPredAggrParams(preferType = "none"),
                                   toxAggrParams = getDefPredAggrParams(preferType = "none"))
    calcTabMax <- as.data.table(fGroupsCompsC, concAggrParams = getDefPredAggrParams(max, preferType = "none"),
                                toxAggrParams = getDefPredAggrParams(max, preferType = "none"))
}

test_that("Basics for calculation", {
    skip_if_not(doSIRIUS)
    
    expect_known_value(list(concentrations(fGroupsFormsC), toxicities(fGroupsFormsC)), testFile("calc-forms"))
    expect_known_value(list(concentrations(fGroupsCompsC), toxicities(fGroupsCompsC)), testFile("calc-comps"))
    checkmate::expect_data_table(concentrations(fGroupsFormsC), min.rows = 1)
    checkmate::expect_data_table(toxicities(fGroupsFormsC), min.rows = 1)
    checkmate::expect_data_table(concentrations(fGroupsCompsC), min.rows = 1)
    checkmate::expect_data_table(toxicities(fGroupsCompsC), min.rows = 1)
    checkmate::expect_subset(c("suspect", "SIRIUS_FP"), concentrations(fGroupsFormsC)$type)
    checkmate::expect_subset(c("suspect", "SIRIUS_FP"), toxicities(fGroupsFormsC)$type)
    checkmate::expect_subset(c("suspect", "SIRIUS_FP", "compound"), concentrations(fGroupsCompsC)$type)
    checkmate::expect_subset(c("suspect", "SIRIUS_FP", "compound"), toxicities(fGroupsCompsC)$type)
    
    expect_setequal("suspect", calcTab$conc_types) # default preferType == "suspect and all suspects should have results
    expect_setequal(c("suspect", "SIRIUS_FP", "compound"), unlist(strsplit(calcTabNoPref$conc_types, ",")))
    expect_lt(mean(calcTabNoPref[["standard-pos-2_conc"]], na.rm = TRUE),
               mean(calcTabMax[["standard-pos-2_conc"]], na.rm = TRUE))
    
    expect_setequal("suspect", calcTab$LC50_types)
    expect_setequal(c("suspect", "SIRIUS_FP", "compound"), unlist(strsplit(calcTabNoPref$LC50_types, ",")))
    expect_lt(mean(calcTabNoPref$LC50, na.rm = TRUE), mean(calcTabMax$LC50, na.rm = TRUE))
})
