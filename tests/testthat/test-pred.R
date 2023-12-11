context("prediction")

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
fGroupsEmpty <- delete(fGroupsComps) # don't use getEmptyTestFGroups(): we want to have a screening object

fGroupsComps <- predictRespFactors(fGroupsComps, calib, eluent, organicModifier = "MeOH", pHAq = 4,
                                   calibConcUnit = "M")
fGroupsComps <- predictTox(fGroupsComps)

fGroupsEmpty <- predictRespFactors(fGroupsEmpty, list(), eluent, organicModifier = "MeOH", pHAq = 4,
                                   calibConcUnit = "M")
fGroupsEmpty <- predictTox(fGroupsEmpty)

doSIRIUS <- TRUE #!is.null(getOption("patRoon.path.SIRIUS")) && nzchar(getOption("patRoon.path.SIRIUS"))
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
    
    formsEmpty <- delete(formsSIR)
    compsEmpty <- delete(compsSIR)
}

getMS2QLM <- function(obj) if (testWithSets()) setObjects(obj)[[1]]@MS2QuantMeta$linModel else obj@MS2QuantMeta$linModel

test_that("Basics for prediction", {
    expect_known_value(screenInfo(fGroupsComps)[, c("RF_SMILES", "LC50_SMILES"), with = FALSE], testFile("pred-scr"))
    checkmate::expect_names(names(screenInfo(fGroupsComps)), must.include = c("RF_SMILES", "LC50_SMILES"))
    checkmate::expect_data_table(screenInfo(fGroupsComps)[, c("RF_SMILES", "LC50_SMILES"), with = FALSE],
                                 types = "numeric")
    checkmate::expect_numeric(screenInfo(filter(fGroupsComps, minRF = 1E6))$RF_SMILES, lower = 1E6)
    checkmate::expect_numeric(screenInfo(filter(fGroupsComps, maxLC50 = 1E5))$LC50_SMILES, upper = 1E5)
    expect_equal(nrow(screenInfo(fGroupsEmpty)), 0)
    
    skip_if_not(doSIRIUS)

    expect_known_value(formsTab[, grep("^(RF|LC50)_", names(formsTab), value = TRUE), with = FALSE], testFile("pred-forms"))
    expect_known_value(compsTab[, grep("^(RF|LC50)_", names(compsTab), value = TRUE), with = FALSE], testFile("pred-comps"))
    
    expect_length(grep("^(RF|LC50)_", names(formsTab)), if (testWithSets()) 4 else 2)
    expect_length(grep("^(RF|LC50)_", names(compsTab)), if (testWithSets()) 6 else 3)
    expect_length(grep("^(RF|LC50)_", compsSIR@scoreTypes), if (testWithSets()) 6 else 3)
    
    expect_error(predictRespFactors(formsEmpty, fGroupsComps, calib, eluent, organicModifier = "MeOH", pHAq = 4,
                                    calibConcUnit = "M", type = "both"), NA)
    expect_error(predictTox(formsEmpty), NA)
    expect_error(predictRespFactors(compsEmpty, fGroupsComps, calib, eluent, organicModifier = "MeOH", pHAq = 4,
                                    calibConcUnit = "M", type = "both"), NA)
    expect_error(predictTox(compsEmpty), NA)
    
    expect_equal(getMS2QLM(fGroupsForms), getMS2QLM(fGroupsComps))
    expect_equal(getMS2QLM(fGroupsForms), getMS2QLM(formsSIR))
    expect_equal(getMS2QLM(fGroupsForms), getMS2QLM(compsSIR))
    
})

calibConcs <- data.table(name = "Chloridazon", "standard-pos" = 100)
calibScr <- getQuantCalibFromScreening(fGroupsComps, calibConcs)
calibScrAvg <- getQuantCalibFromScreening(fGroupsComps, calibConcs, average = TRUE)
calibScrAr <- getQuantCalibFromScreening(fGroupsComps, calibConcs, areas = TRUE)
test_that("getQuantCalibFromScreening()", {
    checkmate::expect_data_table(calibScr, nrows = 2)
    checkmate::expect_names(names(calibScr), permutation.of = c("name", "SMILES", "rt", "conc", "intensity"))
    checkmate::expect_data_table(calibScrAvg, nrows = 1)
    expect_equal(mean(calibScr$intensity), calibScrAvg$intensity)
    checkmate::expect_data_table(calibScrAr, nrows = 2)
    expect_true(all(calibScrAr$intensity > calibScr$intensity))
    expect_equal(calibScr[, -"intensity"], calibScrAr[, -"intensity"], check.attributes = FALSE)
})

if (doSIRIUS)
{
    fGroupsFormsC <- calculateConcs(fGroupsForms, formsSIR)
    fGroupsFormsC <- calculateTox(fGroupsFormsC, formsSIR)
    fGroupsCompsC <- calculateConcs(fGroupsComps, compsSIR)
    fGroupsCompsC <- calculateTox(fGroupsCompsC, compsSIR)
    fGroupsOnlyCompsC <- calculateConcs(getCompFGroups(), compsSIR)
    fGroupsOnlyCompsC <- calculateTox(fGroupsOnlyCompsC, compsSIR)
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

    expect_equal(calculateConcs(fGroupsComps, compsEmpty), fGroupsComps)
    expect_equal(calculateConcs(fGroupsEmpty, compsEmpty), fGroupsEmpty)
    expect_equal(calculateTox(fGroupsComps, compsEmpty), fGroupsComps)
    expect_equal(calculateTox(fGroupsEmpty, compsEmpty), fGroupsEmpty)
    
    expect_setequal("suspect", calcTab$conc_types) # default preferType == "suspect and all suspects should have results
    expect_setequal(c("suspect", "SIRIUS_FP", "compound"), unlist(strsplit(calcTabNoPref$conc_types, ",")))
    expect_lt(mean(calcTabNoPref[["standard-pos-2_conc"]], na.rm = TRUE),
               mean(calcTabMax[["standard-pos-2_conc"]], na.rm = TRUE))
    
    expect_setequal("suspect", calcTab$LC50_types)
    expect_setequal(c("suspect", "SIRIUS_FP", "compound"), unlist(strsplit(calcTabNoPref$LC50_types, ",")))
    expect_lt(mean(calcTabNoPref$LC50, na.rm = TRUE), mean(calcTabMax$LC50, na.rm = TRUE))
    
    # NOTE: use as.data.table here instead of slot data, as these like the filters use aggregated data
    expect_gte(min(as.data.table(filter(fGroupsCompsC, absMinConc = 0.2))[, paste0(analyses(fGroupsCompsC), "_conc"), with = FALSE],
                   na.rm = TRUE), 0.2)
    expect_lt(min(as.data.table(filter(fGroupsCompsC, absMinConc = 0.2, negate = TRUE))[, paste0(analyses(fGroupsCompsC), "_conc"), with = FALSE],
                   na.rm = TRUE), 0.2)
    expect_gte(min(as.data.table(filter(fGroupsCompsC, absMinConcTox = 1E-6), normConcToTox = TRUE)[, paste0(analyses(fGroupsCompsC), "_conc"), with = FALSE],
                   na.rm = TRUE), 1E-6)
    expect_lt(min(as.data.table(filter(fGroupsCompsC, absMinConcTox = 1E-6, negate = TRUE), normConcToTox = TRUE)[, paste0(analyses(fGroupsCompsC), "_conc"), with = FALSE],
                  na.rm = TRUE), 1E-6)
    expect_lte(max(as.data.table(filter(fGroupsCompsC, absMaxTox = 6E4))$LC50, na.rm = TRUE), 6E4)
    expect_gt(max(as.data.table(filter(fGroupsCompsC, absMaxTox = 6E4, negate = TRUE))$LC50, na.rm = TRUE), 6E4)
    
    # NOTE: use fGroupsOnlyCompsC as fGroupsCompsC has suspect results and therefore results for everything
    expect_gt(length(getFeatures(fGroupsOnlyCompsC)), length(getFeatures(filter(fGroupsOnlyCompsC, absMinConc = 0.2))))
    expect_gt(length(getFeatures(filter(fGroupsOnlyCompsC, absMinConc = 0.2))),
              length(getFeatures(filter(fGroupsOnlyCompsC, absMinConc = 0.2, removeNA = TRUE))))
    expect_gt(length(getFeatures(fGroupsOnlyCompsC)), length(getFeatures(filter(fGroupsOnlyCompsC, absMaxTox = 6E4))))
    expect_gt(length(getFeatures(fGroupsOnlyCompsC)), length(getFeatures(filter(fGroupsOnlyCompsC, absMinConcTox = 0.01))))
    expect_gt(length(getFeatures(filter(fGroupsOnlyCompsC, absMinConcTox = 0.01))),
              length(getFeatures(filter(fGroupsOnlyCompsC, absMinConcTox = 0.01, removeNA = TRUE))))
    # UNDONE: can we actually get NA tox values?
    # expect_gt(length(getFeatures(filter(fGroupsOnlyCompsC, absMaxTox = 6E4))),
    #           length(getFeatures(filter(fGroupsOnlyCompsC, absMaxTox = 6E4, removeNA = TRUE))))
})
