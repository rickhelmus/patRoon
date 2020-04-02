context("screening")

susps <- as.data.table(patRoonData::targets)

susps[, InChI := babelConvert(SMILES, "smi", "inchi")]
susps[, neutralMass := getNeutralMassFromSMILES(SMILES)]
susps[, formula := sapply(getMoleculesFromSMILES(SMILES), function(m) rcdk::get.mol2formula(m)@string)]
susps[, adduct := "[M+H]+"]
susps[name %in% c("TBA", "TPA"), adduct := "[M]+"]

fGroups <- getTestFGroups(getTestAnaInfo())
scr <- screenSuspects(fGroups, susps)
scrF <- screenSuspects(getFeatures(fGroups), susps)
scrSMI <- screenSuspects(fGroups, susps[, c("name", "rt", "adduct", "SMILES")])

suspsMissing <- copy(susps)
suspsMissing[1, mz := NA]
suspsMissing[2, neutralMass := NA]
suspsMissing[3, formula := NA_character_]
suspsMissing[4, SMILES := ""]
suspsMissing[5, InChI := ""]
suspsMissingRow <- copy(susps)
suspsMissingRow[2, c("mz", "neutralMass", "formula", "SMILES", "InChI") := NA]

test_that("suspect screening is OK", {
    expect_equal(nrow(scr), nrow(susps))
    expect_length(unique(scrF$name), nrow(susps))
    expect_known_value(scr, testFile("screening"))
    expect_known_value(scrF[, -"feature"], testFile("screening-feat"))
    expect_length(groupFeaturesScreening(fGroups, scr), nrow(susps))

    # check suspects without retention
    expect_gte(nrow(screenSuspects(fGroups, susps[, -3])), nrow(scr))
    expect_gte(nrow(screenSuspects(getFeatures(fGroups), susps[, -3])), nrow(scrF))
    
    # alternative ion mass calculations
    expect_equal(scr, scrSMI, tolerance = 1E-3) # higher tolerance due to inaccurate mz column in patRoonData::targets
    expect_equal(scrSMI, screenSuspects(fGroups, susps[, c("name", "rt", "adduct", "InChI")]))
    expect_equal(scrSMI, screenSuspects(fGroups, susps[, c("name", "rt", "adduct", "neutralMass")]))
    expect_equal(scrSMI, screenSuspects(fGroups, susps[, c("name", "rt", "adduct", "formula")]))

    # same, with missing data (having 2 options for ion mass calculation should be sufficient)
    expect_equal(scrSMI, screenSuspects(fGroups, suspsMissing[, c("name", "rt", "adduct", "mz", "neutralMass")]),
                 tolerance = 1E-3)
    expect_equal(scrSMI, screenSuspects(fGroups, suspsMissing[, c("name", "rt", "adduct", "neutralMass", "formula")]))
    expect_equal(scrSMI, screenSuspects(fGroups, suspsMissing[, c("name", "rt", "adduct", "formula", "SMILES")]))
    expect_equal(scrSMI, screenSuspects(fGroups, suspsMissing[, c("name", "rt", "adduct", "SMILES", "InChI")]))
    
    expect_warning(screenSuspects(fGroups, suspsMissingRow, skipInvalid = TRUE))
    expect_error(screenSuspects(fGroups, suspsMissingRow, skipInvalid = FALSE))
    
    # adduct argument
    expect_equal(screenSuspects(fGroups, susps[name %in% c("TBA", "TPA")]),
                 screenSuspects(fGroups, susps[name %in% c("TBA", "TPA"), -"adduct"], adduct = "[M]+"))
})

TQFile <- file.path(getTestDataPath(), "GlobalResults-TASQ.csv")
TQRes <- fread(TQFile)
fGroupsTQ <- importFeatureGroupsBrukerTASQ(TQFile, getTestAnaInfo())
fGroupsTQ <- filter(fGroupsTQ, blankThreshold = 5, removeBlanks = TRUE)
TQFileNoRT <- file.path(getTestDataPath(), "GlobalResults_noRT-TASQ.csv")
TQResNoRT <- fread(TQFileNoRT)
fGroupsTQNoRT <- importFeatureGroupsBrukerTASQ(TQFileNoRT, getTestAnaInfo())
fGroupsTQNoRT <- filter(fGroupsTQNoRT, blankThreshold = 5, removeBlanks = TRUE)

test_that("TASQ import works", {
    expect_equal(unique(TQRes$Analyte), names(fGroupsTQ))
    expect_known_value(groups(fGroupsTQ), testFile("susp-tasq"))
    expect_known_show(fGroupsTQ, testFile("susp-tasq", text = TRUE))
    
    expect_lt(length(unique(TQResNoRT$Analyte)), length(fGroupsTQNoRT))
    expect_known_value(groups(fGroupsTQNoRT), testFile("susp-tasq-nort"))
    expect_known_show(fGroupsTQNoRT, testFile("susp-tasq-nort", text = TRUE))
    
})
