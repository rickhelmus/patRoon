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
    
    # adduct argument
    expect_equal(screenSuspects(fGroups, susps[name %in% c("TBA", "TPA")]),
                 screenSuspects(fGroups, susps[name %in% c("TBA", "TPA"), -"adduct"], adduct = "[M]+"))
})

TQFile <- file.path(getTestDataPath(), "GlobalResults-TASQ.csv")
TQRes <- fread(TQFile)
fGroupsTQ <- importFeatureGroupsBrukerTASQ(TQFile, getTestAnaInfo())
fGroupsTQ <- filter(fGroupsTQ, blankThreshold = 5, removeBlanks = TRUE)

test_that("TASQ import works", {
    expect_equal(unique(TQRes$Analyte), names(fGroupsTQ))
    expect_known_value(groups(fGroupsTQ), testFile("susp-tasq"))
    expect_known_show(fGroupsTQ, testFile("susp-tasq", text = TRUE))
})
