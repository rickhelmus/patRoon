context("screening")

susps <- as.data.table(patRoonData::targets)

# UNDONE: update patRoonData at some point
susps[, SMILES := c("[nH]1nnc2ccccc12", "CNc1ccccn1", "Oc1ccc2ccccc2n1", "Oc1ccnc2ccccc12",
                    "Oc1ccc2ncccc2c1", "Oc1cccc2cccnc12", "Cc1ccc2[nH]nnc2c1", "CN1N(C(=O)C=C1C)c2ccccc2",
                    "CCNc1nc(Cl)nc(NC(C)C)n1","C1=CC(=C(C(=C1)Cl)C(=O)N)Cl", "Cn1cnc2N(C)C(=O)N(C)C(=O)c12", "NC(=O)N1c2ccccc2C=Cc3ccccc13",
                    "NC1=C(Cl)C(=O)N(N=C1)c2ccccc2", "CCN(CC)C(=O)c1cccc(C)c1",
                    "COCCOCCOC", "CN(C)C(=O)Nc1ccc(Cl)c(Cl)c1", "CC(=O)Nc1ccc(O)cc1",
                    "NC(=O)Nc1ccccc1", "Cc1cc(C)nc(N[S](=O)(=O)c2ccc(N)cc2)n1",
                    "Cc1onc(N[S](=O)(=O)c2ccc(N)cc2)c1", "CCCC[N+](CCCC)(CCCC)CCCC", "CCO[P](=O)(OCC)OCC",
                    "CCCC[NH+](CCCC)CCCC")]
susps[, InChI := babelConvert(SMILES, "smi", "inchi")]
susps[, neutralMass := getNeutralMassFromSMILES(SMILES)]
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
