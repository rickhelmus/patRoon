context("screening")

susps <- as.data.table(patRoonData::targets)

susps[, InChI := babelConvert(SMILES, "smi", "inchi")]
susps[, neutralMass := getNeutralMassFromSMILES(SMILES)]
susps[, formula := convertToFormulaBabel(SMILES, "smi")]
susps[, adduct := "[M+H]+"]
susps[name %in% c("TBA", "TPA"), adduct := "[M]+"]

fGroups <- getTestFGroups(getTestAnaInfo())

getScrInfo <- function(susps, ...) screenInfo(screenSuspects(fGroups, susps, onlyHits = TRUE, ...))

scr <- getScrInfo(susps)
scrSMI <- getScrInfo(susps[, c("name", "rt", "adduct", "SMILES")])

suspsMissing <- copy(susps)
suspsMissing[1, mz := NA]
suspsMissing[2, neutralMass := NA]
suspsMissing[3, formula := NA_character_]
suspsMissing[4, SMILES := ""]
suspsMissing[5, InChI := ""]
suspsMissingRow <- copy(susps)
suspsMissingRow[2, c("mz", "neutralMass", "formula", "SMILES", "InChI") := NA]

suspsFrag <- copy(susps)
suspsFrag[name == "1H-benzotriazole", fragments_mz := "92.0495"] # UNDONE: add qualifiers to patRoonData?
suspsFragForm <- copy(susps)
suspsFragForm[name == "1H-benzotriazole", fragments_formula := "C6H6N"] # UNDONE: add qualifiers to patRoonData?

test_that("suspect screening is OK", {
    expect_equal(nrow(scr), nrow(susps))
    expect_known_value(scr, testFile("screening"))
    expect_length(screenSuspects(fGroups, susps, onlyHits = TRUE), nrow(susps))

    # check suspects without retention
    expect_gte(nrow(getScrInfo(susps[, -3])), nrow(scr))
    
    # alternative ion mass calculations
    expect_equal_scr(scr, scrSMI, tolerance = 1E-3) # higher tolerance due to inaccurate mz column in patRoonData::targets
    expect_equal_scr(scrSMI, getScrInfo(susps[, c("name", "rt", "adduct", "InChI")]))
    expect_equal_scr(scrSMI, getScrInfo(susps[, c("name", "rt", "adduct", "neutralMass")]))
    expect_equal_scr(scrSMI, getScrInfo(susps[, c("name", "rt", "adduct", "formula")]))

    # same, with missing data (having 2 options for ion mass calculation should be sufficient)
    expect_equal_scr(scrSMI, getScrInfo(suspsMissing[, c("name", "rt", "adduct", "mz", "neutralMass")]),
                 tolerance = 1E-3)
    expect_equal_scr(scrSMI, getScrInfo(suspsMissing[, c("name", "rt", "adduct", "neutralMass", "formula")]))
    expect_equal_scr(scrSMI, getScrInfo(suspsMissing[, c("name", "rt", "adduct", "formula", "SMILES")]))
    expect_equal_scr(scrSMI, getScrInfo(suspsMissing[, c("name", "rt", "adduct", "SMILES", "InChI")]))
    
    # adduct argument
    expect_equal_scr(getScrInfo(susps[name %in% c("TBA", "TPA")]),
                     getScrInfo(susps[name %in% c("TBA", "TPA"), -"adduct"], adduct = "[M]+"))

    expect_warning(screenSuspects(fGroups, suspsMissingRow, skipInvalid = TRUE))
    expect_error(screenSuspects(fGroups, suspsMissingRow, skipInvalid = FALSE))
})

# NOTE: try keep this in sync with MF tests for caching purposes
hasMF <- !is.null(getOption("patRoon.path.MetFragCL")) && nzchar(getOption("patRoon.path.MetFragCL"))
if (hasMF)
{
    fGroupsForAnn <- getCompFGroups()
    plists <- generateMSPeakLists(fGroupsForAnn, "mzr")
    compsMF <- callMF(fGroupsForAnn, plists)
    compsMFMoNa <- callMF(fGroupsForAnn, plists, scoreTypes = c("fragScore", "individualMoNAScore"))
    forms <- generateFormulas(fGroupsForAnn, "genform", plists)
    
    fGroupsAnn <- annotateSuspects(fGroupsForAnn)
    fGroupsAnnMF <- annotateSuspects(fGroupsForAnn, MSPeakLists = plists, formulas = forms, compounds = compsMF)
    fGroupsAnnMoNA <- annotateSuspects(fGroupsAnn, MSPeakLists = plists, formulas = forms, compounds = compsMFMoNa)
    fGroupsOnlyForms <- annotateSuspects(fGroupsAnn, MSPeakLists = plists, formulas = forms)
    
    fGroupsAnnFragNoRT <- screenSuspects(fGroupsForAnn, suspsFrag[, -"rt"])
    fGroupsAnnFragNoRT <- annotateSuspects(fGroupsAnnFragNoRT, MSPeakLists = plists)
    fGroupsAnnFrag <- screenSuspects(fGroupsForAnn, suspsFrag)
    fGroupsAnnFrag <- annotateSuspects(fGroupsAnnFrag, MSPeakLists = plists)
    
    idlFrag <- getWorkPath("fragtest.yml")
    genIDLevelRulesFile(idlFrag, exLevels = "3a|c")
    fGroupsAnnFragForm <- screenSuspects(fGroupsForAnn, suspsFragForm[, -"rt"])
    fGroupsAnnFragForm <- annotateSuspects(fGroupsAnnFragForm, MSPeakLists = plists, compounds = compsMF,
                                           IDFile = idlFrag)
}

minIDLevel <- function(ann) min(numericIDLevel(screenInfo(ann)$estIDLevel))
test_that("Suspect annotation works", {
    skip_if_not(hasMF)
    
    expect_equal(minIDLevel(fGroupsAnn), 5)
    expect_equal(minIDLevel(fGroupsOnlyForms), 4)
    expect_equal(minIDLevel(fGroupsAnnMF), 3)
    expect_equal(minIDLevel(fGroupsAnnFragNoRT), 3)
    expect_equal(minIDLevel(fGroupsAnnFragForm), 3)
    expect_equal(minIDLevel(fGroupsAnnMoNA), 2)
    expect_equal(minIDLevel(fGroupsAnnFrag), 1)
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
    expect_known_value(groupTable(fGroupsTQ), testFile("susp-tasq"))
    expect_known_value(groupInfo(fGroupsTQ), testFile("susp-tasq-gi"))
    expect_known_show(fGroupsTQ, testFile("susp-tasq", text = TRUE))
    
    expect_lt(length(unique(TQResNoRT$Analyte)), length(fGroupsTQNoRT))
    expect_known_value(groupTable(fGroupsTQNoRT), testFile("susp-tasq-nort"))
    expect_known_value(groupInfo(fGroupsTQNoRT), testFile("susp-tasq-nort-gi"))
    expect_known_show(fGroupsTQNoRT, testFile("susp-tasq-nort", text = TRUE))
})
