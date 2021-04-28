context("screening")

susps <- as.data.table(patRoonData::suspectsPos)

susps[, InChI := babelConvert(SMILES, "smi", "inchi")]
susps[, neutralMass := getNeutralMassFromSMILES(SMILES)]
susps[, formula := convertToFormulaBabel(SMILES, "smi")]
susps[, adduct := fifelse(name == "Aldicarb", "[M+Na]+", "[M+H]+")]

fGroups <- getTestFGroups(doFilter = FALSE, minFWHM = 1)
fGroupsScr <- doScreen(fGroups, susps, onlyHits = TRUE)
fGroupsScrNoRT <- doScreen(fGroups, susps[, -"rt"], onlyHits = TRUE)
getScrInfo <- function(susps, ...) screenInfo(doScreen(fGroups, susps, onlyHits = TRUE, ...))

scr <- getScrInfo(susps)
scrSMI <- getScrInfo(susps[, c("name", "rt", "SMILES", "adduct")])

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
    expect_setequal(scr$name, susps$name)
    expect_known_value(scr, testFile("screening"))
    expect_setequal(names(fGroupsScr), screenInfo(fGroupsScr)$group)

    # check suspects without retention
    expect_gte(nrow(getScrInfo(susps[, -3])), nrow(scr))
    
    # valid suspect names
    withr::with_options(list(patRoon.cache.mode = "none"), {
        expect_warning(doScreen(fGroups, data.table(name = "test.", SMILES = "C1=CC=C(C=C1)C(=O)O")))
        expect_error(doScreen(fGroups, data.table(name = "", SMILES = "C1=CC=C(C=C1)C(=O)O")))
    })
    
    # alternative ion mass calculations
    expect_equal_scr(scr, scrSMI, tolerance = 1E-3) # higher tolerance due to inaccurate mz column in patRoonData
    expect_equal_scr(scrSMI, getScrInfo(susps[, c("name", "rt", "InChI", "adduct")]))
    expect_equal_scr(scrSMI, getScrInfo(susps[, c("name", "rt", "neutralMass", "adduct")]))
    expect_equal_scr(scrSMI, getScrInfo(susps[, c("name", "rt", "formula", "adduct")]))

    # same, with missing data (having 2 options for ion mass calculation should be sufficient)
    expect_equal_scr(scrSMI, getScrInfo(suspsMissing[, c("name", "rt", "neutralMass", "formula", "adduct")]))
    expect_equal_scr(scrSMI, getScrInfo(suspsMissing[, c("name", "rt", "formula", "SMILES", "adduct")]))
    expect_equal_scr(scrSMI, getScrInfo(suspsMissing[, c("name", "rt", "SMILES", "InChI", "adduct")]))
    
    expect_warning(doScreen(fGroups, suspsMissingRow, skipInvalid = TRUE))
    expect_error(doScreen(fGroups, suspsMissingRow, skipInvalid = FALSE))
    
    # subsetting
    expect_length(fGroupsScr[, suspects = susps$name], length(fGroupsScr))
    expect_length(fGroupsScr[, suspects = susps$name[1:2]], 2)
    expect_length(fGroupsScr[, suspects = "doesnotexist"], 0)
    
    expect_equal(nrow(as.data.table(fGroupsScr, collapseSuspects = ",")), length(fGroupsScr))
    expect_gt(nrow(as.data.table(fGroupsScrNoRT, collapseSuspects = ",")),
              nrow(as.data.table(fGroupsScr, collapseSuspects = ",")))
    
    skip_if(testWithSets())

    expect_equal_scr(scrSMI, getScrInfo(suspsMissing[, c("name", "rt", "mz", "neutralMass", "adduct")]),
                     tolerance = 1E-3)
    
    # adduct argument
    # UNDONE: do something with selectIons?
    expect_equal_scr(getScrInfo(susps[name == "Aldicarb"]),
                     getScrInfo(susps[name == "Aldicarb", -"adduct"], adduct = "[M+Na]+"))
    
})

doAnnot <- function(...) annotateSuspects(..., logPath = NULL) # disable logging as it slow downs testing

# NOTE: try keep this in sync with MF tests for caching purposes
hasMF <- !is.null(getOption("patRoon.path.MetFragCL")) && nzchar(getOption("patRoon.path.MetFragCL"))
if (hasMF)
{
    plists <- generateMSPeakLists(fGroupsScr, "mzr")
    compsMF <- callMF(fGroupsScr, plists, db = file.path(getTestDataPath(), "test-mf-db-isomers.csv"))
    compsMFMoNa <- callMF(fGroupsScr, plists, scoreTypes = c("fragScore", "individualMoNAScore"),
                          db = file.path(getTestDataPath(), "test-mf-db-isomers.csv"))
    forms <- doGenForms(fGroupsScr, plists, "genform", calculateFeatures = FALSE)
    
    fGroupsAnnNothing <- doAnnot(fGroupsScr)
    fGroupsAnnMF <- doAnnot(fGroupsScr, MSPeakLists = plists, formulas = forms, compounds = compsMF)
    fGroupsAnnMFJ <- doAnnot(fGroupsScr, MSPeakLists = plists, formulas = forms, compounds = compsMF,
                                      simMSMSMethod = "jaccard")
    fGroupsAnnMoNA <- doAnnot(fGroupsScr, MSPeakLists = plists, formulas = forms, compounds = compsMFMoNa)
    fGroupsOnlyForms <- doAnnot(fGroupsScr, MSPeakLists = plists, formulas = forms)
    fGroupsAnnNoRT <- doAnnot(fGroupsScrNoRT, MSPeakLists = plists, formulas = forms, compounds = compsMFMoNa)
    
    fGroupsAnnFragNoRT <- doScreen(fGroupsScr, suspsFrag[, -"rt"], onlyHits = TRUE)
    fGroupsAnnFragNoRT <- doAnnot(fGroupsAnnFragNoRT, MSPeakLists = plists)
    fGroupsAnnFrag <- doScreen(fGroupsScr, suspsFrag, onlyHits = TRUE)
    fGroupsAnnFrag <- doAnnot(fGroupsAnnFrag, MSPeakLists = plists)
    
    idlFrag <- getWorkPath("fragtest.yml")
    genIDLevelRulesFile(idlFrag, exLevels = "3a|c")
    fGroupsAnnFragFormNoRT <- doScreen(fGroups, suspsFragForm[, -"rt"], onlyHits = TRUE)
    fGroupsAnnFragForm <- doAnnot(fGroupsAnnFragFormNoRT, MSPeakLists = plists, formulas = forms,
                                  compounds = compsMF, IDFile = idlFrag)
}

if (testWithSets())
{
    getAllSuspVals <- function(ann, col)
    {
        unlist(screenInfo(ann)[, getAllSuspSetCols(col, names(screenInfo(ann)), sets(ann)), with = FALSE])
    }
} else
{
    getAllSuspVals <- function(ann, col) screenInfo(ann)[[col]]
}

minIDLevel <- function(ann) min(numericIDLevel(getAllSuspVals(ann, "estIDLevel")), na.rm = TRUE)
maxIDLevel <- function(ann) max(numericIDLevel(getAllSuspVals(ann, "estIDLevel")), na.rm = TRUE)
getMinScrCol <- function(ann, col) min(getAllSuspVals(ann, col), na.rm = TRUE)
getMaxScrCol <- function(ann, col) max(getAllSuspVals(ann, col), na.rm = TRUE)

test_that("Suspect annotation works", {
    skip_if_not(hasMF)
    
    expect_known_value(screenInfo(fGroupsAnnMF), testFile("screen-ann-MF"))
    expect_known_value(screenInfo(fGroupsAnnMFJ), testFile("screen-ann-MF-J"))
    
    expect_equal(minIDLevel(fGroupsAnnNothing), 5)
    expect_equal(minIDLevel(fGroupsOnlyForms), 4)
    expect_equal(minIDLevel(fGroupsAnnMF), 3)
    expect_equal(minIDLevel(fGroupsAnnFragNoRT), 3)
    expect_equal(minIDLevel(fGroupsAnnFragForm), 3)
    expect_equal(minIDLevel(fGroupsAnnMoNA), 2)
    expect_equal(minIDLevel(fGroupsAnnFrag), 1)
    
    expect_equal(fGroupsAnnFrag, doAnnot(fGroupsAnnFrag, MSPeakLists = plists, checkFragments = "mz"))
    expect_false(isTRUE(all.equal(fGroupsAnnFrag, doAnnot(fGroupsAnnFrag, MSPeakLists = plists, checkFragments = "formula"))))
    expect_true(all(is.na(screenInfo(fGroupsAnnMF)$annSimBoth) |
                        screenInfo(fGroupsAnnMF)$annSimBoth >= pmax(screenInfo(fGroupsAnnMF)$annSimForm, screenInfo(fGroupsAnnMF)$annSimComp, na.rm = TRUE)))
})

# take fGroupsAnnNoRT: doesn't have rt in susp list, so has double hits
selectedHitsInt <- filter(fGroupsAnnNoRT, selectHitsBy = "intensity", onlyHits = TRUE)
selectedHitsLev <- filter(fGroupsAnnNoRT, selectHitsBy = "level", onlyHits = TRUE)
selectedFGroupsLev <- filter(fGroupsAnnNoRT, selectBestFGroups = TRUE, onlyHits = TRUE)

test_that("Screen filters", {
    expect_known_value(as.data.table(selectedHitsInt, collapseSuspects = NULL), testFile("screen-ann-sel-hits_int"))
    expect_known_value(as.data.table(selectedHitsLev, collapseSuspects = NULL), testFile("screen-ann-sel-hits_lev"))
    expect_known_value(as.data.table(selectedFGroupsLev, collapseSuspects = NULL), testFile("screen-ann-sel-groups"))
    
    # mean intensities should be increased when selecting a group for this suspect
    expect_lt(mean(as.matrix(as.data.table(fGroupsAnnNoRT, collapseSuspects = NULL)[susp_name == "DEET",
                                                                                    analyses(fGroupsAnnNoRT),
                                                                                    with = FALSE])),
              mean(as.matrix(as.data.table(selectedHitsInt, collapseSuspects = NULL)[susp_name == "DEET",
                                                                                     analyses(selectedHitsInt),
                                                                                     with = FALSE])))
    expect_gt(maxIDLevel(fGroupsAnnNoRT[, suspects = "DEET"]), maxIDLevel(selectedHitsLev[, suspects = "DEET"]))
    # UNDONE: these are not really good examples as the ID level is the same for all duplicates...
    # for now just verify that all groups are unique
    expect_equal(anyDuplicated(screenInfo(selectedFGroupsLev)$group), 0)
    
    expect_lt(length(selectedHitsInt), length(fGroupsAnnNoRT))
    expect_lt(length(selectedHitsLev), length(fGroupsAnnNoRT))
    expect_lt(length(selectedFGroupsLev), nrow(screenInfo(fGroupsAnnNoRT)))
    
    expect_equal(maxIDLevel(filter(fGroupsAnnNoRT, maxLevel = 3)), 3)
    expect_lte(getMaxScrCol(filter(fGroupsAnnNoRT, maxFormRank = 3), "formRank"), 3)
    expect_lte(getMaxScrCol(filter(fGroupsAnnNoRT, maxCompRank = 3), "compRank"), 3)
    expect_gte(getMinScrCol(filter(fGroupsAnnNoRT, minAnnSimForm = 0.9), "annSimForm"), 0.9)
    expect_gte(getMinScrCol(filter(fGroupsAnnNoRT, minAnnSimComp = 0.9), "annSimComp"), 0.9)
    expect_gte(getMinScrCol(filter(fGroupsAnnNoRT, minAnnSimBoth = 0.9), "annSimBoth"), 0.9)
    expect_equal(getMinScrCol(filter(fGroupsAnnFrag, absMinFragMatches = 1), "maxFragMatches"), 1)
    expect_equal(getMinScrCol(filter(fGroupsAnnFragForm, absMinFragMatches = 1), "maxFragMatches"), 1)
    expect_equal(getMinScrCol(filter(fGroupsAnnFrag, relMinFragMatches = 1), "maxFragMatches"), 1)
})

selectedNegHitsInt <- filter(fGroupsAnnNoRT, selectHitsBy = "intensity", onlyHits = FALSE, negate = TRUE)
selectedNegHitsLev <- filter(fGroupsAnnNoRT, selectHitsBy = "level", onlyHits = FALSE, negate = TRUE)
selectedNegFGroupsLev <- filter(fGroupsAnnNoRT, selectBestFGroups = TRUE, onlyHits = FALSE, negate = TRUE)

test_that("Negated screen filters", {
    expect_known_value(as.data.table(selectedNegHitsInt, collapseSuspects = NULL), testFile("screen-ann-sel-neg-hits_int"))
    expect_known_value(as.data.table(selectedNegHitsLev, collapseSuspects = NULL), testFile("screen-ann-sel-neg-hits_lev"))
    expect_known_value(as.data.table(selectedNegFGroupsLev, collapseSuspects = NULL), testFile("screen-ann-sel-neg-groups"))
    
    # as above, but opposite
    expect_gt(mean(as.matrix(as.data.table(fGroupsAnnNoRT, collapseSuspects = NULL)[susp_name == "DEET",
                                                                                    analyses(fGroupsAnnNoRT),
                                                                                    with = FALSE])),
              mean(as.matrix(as.data.table(selectedNegHitsInt, collapseSuspects = NULL)[susp_name == "DEET",
                                                                                        analyses(selectedNegHitsInt),
                                                                                        with = FALSE])))
    expect_lt(minIDLevel(fGroupsAnnNoRT[, suspects = "DEET"]), minIDLevel(selectedNegHitsLev[, suspects = "DEET"]))
    # UNDONE: these are not really good examples as the ID level is the same for all duplicates...
    # for now just verify that all groups are unique
    expect_equal(anyDuplicated(screenInfo(selectedNegFGroupsLev)$group), 0)
    
    expect_length(selectedNegHitsInt, length(selectedHitsInt))
    expect_length(selectedNegHitsLev, length(selectedHitsLev))
    expect_equal(nrow(screenInfo(selectedNegFGroupsLev)), nrow(screenInfo(selectedFGroupsLev)))
    
    expect_gt(maxIDLevel(filter(fGroupsAnnNoRT, maxLevel = 3, negate = TRUE)), 3)
    expect_gt(getMaxScrCol(filter(fGroupsAnnNoRT, maxFormRank = 3, negate = TRUE), "formRank"), 3)
    expect_gt(getMaxScrCol(filter(fGroupsAnnNoRT, maxCompRank = 2, negate = TRUE), "compRank"), 2)
    expect_lt(getMinScrCol(filter(fGroupsAnnNoRT, minAnnSimForm = 0.9, negate = TRUE), "annSimForm"), 0.9)
    expect_lt(getMinScrCol(filter(fGroupsAnnNoRT, minAnnSimComp = 0.9, negate = TRUE), "annSimComp"), 0.9)
    expect_lt(getMinScrCol(filter(fGroupsAnnNoRT, minAnnSimBoth = 0.9, negate = TRUE), "annSimBoth"), 0.9)
    expect_true(all(is.na(screenInfo(filter(fGroupsAnnFrag, absMinFragMatches = 1, negate = TRUE))$maxFragMatches)))
    expect_true(all(is.na(screenInfo(filter(fGroupsAnnFragForm, absMinFragMatches = 1, negate = TRUE))$maxFragMatches)))
    expect_true(all(is.na(screenInfo(filter(fGroupsAnnFrag, relMinFragMatches = 1, negate = TRUE))$maxFragMatches)))
})

fGroupsEmpty <- groupFeatures(findFeatures(getTestAnaInfo(), "openms", noiseThrInt = 1E9), "openms")
suspsEmpty <- data.table(name = "doesnotexist", SMILES = "C")
fGroupsScrEmpty <- doScreen(fGroups, suspsEmpty)

if (hasMF)
    fGroupsScrAnnEmpty <- doAnnot(fGroupsScrEmpty, MSPeakLists = plists, formulas = forms, compounds = compsMF)

test_that("Empty objects", {
    expect_length(doScreen(fGroupsEmpty, suspsEmpty), 0)
    expect_length(fGroupsScrEmpty, length(fGroups))
    expect_equal(nrow(as.data.table(fGroupsScrEmpty, onlyHits = TRUE)), 0)
    
    expect_length(filter(fGroupsScrEmpty, onlyHits = TRUE), 0)
    expect_length(filter(fGroupsScrEmpty, onlyHits = TRUE, negate = TRUE), length(fGroups))

    skip_if_not(hasMF)
    expect_length(fGroupsScrAnnEmpty, length(fGroups))
    expect_equal(nrow(as.data.table(fGroupsScrAnnEmpty, onlyHits = TRUE)), 0)
    expect_length(filter(fGroupsScrAnnEmpty, selectHitsBy = "intensity", onlyHits = TRUE), 0)
    expect_length(filter(fGroupsScrAnnEmpty, selectHitsBy = "level", onlyHits = TRUE), 0)
    expect_length(filter(fGroupsScrAnnEmpty, selectBestFGroups = TRUE, onlyHits = TRUE), 0)
    expect_length(filter(fGroupsScrAnnEmpty, minAnnSimForm = 0.0, onlyHits = TRUE), 0)
})

csvSuspCols <- c("susp_name", "susp_compRank", "susp_annSimBoth", "susp_estIDLevel")
if (testWithSets())
    csvSuspCols <- getAllSuspSetCols(csvSuspCols, names(screenInfo(fGroupsAnnNoRT)), sets(fGroupsAnnNoRT))

test_that("reporting works", {
    skip_if_not(hasMF)
    
    expect_file(reportCSV(fGroupsAnnNoRT, getWorkPath()),
                getWorkPath(paste0(class(fGroupsAnnNoRT), ".csv")))
    checkmate::expect_names(names(fread(getWorkPath(paste0(class(fGroupsAnnNoRT), ".csv")))),
                            must.include = csvSuspCols)
    expect_file(reportPDF(fGroupsAnnNoRT, getWorkPath()), getWorkPath(paste0(class(fGroupsAnnNoRT), ".pdf")))
    expect_reportHTML(makeReportHTML(fGroupsAnnNoRT))

    expect_error(reportCSV(fGroupsScrAnnEmpty[, 1:10], getWorkPath()), NA)
    expect_error(reportPDF(fGroupsScrAnnEmpty[, 1:10], getWorkPath()), NA)
    expect_error(makeReportHTML(fGroupsScrAnnEmpty[, 1:10]), NA)
})

test_that("sets functionality", {
    skip_if_not(testWithSets())
    
    # some tests from feature groups to ensure proper subsetting/unsetting
    expect_equal(analysisInfo(unset(fGroupsScr, "positive")), getTestAnaInfoPos())
    expect_equal(analysisInfo(fGroupsScr[, sets = "positive"])[, 1:4], getTestAnaInfoPos())
    expect_setequal(unique(annotations(fGroupsScr)$adduct), c("[M+H]+", "[M-H]-"))
    expect_equal(fGroupsScr, fGroupsScr[, sets = sets(fGroupsScr)])
    expect_length(fGroupsScr[, sets = character()], 0)
    expect_equal(sets(filter(fGroupsScr, sets = "positive", negate = TRUE)), "negative")
})

getTQAnalytes <- function(path)
{
    ret <- fread(path)[["Analyte Name"]]
    return(ret[nzchar(ret)]) # omit first empty line
}

TQFile <- file.path(getTestDataPath(), "GlobalResults-TASQ-pos.csv")
fGroupsTQ <- importFeatureGroupsBrukerTASQ(TQFile, getTestAnaInfo(), clusterRTWindow = 14) # HACK: increase RT window a bit to avoid split groups
TQFileNoRT <- file.path(getTestDataPath(), "GlobalResults-TASQ_noRT-pos.csv")
fGroupsTQNoRT <- importFeatureGroupsBrukerTASQ(TQFileNoRT, getTestAnaInfo())
fGroupsTQNoRT <- filter(fGroupsTQNoRT, blankThreshold = 5, removeBlanks = TRUE)

test_that("TASQ import works", {
    expect_equal(unique(getTQAnalytes(TQFile)), names(fGroupsTQ))
    expect_known_value(groupTable(fGroupsTQ), testFile("susp-tasq"))
    expect_known_value(groupInfo(fGroupsTQ), testFile("susp-tasq-gi"))
    expect_known_show(fGroupsTQ, testFile("susp-tasq", text = TRUE))
    
    expect_lt(length(unique(getTQAnalytes(TQFileNoRT))), length(fGroupsTQNoRT))
    expect_known_value(groupTable(fGroupsTQNoRT), testFile("susp-tasq-nort"))
    expect_known_value(groupInfo(fGroupsTQNoRT), testFile("susp-tasq-nort-gi"))
    expect_known_show(fGroupsTQNoRT, testFile("susp-tasq-nort", text = TRUE))
})
