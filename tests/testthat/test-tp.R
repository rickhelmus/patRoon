context("TPs")

fGroups <- getCompFGroups()

suspL <- patRoonData::suspectsPos[patRoonData::suspectsPos$name %in% screenInfo(fGroups)$name, ]

TPsLogic <- doGenLogicTPs(fGroups)
TPsLogicCustom <- doGenLogicTPs(fGroups, transformations = data.table(transformation = "test", add = "C",
                                                                      sub = "", retDir = 1))
TPsLibPC <- generateTPs("library")
TPsLibSusp <- generateTPs("library", suspL)
TPsLibScr <- generateTPs("library", fGroups)
TPsLibCustom <- generateTPs("library",
                            TPLibrary = data.table(parent_name = "Diuron", parent_SMILES = "CN(C)C(O)=NC1C=C(Cl)C(Cl)=CC=1",
                                                   TP_name = "Linuron", TP_SMILES = "CN(C(=O)NC1=CC(=C(C=C1)Cl)Cl)OC"))
TPsBTSusp <- generateTPs("biotransformer", suspL)
TPsBTScr <- generateTPs("biotransformer", fGroups)
TPsBTSuspMore <- generateTPs("biotransformer", patRoonData::suspectsPos[1:25, ]) # for filter tests

fGroupsEmpty <- getEmptyTestFGroups()
fGroupsScrEmpty <- doScreen(fGroupsEmpty, data.table(name = "doesnotexist", SMILES = "C", mz = 12))
TPsLogicEmpty <- doGenLogicTPs(fGroupsEmpty)
TPsLibEmpty <- generateTPs("library", fGroupsScrEmpty)
TPsBTEmpty <- generateTPs("biotransformer", fGroupsScrEmpty)

doMetFrag <- !is.null(getOption("patRoon.path.MetFragCL")) && nzchar(getOption("patRoon.path.MetFragCL"))

if (doMetFrag)
{
    plists <- generateMSPeakLists(fGroups, "mzr")
    plistsEmpty <- plists[FALSE, reAverage = TRUE]
    plistsEmptyMS <- removeMSPlists(plists, "MS")
    
    compsMF <- callMF(fGroups, plists)
    TPsLibComp <- generateTPs("library", compsMF)
    TPsBTComp <- generateTPs("biotransformer", compsMF)
}

test_that("verify TP generation", {
    expect_known_value(TPsLogic, testFile("tp-logic"))
    expect_known_value(TPsLibPC, testFile("tp-lib_pc"))
    expect_known_value(TPsLibSusp, testFile("tp-lib_susp"))
    expect_known_value(TPsLibScr, testFile("tp-lib_scr"))
    expect_known_value(TPsBTSusp, testFile("tp-bt_susp"))
    expect_known_value(TPsBTScr, testFile("tp-bt_scr"))

    expect_known_show(TPsLogic, testFile("tp-logic", text = TRUE))
    expect_known_show(TPsLibPC, testFile("tp-lib_pc", text = TRUE))
    expect_known_show(TPsLibSusp, testFile("tp-lib_susp", text = TRUE))
    expect_known_show(TPsLibScr, testFile("tp-lib_scr", text = TRUE))
    expect_known_show(TPsBTSusp, testFile("tp-bt_susp", text = TRUE))
    expect_known_show(TPsBTScr, testFile("tp-bt_scr", text = TRUE))

    expect_setequal(parents(TPsLogic)$name, names(fGroups))
    checkmate::expect_names(names(parents(TPsLogic)), permutation.of = c("name", "rt", "neutralMass"))
    expect_true(all(as.data.table(doGenLogicTPs(fGroups, minMass = 100))$neutralMass >= 100))
    expect_length(TPsLogicCustom, length(fGroups))

    expect_setequal(parents(TPsLibPC)$name, patRoon:::PubChemTransformations$parent_name)
    checkmate::expect_names(parents(TPsLibSusp)$name, subset.of = suspL$name)
    checkmate::expect_names(parents(TPsLibScr)$name, subset.of = screenInfo(fGroups)$name)
    checkmate::expect_names(names(parents(TPsLibPC)), must.include = c("name", "SMILES", "InChI", "InChIKey",
                                                                       "formula", "neutralMass"))
    expect_equal(parents(TPsLibCustom)$name, "Diuron")
    expect_match(TPsLibCustom[[1]]$name, "Linuron")

    checkmate::expect_names(parents(TPsBTSusp)$name, subset.of = suspL$name)
    checkmate::expect_names(parents(TPsBTScr)$name, subset.of = screenInfo(fGroups)$name)
    checkmate::expect_names(names(parents(TPsBTSusp)), must.include = c("name", "SMILES", "InChI", "InChIKey",
                                                                        "formula", "neutralMass"))
    
    expect_length(TPsLogicEmpty, 0)
    expect_length(TPsLibEmpty, 0)
    expect_length(TPsBTEmpty, 0)
    
    skip_if_not(doMetFrag)

    expect_known_value(TPsLibComp, testFile("tp-lib_comp"))
    expect_known_value(TPsBTComp, testFile("tp-bt_comp"))
    expect_known_show(TPsLibComp, testFile("tp-lib_comp", text = TRUE))
    expect_known_show(TPsBTComp, testFile("tp-bt_comp", text = TRUE))

    checkmate::expect_names(parents(TPsLibComp)$name, subset.of = as.data.table(compsMF)$identifier)
    checkmate::expect_names(parents(TPsBTComp)$name, subset.of = as.data.table(compsMF)$identifier)
})

assertSusp <- function(...) patRoon:::assertSuspectList(..., needsAdduct = FALSE, skipInvalid = FALSE)
testMFDB <- function(...)
{
    outf <- tempfile(fileext = ".csv")
    convertToMFDB(..., out = outf, includeParents = FALSE); db <- fread(outf)
    cols <- c("Identifier", "MolecularFormula", "MonoisotopicMass",
              "SMILES", "InChI", "InChIKey", "InChIKey1", "ALogP", "LogP", "parent", "transformation",
              "enzyme", "evidencedoi")
    
    checkmate::expect_data_table(MFDBNP, any.missing = FALSE, nrows = length(TPsLibScr))
    checkmate::expect_names(names(MFDBNP), subset.of = cols)
    
    convertToMFDB(..., out = outf, includeParents = TRUE); db <- fread(outf)
    checkmate::expect_data_table(MFDBP, nrows = length(TPsLibScr) + nrow(parents(TPsLibScr)))
    checkmate::expect_names(names(MFDBP), subset.of = cols)
}

getMFDB <- function(...) { outf <- tempfile(fileext = ".csv"); convertToMFDB(..., out = outf); fread(outf) }
MFDBNP <- getMFDB(TPsLibScr, includeParents = FALSE)
MFDBP <- getMFDB(TPsLibScr, includeParents = TRUE)
TPsBTSuspFEF <- filter(TPsBTSuspMore, removeEqualFormulas = TRUE)
TPsBTSuspFEFN <- filter(TPsBTSuspMore, removeEqualFormulas = TRUE, negate = TRUE)
pnames <- parents(TPsLogic)$name
test_that("basic usage", {
    expect_length(TPsLogic["nope"], 0)
    expect_equivalent(names(TPsLogic[1:2]), pnames[1:2])
    expect_equivalent(names(TPsLogic[pnames[2:3]]), pnames[2:3])
    expect_equivalent(names(TPsLogic[c(TRUE, FALSE)]), pnames[c(TRUE, FALSE)])
    expect_equal(length(TPsLogic[FALSE]), 0)
    expect_length(TPsLogicEmpty[1:5], 0)
    
    expect_equivalent(TPsLogic[[2]], products(TPsLogic)[[2]])
    expect_equivalent(TPsLogic[[names(TPsLogic)[2]]], products(TPsLogic)[[2]])
    expect_equivalent(callDollar(TPsLogic, names(TPsLogic)[2]), TPsLogic[[2]])
    
    expect_equal(nrow(as.data.table(TPsLogic)), length(TPsLogic))
    
    expect_error(assertSusp(convertToSuspects(TPsLogic)), NA)
    expect_error(assertSusp(convertToSuspects(TPsLibScr, includeParents = TRUE)), NA)
    expect_error(assertSusp(convertToSuspects(TPsBTSusp, includeParents = TRUE)), NA)
    expect_error(convertToSuspects(TPsLogicEmpty), "create")
    
    testMFDB(TPsLibScr)
    testMFDB(TPsBTSusp)
    expect_error(convertToMFDB(TPsBTEmpty, tempfile(fileext = ".csv")), "create")
    
    expect_false(any(mapply(parents(TPsBTSuspFEF)$formula, products(TPsBTSuspFEF), FUN = function(f, p) any(f %in% p$formula))))
    expect_true(all(mapply(parents(TPsBTSuspFEFN)$formula, products(TPsBTSuspFEFN), FUN = function(f, p) all(f %in% p$formula))))
    expect_gte(min(as.data.table(filter(TPsBTSuspMore, minSimilarity = 0.5))$similarity), 0.5)
    expect_lt(max(as.data.table(filter(TPsBTSuspMore, minSimilarity = 0.5, negate = TRUE))$similarity), 0.5)
    
    skip_if_not(doMetFrag)
    
    expect_error(assertSusp(convertToSuspects(TPsLibComp, includeParents = TRUE)), NA)
    expect_error(assertSusp(convertToSuspects(TPsBTComp, includeParents = TRUE)), NA)
    testMFDB(TPsLibComp)
    testMFDB(TPsBTComp)
})

fGroupsMore <- getTestFGroups(getTestAnaInfoAnn())
componTPsNone <- generateComponents(fGroupsMore[, 1:50], "tp", TPs = NULL)
componTPsNoneTPDiff <- generateComponents(fGroupsMore[, 1:25], "tp", fGroupsMore[, 26:50], TPs = NULL)

componTPsLib <- generateComponents(doScreen(fGroupsMore, convertToSuspects(TPsLibSusp, includeParents = TRUE)), "tp",
                                   TPs = TPsLibSusp)

TPsLogicMore <- doGenLogicTPs(fGroupsMore[, 1:50])
fGroupsTPsLogic <- doScreen(fGroupsMore, convertToSuspects(TPsLogicMore), onlyHits = TRUE)
componTPsLogic <- generateComponents(fGroupsTPsLogic, "tp", TPs = TPsLogicMore)

fGroupsMoreScr <- doScreen(fGroupsMore, convertToSuspects(TPsBTSusp, includeParents = TRUE), onlyHits = TRUE)
componTPsBT <- generateComponents(fGroupsMoreScr, "tp", TPs = TPsBTSusp)

if (doMetFrag)
{
    fGroupsAnn <- doScreen(fGroupsMore, convertToSuspects(TPsLibSusp, includeParents = TRUE), onlyHits = TRUE)
    plistsAnn <- generateMSPeakLists(fGroupsAnn, "mzr")
    formsAnn <- doGenForms(fGroupsAnn, plistsAnn, "genform", elements = "CHNOPSClF", calculateFeatures = FALSE)
    TPDB <- tempfile(fileext = ".csv")
    convertToMFDB(TPsLibSusp, TPDB, includeParents = TRUE)
    compsAnn <- doGenComps(fGroupsAnn, plistsAnn, "metfrag", database = "csv", extraOpts = list(LocalDatabasePath = TPDB))
    componTPsAnn <- generateComponents(fGroupsAnn, "tp", TPs = TPsLibSusp, MSPeakLists = plistsAnn, formulas = formsAnn,
                                       compounds = compsAnn)
    componTPsAnnPL <- generateComponents(fGroupsAnn, "tp", TPs = TPsLibSusp, MSPeakLists = plistsAnn)
}

test_that("TP componentization", {
    expect_known_value(componTPsNone, testFile("tp-compon-none"))
    expect_known_value(componTPsLogic, testFile("tp-compon-logic"))
    expect_known_value(componTPsLib, testFile("tp-compon-lib"))
    expect_known_value(componTPsBT, testFile("tp-compon-bt"))
    
    expect_known_show(componTPsNone, testFile("tp-compon-none", text = TRUE))
    expect_known_show(componTPsLogic, testFile("tp-compon-logic", text = TRUE))
    expect_known_show(componTPsLib, testFile("tp-compon-lib", text = TRUE))
    expect_known_show(componTPsBT, testFile("tp-compon-bt", text = TRUE))
    
    expect_equal(nrow(as.data.table(componTPsNone)), length(groupNames(componTPsNone))^2)
    expect_length(intersect(componentInfo(componTPsNoneTPDiff)$parent_group, groupNames(componTPsNoneTPDiff)), 0)
    expect_equal(componTPsNoneTPDiff, generateComponents(fGroupsMore[, 1:25], "tp", fGroupsMore[, 26:50],
                                                         ignoreParents = TRUE, TPs = NULL))
    expect_warning(generateComponents(fGroupsMore[, 1:25], "tp", ignoreParents = TRUE, TPs = NULL))
    expect_warning(generateComponents(fGroupsMore[, 1:25], "tp", fGroupsMore[, 10:15], ignoreParents = TRUE, TPs = NULL))
    expect_setequal(groupNames(generateComponents(fGroupsMore[, 1:25], "tp", fGroupsMore[, 20:30], ignoreParents = TRUE,
                                                  TPs = NULL)), names(fGroupsMore)[26:30])
    
    expect_setequal(as.data.table(componTPsLogic)$retDir, c(-1, 0, 1))
    expect_true(all(as.data.table(generateComponents(doScreen(fGroupsMore,
                                                              convertToSuspects(TPsLogicMore, includeParents = TRUE)),
                                                     "tp", TPs = TPsLogicMore, minRTDiff = 1E5))$retDir == 0))
    
    expect_length(generateComponents(fGroups, "tp", TPs = TPsLogicEmpty), 0)
    expect_length(generateComponents(fGroupsScrEmpty, "tp", TPs = TPsLogicEmpty), 0)
    expect_length(generateComponents(fGroupsEmpty, "tp", TPs = NULL), 0)
    
    skip_if_not(doMetFrag)
    
    expect_known_value(componTPsAnn, testFile("tp-compon-ann"))
    expect_known_value(componTPsAnnPL, testFile("tp-compon-ann_pl"))
    
    expect_known_show(componTPsAnn, testFile("tp-compon-ann", text = TRUE))
    expect_known_show(componTPsAnnPL, testFile("tp-compon-ann_pl", text = TRUE))
    
    checkmate::expect_names(names(as.data.table(componTPsAnn)), must.include = c("specSimilarity", "specSimilarityPrec",
                                                                                 "specSimilarityBoth",
                                                                                 "fragmentMatches",
                                                                                 "neutralLossMatches"))
    checkmate::expect_names(names(as.data.table(componTPsAnnPL)), must.include = c("specSimilarity", "specSimilarityPrec",
                                                                                   "specSimilarityBoth"))
    
    expect_error(generateComponents(fGroupsMore[, 1:50], "tp", TPs = NULL, MSPeakLists = plistsEmpty), NA)
    expect_error(generateComponents(fGroupsMore[, 1:50], "tp", TPs = NULL, MSPeakLists = plistsEmptyMS), NA)
})

getMinCompTblVal <- function(x, field) min(as.data.table(x)[[field]], na.rm = TRUE)
getMaxCompTblVal <- function(x, field) max(as.data.table(x)[[field]], na.rm = TRUE)

componTPsRetF <- filter(componTPsLogic, retDirMatch = TRUE)
componTPsRetFN <- filter(componTPsLogic, retDirMatch = TRUE, negate = TRUE)

if (doMetFrag)
{
    # HACK: add unlikely element so we can easily test formula filter below
    componTPsLogicMod <- componTPsLogic
    componTPsLogicMod@components[[1]] <- copy(componTPsLogicMod@components[[1]])
    componTPsLogicMod@components[[1]][, trans_add := "Cl"]
    
    plistsLogicAnn <- generateMSPeakLists(fGroupsTPsLogic, "mzr")
    formsLogicAnn <- doGenForms(fGroupsTPsLogic, plistsLogicAnn, "genform", calculateFeatures = FALSE)
}

test_that("TP component usage", {
    expect_HTML(plotGraph(componTPsLogic, onlyLinked = FALSE))
    
    expect_equal(as.data.table(componTPsLogic)[TP_retDir == retDir | TP_retDir == 0 | retDir == 0, -"size"],
                 as.data.table(componTPsRetF)[, -"size"])
    expect_equal(as.data.table(componTPsLogic)[TP_retDir != retDir | TP_retDir == 0 | retDir == 0, -"size"],
                 as.data.table(componTPsRetFN)[, -"size"])
    
    expect_reportHTML(makeReportHTML(fGroupsMoreScr, components = componTPsBT))
    
    skip_if_not(doMetFrag)
    
    expect_lt(nrow(as.data.table(filter(componTPsLogicMod, formulas = formsLogicAnn))),
              nrow(as.data.table(componTPsLogicMod)))
    
    expect_gte(getMinCompTblVal(filter(componTPsAnn, minSpecSim = 0.2), "specSimilarity"), 0.2)
    expect_gte(getMinCompTblVal(filter(componTPsAnn, minSpecSimPrec = 0.2), "specSimilarityPrec"), 0.2)
    expect_gte(getMinCompTblVal(filter(componTPsAnn, minSpecSimBoth = 0.2), "specSimilarityBoth"), 0.2)
    expect_gte(getMinCompTblVal(filter(componTPsAnn, minFragMatches = 2), "fragmentMatches"), 2)
    expect_gte(getMinCompTblVal(filter(componTPsAnn, minNLMatches = 2), "neutralLossMatches"), 2)
    
    expect_lt(getMaxCompTblVal(filter(componTPsAnn, minSpecSim = 0.2, negate = TRUE), "specSimilarity"), 0.2)
    expect_lt(getMaxCompTblVal(filter(componTPsAnn, minSpecSimPrec = 0.2, negate = TRUE), "specSimilarityPrec"), 0.2)
    expect_lt(getMaxCompTblVal(filter(componTPsAnn, minSpecSimBoth = 0.3, negate = TRUE), "specSimilarityBoth"), 0.3)
    expect_lt(getMaxCompTblVal(filter(componTPsAnn, minFragMatches = 2, negate = TRUE), "fragmentMatches"), 2)
    expect_lt(getMaxCompTblVal(filter(componTPsAnn, minNLMatches = 4, negate = TRUE), "neutralLossMatches"), 4)
})
