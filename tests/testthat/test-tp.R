# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

context("TPs")

fGroups <- getCompFGroups()

suspL <- patRoonData::suspectsPos[patRoonData::suspectsPos$name %in% screenInfo(fGroups)$name, ]

TPsLogic <- doGenLogicTPs(fGroups)
TPsLogicCustom <- doGenLogicTPs(fGroups, transformations = data.table(transformation = "test", add = "C",
                                                                      sub = "", retDir = 1))
TPsLibPC <- generateTPs("library")
TPsLibPCGen2 <- generateTPs("library", generations = 2)
TPsLibSusp <- generateTPs("library", suspL)
TPsLibScr <- generateTPs("library", fGroups)
TPCustLib <- data.table(parent_name = "Diuron", parent_SMILES = "CN(C)C(O)=NC1C=C(Cl)C(Cl)=CC=1",
                        TP_name = "Linuron", TP_SMILES = "CN(C(=O)NC1=CC(=C(C=C1)Cl)Cl)OC")
TPsLibCustom <- generateTPs("library", TPLibrary = TPCustLib)

TPCustFormLib <- data.table(parent_name = "Diuron", parent_formula = "C9H10Cl2N2O",
                            TP_name = "Monuron", TP_formula = "C9H11ClN2O")
TPsLibForm <- generateTPs("library_formula", TPLibrary = TPCustFormLib)
TPsLibFormGen2 <- generateTPs("library_formula", TPLibrary = TPCustFormLib, generations = 2)
TPsLibFormSusp <- generateTPs("library_formula", TPLibrary = TPCustFormLib, suspL)
TPsLibFormScr <- generateTPs("library_formula", TPLibrary = TPCustFormLib, fGroups)

TPsBTSusp <- generateTPs("biotransformer", suspL)
TPsBTScr <- generateTPs("biotransformer", fGroups)
TPsBTSuspMore <- generateTPs("biotransformer", patRoonData::suspectsPos[1:25, ], calcSims = TRUE) # for filter tests

doCTS <- FALSE
if (doCTS)
{
    TPsCTSSusp <- generateTPs("cts", suspL, "photolysis_unranked")
    TPsCTSScr <- generateTPs("cts", fGroups, "photolysis_unranked")
}

fGroupsEmpty <- getEmptyTestFGroups()
fGroupsScrEmpty <- doScreen(fGroupsEmpty, data.table(name = "doesnotexist", SMILES = "C", mz = 12))
TPsLogicEmpty <- doGenLogicTPs(fGroupsEmpty)
TPsLibEmpty <- generateTPs("library", fGroupsScrEmpty)
TPsLibFormEmpty <- generateTPs("library_formula", TPLibrary = TPCustFormLib, fGroupsScrEmpty)
TPsBTEmpty <- generateTPs("biotransformer", fGroupsScrEmpty)
if (doCTS)
    TPsCTSEmpty <- generateTPs("cts", fGroupsScrEmpty, "photolysis_unranked")

doMetFrag <- TRUE # !is.null(getOption("patRoon.path.MetFragCL")) && nzchar(getOption("patRoon.path.MetFragCL"))

if (doMetFrag)
{
    plists <- generateMSPeakLists(fGroups, "mzr")
    plistsEmpty <- plists[FALSE, reAverage = TRUE]
    plistsEmptyMS <- removeMSPlists(plists, "MS")
    
    compsMF <- callMF(fGroups, plists)
    TPsLibComp <- generateTPs("library", compsMF)
    TPsLibCompEmpty <- generateTPs("library", delete(compsMF))
    TPsBTComp <- generateTPs("biotransformer", compsMF)
    TPsBTCompEmpty <- generateTPs("biotransformer", delete(compsMF))
    if (doCTS)
    {
        TPsCTSComp <- generateTPs("cts", compsMF, "photolysis_unranked")
        TPsCTSCompEmpty <- generateTPs("cts", delete(compsMF), "photolysis_unranked")
    }
}

test_that("verify TP generation", {
    expect_known_value(TPsLogic, testFile("tp-logic"))
    expect_known_value(TPsLibPC, testFile("tp-lib_pc"))
    expect_known_value(TPsLibSusp, testFile("tp-lib_susp"))
    expect_known_value(TPsLibScr, testFile("tp-lib_scr"))
    expect_known_value(TPsLibCustom, testFile("tp-lib_custom"))
    expect_known_value(TPsLibForm, testFile("tp-lib_form"))
    expect_known_value(TPsLibFormSusp, testFile("tp-lib_form_susp"))
    expect_known_value(TPsLibFormScr, testFile("tp-lib_form_scr"))
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
    expect_range(as.data.table(TPsLibPCGen2)$generation, 1:2)
    expect_equal(parents(TPsLibCustom)$name, "Diuron")
    expect_match(TPsLibCustom[[1]]$name_lib, "Linuron")
    expect_equal(TPsLibCustom, generateTPs("library", TPLibrary = TPCustLib, matchParentsBy = "InChI"))
    expect_equal(TPsLibPCGen2, generateTPs("library", generations = 2, matchGenerationsBy = "InChI"))

    checkmate::expect_names(parents(TPsBTSusp)$name, subset.of = suspL$name)
    checkmate::expect_names(parents(TPsBTScr)$name, subset.of = screenInfo(fGroups)$name)
    checkmate::expect_names(names(parents(TPsBTSusp)), must.include = c("name", "SMILES", "InChI", "InChIKey",
                                                                        "formula", "neutralMass"))
    
    expect_length(TPsLogicEmpty, 0)
    expect_length(TPsLibEmpty, 0)
    expect_length(TPsLibFormEmpty, 0)
    expect_length(TPsBTEmpty, 0)
    
    skip_if_not(doMetFrag)

    expect_known_value(TPsLibComp, testFile("tp-lib_comp"))
    expect_known_value(TPsBTComp, testFile("tp-bt_comp"))
    expect_known_show(TPsLibComp, testFile("tp-lib_comp", text = TRUE))
    expect_known_show(TPsBTComp, testFile("tp-bt_comp", text = TRUE))

    checkmate::expect_names(parents(TPsLibComp)$name, subset.of = as.data.table(compsMF)$identifier)
    checkmate::expect_names(parents(TPsBTComp)$name, subset.of = as.data.table(compsMF)$identifier)
    expect_length(TPsLibCompEmpty, 0)
    expect_length(TPsBTCompEmpty, 0)
    
    skip_if_not(doCTS)
    expect_known_value(TPsCTSSusp, testFile("tp-cts_susp"))
    expect_known_value(TPsCTSScr, testFile("tp-cts_scr"))
    expect_known_show(TPsCTSSusp, testFile("tp-cts_susp", text = TRUE))
    expect_known_show(TPsCTSScr, testFile("tp-cts_scr", text = TRUE))
    expect_length(TPsCTSEmpty, 0)
    expect_known_value(TPsCTSComp, testFile("tp-cts_comp"))
    expect_known_show(TPsCTSComp, testFile("tp-cts_comp", text = TRUE))
    checkmate::expect_names(parents(TPsCTSComp)$name, subset.of = as.data.table(compsMF)$identifier)
    expect_length(TPsCTSCompEmpty, 0)
})

assertSusp <- function(...) patRoon:::assertSuspectList(..., needsAdduct = FALSE, skipInvalid = FALSE)
testMFDB <- function(TPs)
{
    getMFDB <- function(...) { outf <- tempfile(fileext = ".csv"); convertToMFDB(..., out = outf); fread(outf) }
    
    MFDBNP <- getMFDB(TPs, includeParents = FALSE)
    MFDBP <- getMFDB(TPs, includeParents = TRUE)
    
    cols <- c("Identifier", "MolecularFormula", "MonoisotopicMass",
              "SMILES", "InChI", "InChIKey", "InChIKey1", "molNeutralized", "ALogP", "LogP", "XLogP", "parent",
              "transformation", "enzyme", "evidencedoi")
    
    unTPN <- uniqueN(as.data.table(TPs), by = "InChIKey")
    
    checkmate::expect_data_table(MFDBNP, any.missing = FALSE, nrows = unTPN)
    checkmate::expect_names(names(MFDBNP), subset.of = cols)
    
    checkmate::expect_data_table(MFDBP, nrows = unTPN + nrow(parents(TPs)))
    checkmate::expect_names(names(MFDBP), subset.of = cols)
}

TPsBTSuspFPI <- filter(TPsBTSuspMore, removeParentIsomers = TRUE)
TPsBTSuspFPIN <- filter(TPsBTSuspMore, removeParentIsomers = TRUE, negate = TRUE)
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
    
    expect_length(delete(TPsLogic, i = names(TPsLogic)), 0)
    expect_false(delete(TPsLogic, j = 1)[[1]]$name[1] == TPsLogic[[1]]$name[1])
    expect_length(delete(TPsLogic, j = 3:4), length(TPsLogic) - (length(names(TPsLogic)) * 2))
    expect_false(delete(TPsLogic, j = function(...) 1)[[1]]$name[1] == TPsLogic[[1]]$name[1])
    expect_length(delete(TPsLogic, j = function(...) 3:4), length(TPsLogic) - (length(names(TPsLogic)) * 2))
    expect_length(delete(TPsLogic, j = function(...) TRUE), 0)
    expect_equal(delete(TPsLogic, i = character()), TPsLogic)
    expect_equal(delete(TPsLogic, j = integer()), TPsLogic)
    expect_length(delete(TPsLogic), 0)
    expect_length(delete(TPsLibEmpty), 0)
    
    expect_setequal(as.data.table(filter(TPsLibPC, properties = list(retDir = c(-1, 1))))$retDir, c(-1, 1))
    expect_setequal(as.data.table(filter(TPsLibPC, properties = list(retDir = c(-1, 1)), negate = TRUE))$retDir, 0)
    expect_equal(anyDuplicated(as.data.table(filter(TPsBTSuspMore, removeDuplicates = TRUE)), by = c("SMILES", "parent")), 0)
    expect_gt(anyDuplicated(as.data.table(filter(TPsBTSuspMore, removeDuplicates = TRUE, negate = TRUE)),
                              by = c("SMILES", "parent")), 0)
    expect_false(any(mapply(parents(TPsBTSuspFPI)$formula, products(TPsBTSuspFPI), FUN = function(f, p) any(f %in% p$formula))))
    expect_true(all(mapply(parents(TPsBTSuspFPIN)$formula, products(TPsBTSuspFPIN), FUN = function(f, p) all(f %in% p$formula))))
    expect_equal(anyDuplicated(as.data.table(filter(TPsBTSuspMore, removeTPIsomers = TRUE)), by = c("formula", "parent")), 0)
    expect_gt(anyDuplicated(as.data.table(filter(TPsBTSuspMore, removeTPIsomers = TRUE, negate = TRUE)),
                              by = c("formula", "parent")), 0)
    expect_gte(min(as.data.table(filter(TPsBTSuspMore, minSimilarity = 0.5))$similarity), 0.5)
    expect_lt(max(as.data.table(filter(TPsBTSuspMore, minSimilarity = 0.5, negate = TRUE))$similarity), 0.5)
    
    skip_if_not(doMetFrag)
    
    expect_error(assertSusp(convertToSuspects(TPsLibComp, includeParents = TRUE)), NA)
    expect_error(assertSusp(convertToSuspects(TPsBTComp, includeParents = TRUE)), NA)
    testMFDB(TPsLibComp)
    testMFDB(TPsBTComp)
    
    skip_if_not(doCTS)
    expect_error(assertSusp(convertToSuspects(TPsCTSSusp, includeParents = TRUE)), NA)
    testMFDB(TPsCTSSusp)
    expect_error(assertSusp(convertToSuspects(TPsCTSComp, includeParents = TRUE)), NA)
    testMFDB(TPsCTSComp)
})

TPFormLib <- genFormulaTPLibrary(patRoonData::suspectsPos[1:4, ])
TPFormLibG2 <- genFormulaTPLibrary(patRoonData::suspectsPos[1:4, ], generations = 2)
genTPsFormLibG2 <- generateTPs("library_formula", parents = patRoonData::suspectsPos[1:4, ], TPLibrary = TPFormLibG2,
                               generations = 2)
test_that("genFormulaTPLibrary works", {
    checkmate::expect_data_table(TPFormLib, any.missing = FALSE)
    checkmate::expect_names(names(TPFormLib), must.include = c("parent_name", "parent_formula", "parent_neutralMass",
                                                               "TP_name", "TP_formula", "TP_neutralMass"))
    expect_gt(nrow(TPFormLibG2), nrow(TPFormLib))
    expect_setequal(TPFormLibG2$generation, c(1, 2))
    checkmate::expect_subset(TPFormLibG2[generation == 2]$parent_name, TPFormLibG2[generation == 1]$TP_name)
    
    expect_length(generateTPs("library_formula", TPLibrary = TPFormLib), nrow(TPFormLib))
    expect_setequal(as.data.table(genTPsFormLibG2)$generation, c(1, 2))
})

TPsCons <- consensus(TPsLibScr, TPsBTScr)
collapsedTPLen <- function(TPs) sum(sapply(products(TPs), function(x) uniqueN(x$InChIKey)))
test_that("consensus works", {
    expect_length(consensus(TPsLibScr, TPsBTEmpty), collapsedTPLen(TPsLibScr))
    
    expect_known_value(TPsCons, testFile("TPs-cons"))
    expect_known_show(TPsCons, testFile("TPs-cons", text = TRUE))
    expect_setequal(names(consensus(TPsLibScr, TPsBTScr)), union(names(TPsLibScr), names(TPsBTScr)))
    expect_lt(length(consensus(TPsLibScr, TPsBTScr, relMinAbundance = 1)), collapsedTPLen(TPsCons))
    expect_length(consensus(TPsLibEmpty, TPsBTEmpty), 0)
    
    expect_equal(sum(lengths(list(consensus(TPsLibScr, TPsBTScr, uniqueFrom = 1),
                                  consensus(TPsLibScr, TPsBTScr, uniqueFrom = 2),
                                  consensus(TPsLibScr, TPsBTScr, relMinAbundance = 1)))),
                 collapsedTPLen(TPsCons))
    expect_equal(sum(lengths(list(consensus(TPsLibScr, TPsBTScr, uniqueFrom = 1:2, uniqueOuter = TRUE),
                                  consensus(TPsLibScr, TPsBTScr, relMinAbundance = 1)))),
                 collapsedTPLen(TPsCons))
    expect_length(consensus(TPsLibScr, TPsBTScr, uniqueFrom = 1:2), collapsedTPLen(TPsCons))
    expect_lt(length(consensus(TPsLibScr, TPsBTScr, uniqueFrom = 1:2, uniqueOuter = TRUE)), collapsedTPLen(TPsCons))
    expect_length(consensus(TPsLibEmpty, TPsBTEmpty, uniqueFrom = 1), 0)
    expect_length(consensus(TPsLibEmpty, TPsBTEmpty, uniqueFrom = 1, uniqueOuter = TRUE), 0)
})

test_that("plotting works", {
    expect_HTML(plotGraph(TPsLibScr, which = 1))
    expect_HTML(plotGraph(TPsLibFormScr, which = 1))
    
    expect_doppel("venn", function() plotVenn(TPsLibScr, TPsBTScr))
    expect_error(plotVenn(TPsLibEmpty, TPsBTEmpty))
    expect_equal(expect_plot(plotVenn(TPsLibScr, TPsBTScr))$areas[2], collapsedTPLen(TPsBTScr))
    expect_equal(expect_plot(plotVenn(TPsLibScr, TPsBTEmpty))$areas[1], collapsedTPLen(TPsLibScr))
    expect_equal(expect_plot(plotVenn(TPsLibEmpty, TPsBTScr))$areas[2], collapsedTPLen(TPsBTScr))
    expect_equal(expect_plot(plotVenn(TPsLibScr, TPsBTEmpty))$intersectionCounts, 0)
    
    expect_ggplot(plotUpSet(TPsLibScr, TPsBTScr))
    expect_error(plotUpSet(TPsLibEmpty, TPsBTEmpty))
    expect_error(plotUpSet(TPsLibScr, TPsBTEmpty))
    
    expect_equal(expect_plot(plotVenn(TPsLibScr, TPsBTScr))$intersectionCounts,
                 length(consensus(TPsLibScr, TPsBTScr, relMinAbundance = 1)))
    
})


fGroupsMore <- getTestFGroups(getTestAnaInfoAnn())
componTPsNone <- generateComponents(fGroupsMore[, 1:50], "tp", TPs = NULL)
componTPsNoneTPDiff <- generateComponents(fGroupsMore[, 1:25], "tp", fGroupsMore[, 26:50], TPs = NULL)

componTPsLib <- generateComponents(doScreen(fGroupsMore, convertToSuspects(TPsLibSusp, includeParents = TRUE)), "tp",
                                   TPs = TPsLibSusp)
componTPsLibForm <- generateComponents(doScreen(fGroupsMore, convertToSuspects(TPsLibFormSusp, includeParents = TRUE)),
                                       "tp", TPs = TPsLibFormSusp)

TPsLogicMore <- doGenLogicTPs(fGroupsMore[, 1:50])
fGroupsTPsLogic <- doScreen(fGroupsMore, convertToSuspects(TPsLogicMore), onlyHits = TRUE)
componTPsLogic <- generateComponents(fGroupsTPsLogic, "tp", TPs = TPsLogicMore)

fGroupsMoreScr <- doScreen(fGroupsMore, convertToSuspects(TPsBTSusp, includeParents = TRUE), onlyHits = TRUE)
componTPsBT <- generateComponents(fGroupsMoreScr, "tp", TPs = TPsBTSusp)

if (doCTS)
    componTPsCTS <- generateComponents(fGroupsMoreScr, "tp", TPs = TPsCTSSusp)

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
    expect_known_value(componTPsLibForm, testFile("tp-compon-lib_form"))
    expect_known_value(componTPsBT, testFile("tp-compon-bt"))
    
    expect_known_show(componTPsNone, testFile("tp-compon-none", text = TRUE))
    expect_known_show(componTPsLogic, testFile("tp-compon-logic", text = TRUE))
    expect_known_show(componTPsLib, testFile("tp-compon-lib", text = TRUE))
    expect_known_show(componTPsLibForm, testFile("tp-compon-lib_form", text = TRUE))
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
    
    skip_if_not(doCTS)
    expect_known_value(componTPsCTS, testFile("tp-compon-cts"))
    expect_known_show(componTPsCTS, testFile("tp-compon-cts", text = TRUE))
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
    expect_HTML(plotGraph(TPsLibSusp, which = componentInfo(componTPsLib)$parent_name[1], components = componTPsLib))
    expect_HTML(plotGraph(TPsLibFormSusp, which = componentInfo(componTPsLibForm)$parent_name[1],
                          components = componTPsLibForm))
    
    expect_equal(as.data.table(componTPsLogic)[TP_retDir == retDir | TP_retDir == 0 | retDir == 0, -"size"],
                 as.data.table(componTPsRetF)[, -"size"])
    expect_equal(as.data.table(componTPsLogic)[TP_retDir != retDir | TP_retDir == 0 | retDir == 0, -"size"],
                 as.data.table(componTPsRetFN)[, -"size"])
    
    expect_reportHTML(makeReportHTML(fGroupsMoreScr, components = componTPsBT, TPs = TPsBTSusp))
    
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
