context("TPs")

fGroups <- getCompFGroups()

suspL <- patRoonData::suspectsPos[patRoonData::suspectsPos$name %in% screenInfo(fGroups)$name, ]

TPsLogic <- generateTPs("logic", fGroups)
TPsLogicCustom <- generateTPs("logic", fGroups, transformations = data.table(transformation = "test", add = "C",
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
TPsLogicEmpty <- generateTPs("logic", fGroupsEmpty)

plists <- generateMSPeakLists(fGroups, "mzr")
plistsEmpty <- plists[FALSE, reAverage = TRUE]
plistsEmptyMS <- removeMSPlists(plists, "MS")

doMetFrag <- !is.null(getOption("patRoon.path.MetFragCL")) && nzchar(getOption("patRoon.path.MetFragCL"))

if (doMetFrag)
{
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
    expect_true(all(as.data.table(generateTPs("logic", fGroups, minMass = 100))$neutralMass >= 100))
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
    expect_length(generateTPs("library", fGroupsScrEmpty), 0)
    expect_length(generateTPs("biotransformer", fGroupsScrEmpty), 0)
    
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
              "Precursor MonoisotopicMass", "SMILES", "InChI", "InChIKey",
              "InChIKey1", "ALogP")
    
    checkmate::expect_data_table(MFDBNP, any.missing = FALSE, nrows = length(TPsLibScr))
    checkmate::expect_names(names(MFDBNP), subset.of = cols)
    
    convertToMFDB(..., out = outf, includeParents = TRUE); db <- fread(outf)
    checkmate::expect_data_table(MFDBP, any.missing = FALSE, nrows = length(TPsLibScr) + nrow(parents(TPsLibScr)))
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
    expect_error(assertSusp(convertToSuspects(TPsLibScr)), NA)
    expect_error(assertSusp(convertToSuspects(TPsBTSusp)), NA)
    
    testMFDB(TPsLibScr)
    testMFDB(TPsBTSusp)
    
    expect_false(any(mapply(parents(TPsBTSuspFEF)$formula, products(TPsBTSuspFEF), FUN = function(f, p) any(f %in% p$formula))))
    expect_true(all(mapply(parents(TPsBTSuspFEFN)$formula, products(TPsBTSuspFEFN), FUN = function(f, p) all(f %in% p$formula))))
    expect_gte(min(as.data.table(filter(TPsBTSuspMore, minSimilarity = 0.5))$similarity), 0.5)
    expect_lt(max(as.data.table(filter(TPsBTSuspMore, minSimilarity = 0.5, negate = TRUE))$similarity), 0.5)
    
    skip_if_not(doMetFrag)
    
    expect_error(assertSusp(convertToSuspects(TPsLibComp)), NA)
    expect_error(assertSusp(convertToSuspects(TPsBTComp)), NA)
})
