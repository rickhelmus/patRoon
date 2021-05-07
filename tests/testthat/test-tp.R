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

plists <- generateMSPeakLists(fGroups, "mzr")
plistsEmpty <- plists[FALSE, reAverage = TRUE]
plistsEmptyMS <- removeMSPlists(plists, "MS")
fGroupsEmpty <- getEmptyTestFGroups()
fGroupsScrEmpty <- doScreen(fGroupsEmpty, data.table(name = "doesnotexist", SMILES = "C", mz = 12))

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
    
    expect_length(generateTPs("logic", fGroupsEmpty), 0)
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
