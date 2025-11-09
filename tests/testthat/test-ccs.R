# SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

local_edition(3) # for snapshots

test_that("CCS conversion works", {
    susps <- as.data.table(patRoonDataIMS::suspectsPos); setnames(susps, "mobility_[M+H]+", "mobility_susp")
    suspsExp <- expandSuspMobilities(susps)
    mzs <- suspsExp$mz
    mobs <- suspsExp$mobility
    CCSs <- convertMobilityToCCS(mobs, mzs, getCCSParams("mason-schamp_1/k"))
    
    # copied from test_data/agilent_calib/AcqData/OverrideImsCal.xml
    agilentCalib <- list(massGas = 28.006148, TFix = -0.10004857353171559, beta = 0.1300525995297081)
    expect_equal(convertMobilityToCCS(mobs, mzs, getCCSParams("agilent", calibrant = file.path(getTestDataPath(), "agilent_calib"))),
                 convertMobilityToCCS(mobs, mzs, getCCSParams("agilent", calibrant = file.path(getTestDataPath(), "agilent_calib", "AcqData", "OverrideImsCal.xml"))))
    expect_equal(convertMobilityToCCS(mobs, mzs, getCCSParams("agilent", calibrant = file.path(getTestDataPath(), "agilent_calib"))),
                 convertMobilityToCCS(mobs, mzs, getCCSParams("agilent", calibrant = agilentCalib)))
    expect_equal(convertCCSToMobility(CCSs, mzs, getCCSParams("agilent", calibrant = file.path(getTestDataPath(), "agilent_calib"))),
                 convertCCSToMobility(CCSs, mzs, getCCSParams("agilent", calibrant = file.path(getTestDataPath(), "agilent_calib", "AcqData", "OverrideImsCal.xml"))))
    expect_equal(convertCCSToMobility(CCSs, mzs, getCCSParams("agilent", calibrant = file.path(getTestDataPath(), "agilent_calib"))),
                 convertCCSToMobility(CCSs, mzs, getCCSParams("agilent", calibrant = agilentCalib)))
    
    testRoundtrip <- function(CCSParams)
    {
        m2C <- convertMobilityToCCS(mobs, mzs, CCSParams)
        c2M <- convertCCSToMobility(CCSs, mzs, CCSParams)
        expect_equal(convertMobilityToCCS(c2M, mzs, CCSParams), CCSs)
        expect_equal(convertCCSToMobility(m2C, mzs, CCSParams), mobs)
    }
    testRoundtrip(getCCSParams("agilent", calibrant = agilentCalib))
    testRoundtrip(getCCSParams("mason-schamp_k"))
    testRoundtrip(getCCSParams("mason-schamp_1/k"))
    
    testDifferent <- function(CCSParams1, CCSParams2)
    {
        expect_false(any(convertMobilityToCCS(mobs, mzs, CCSParams1) == convertMobilityToCCS(mobs, mzs, CCSParams2)))
        expect_false(any(convertCCSToMobility(CCSs, mzs, CCSParams1) == convertCCSToMobility(CCSs, mzs, CCSParams2)))
    }
    testDifferent(getCCSParams("agilent", calibrant = agilentCalib),
                  getCCSParams("agilent", calibrant = modifyList(agilentCalib, list(massGas = 32))))
    testDifferent(getCCSParams("agilent", calibrant = agilentCalib, defaultCharge = 1),
                  getCCSParams("agilent", calibrant = agilentCalib, defaultCharge = 2))
    
    testDifferent(getCCSParams("mason-schamp_k", massGas = 28), getCCSParams("mason-schamp_k", massGas = 32))
    testDifferent(getCCSParams("mason-schamp_k", defaultCharge = 1), getCCSParams("mason-schamp_k", defaultCharge = 2))
    testDifferent(getCCSParams("mason-schamp_k", MasonSchampConst = 18500), getCCSParams("mason-schamp_k", MasonSchampConst = 18600))
    testDifferent(getCCSParams("mason-schamp_1/k", massGas = 28), getCCSParams("mason-schamp_k", massGas = 32))
    testDifferent(getCCSParams("mason-schamp_1/k", defaultCharge = 1), getCCSParams("mason-schamp_k", defaultCharge = 2))
    testDifferent(getCCSParams("mason-schamp_1/k", MasonSchampConst = 18500), getCCSParams("mason-schamp_k", MasonSchampConst = 18600))
    
    expect_equal(convertMobilityToCCS(mobs, mzs, getCCSParams("mason-schamp_k", defaultCharge = 2)),
                 convertMobilityToCCS(mobs, mzs, getCCSParams("mason-schamp_k", defaultCharge = 1), charge = rep(2, length(mobs))))
    expect_equal(convertCCSToMobility(CCSs, mzs, getCCSParams("mason-schamp_k", defaultCharge = 2)),
                 convertCCSToMobility(CCSs, mzs, getCCSParams("mason-schamp_k", defaultCharge = 1), charge = rep(2, length(CCSs))))
    
    skip_on_os("mac")
    testRoundtrip(getCCSParams("bruker"))
    testDifferent(getCCSParams("bruker", defaultCharge = 1), getCCSParams("bruker", defaultCharge = 2))
    expect_equal(convertMobilityToCCS(mobs, mzs, getCCSParams("bruker")),
                 convertMobilityToCCS(mobs, mzs, getCCSParams("mason-schamp_1/k")), tolerance = 1)
    expect_equal(convertCCSToMobility(CCSs, mzs, getCCSParams("bruker")),
                 convertCCSToMobility(CCSs, mzs, getCCSParams("mason-schamp_1/k")), tolerance = 1e-3)
})

test_that("assignMobilities() for suspects", {
    susps <- as.data.table(patRoonDataIMS::suspectsPos)
    susps <- assignMobilities(susps, CCSParam = getCCSParams("mason-schamp_1/k"))
    # NOTE: row 1 (hydroxy carbamazepine) is not found in PCL
    suspsNoMob <- susps[2:6, -c("mobility_[M+H]+", "CCS_[M+H]+")]
    
    suspsPredC3SDB <- assignMobilities(suspsNoMob, "c3sdb")
    suspsPredPCL <- assignMobilities(suspsNoMob, "pubchemlite")
    suspsPredSelf <- assignMobilities(suspsNoMob, susps)
    
    expect_snapshot_value(list(
        suspsPredC3SDB,
        suspsPredPCL,
        suspsPredSelf
    ), style = "serialize")
    
    verifyPred <- function(suspsPred)
    {
        checkmate::expect_data_table(suspsPred, nrows = nrow(suspsNoMob))
        checkmate::expect_names(names(suspsPred), must.include = "CCS_[M+H]+")
        if (is.character(suspsPred[["CCS_[M+H]+"]]))
            checkmate::expect_character(suspsPred[["CCS_[M+H]+"]], any.missing = FALSE)
        else
            checkmate::expect_numeric(suspsPred[["CCS_[M+H]+"]], any.missing = FALSE, lower = 0)
    }

    verifyPred(suspsPredC3SDB)
    verifyPred(suspsPredPCL)
    verifyPred(suspsPredSelf)
    
    expect_equal(assignMobilities(suspsNoMob, NULL), suspsNoMob)
    expect_equal(assignMobilities(suspsNoMob, susps, matchFromBy = "SMILES"), suspsPredSelf)
    expect_equal(assignMobilities(suspsPredC3SDB, susps[, -"mobility_[M+H]+"], overwrite = FALSE), suspsPredC3SDB)
    expect_false(isTRUE(all.equal(assignMobilities(suspsPredC3SDB, susps, overwrite = TRUE), suspsPredC3SDB)))
    expect_equal(assignMobilities(suspsPredC3SDB, susps, overwrite = TRUE)[["CCS_[M+H]+"]], suspsPredSelf[["CCS_[M+H]+"]])
    # NOTE: susps CCSs are characters, so type will change
    expect_equal(assignMobilities(suspsPredC3SDB, susps[1:2], overwrite = TRUE)[["CCS_[M+H]+"]][-(1:2)],
                 as.character(suspsPredC3SDB[["CCS_[M+H]+"]][-(1:2)]))
    # check if NA values in from are not overwritten
    suspsPredC3SDBNA <- copy(suspsPredC3SDB); suspsPredC3SDBNA[2, "CCS_[M+H]+" := NA]
    checkmate::expect_numeric(assignMobilities(suspsPredC3SDB, suspsPredC3SDBNA, overwrite = TRUE)[["CCS_[M+H]+"]],
                              lower = 0, any.missing = FALSE)
    
    checkmate::expect_names(names(assignMobilities(suspsNoMob, "c3sdb", adducts = c("[M+H]+", "[M+K]+"))),
                            must.include = c("CCS_[M+H]+", "CCS_[M+K]+"))
    suspsNoMobAdd <- copy(suspsNoMob)
    suspsNoMobAdd[, adduct := "[M+H]+"]
    suspsNoMobAdd[2, adduct := "[M+K]+"]
    suspsPredAddAO <- assignMobilities(suspsNoMobAdd, "c3sdb", predictAdductOnly = TRUE)
    checkmate::expect_names(names(suspsPredAddAO), must.include = c("CCS_[M+H]+", "CCS_[M+K]+"))
    checkmate::expect_scalar_na(suspsPredAddAO[["CCS_[M+H]+"]][2])
    checkmate::expect_number(suspsPredAddAO[["CCS_[M+K]+"]][2], lower = 0, na.ok = FALSE)
    expect_true(all(!is.na(suspsPredAddAO[["CCS_[M+H]+"]][-2])))
    expect_true(all(is.na(suspsPredAddAO[["CCS_[M+K]+"]][-2])))
    suspsPredAddNAO <- assignMobilities(suspsNoMobAdd, "c3sdb", adducts = c("[M+H]+", "[M+K]+"), predictAdductOnly = FALSE)
    checkmate::expect_numeric(suspsPredAddNAO[["CCS_[M+H]+"]], lower = 0, any.missing = FALSE)
    checkmate::expect_numeric(suspsPredAddNAO[["CCS_[M+K]+"]], lower = 0, any.missing = FALSE)
    suspsPredAddAOOne <- assignMobilities(suspsNoMobAdd, "c3sdb", predictAdductOnly = TRUE, adducts = character())
    expect_equal(suspsPredAddAOOne, suspsPredAddAO)
    
    suspsUnSpec <- copy(susps); setnames(suspsUnSpec, c("mobility_[M+H]+", "CCS_[M+H]+"), c("mobility", "CCS"))
    checkmate::expect_names(names(assignMobilities(suspsNoMob, suspsUnSpec, adducts = NA)),
                            must.include = c("mobility", "CCS"))

    CCSParams <- getCCSParams("mason-schamp_1/k")
    expect_equal(assignMobilities(suspsNoMob, CCSParams = NULL), suspsNoMob)
    expect_equal(assignMobilities(suspsNoMob, CCSParams = CCSParams), suspsNoMob)
    checkmate::expect_character(assignMobilities(susps[, -"mobility_[M+H]+"], CCSParams = CCSParams)[["mobility_[M+H]+"]], any.missing = FALSE)
    checkmate::expect_character(assignMobilities(susps[, -"CCS_[M+H]+"], CCSParams = CCSParams)[["CCS_[M+H]+"]], any.missing = FALSE)
    suspsMobNA <- copy(susps); suspsMobNA[2, "mobility_[M+H]+" := NA]
    checkmate::expect_character(assignMobilities(suspsMobNA, CCSParams = CCSParams)[["mobility_[M+H]+"]], any.missing = FALSE)
    suspsMobAllNA <- copy(susps); suspsMobAllNA[, "mobility_[M+H]+" := NA]
    checkmate::expect_character(assignMobilities(suspsMobAllNA, CCSParams = CCSParams)[["mobility_[M+H]+"]], any.missing = FALSE)
    # first suspect has two CCSs
    suspsMob1st <- assignMobilities(susps[1, -"mobility_[M+H]+"], CCSParams = CCSParams)
    checkmate::expect_data_table(expandSuspMobilities(setnames(suspsMob1st, "CCS_[M+H]+", "CCS_susp")), nrows = 2)
    suspsCh <- copy(susps); setnames(suspsCh, "CCS_[M+H]+", "CCS_[M+H]2+")
    suspsNoMobCh <- assignMobilities(suspsNoMob, suspsCh, adducts = "[M+H]2+", CCSParams = CCSParams)
    suspsNoMobCh2 <- copy(suspsNoMob); suspsNoMobCh2[, mz := calculateMasses(neutralMass, as.adduct("[M+H]2+"), "mz")]
    suspsChUnSpec <- copy(suspsCh); setnames(suspsChUnSpec, "CCS_[M+H]2+", "CCS")
    suspsNoMobCh2 <- assignMobilities(suspsNoMobCh2, suspsChUnSpec, adducts = NA,
                                      CCSParams = modifyList(CCSParams, list(defaultCharge = 2)))
    expect_equal(suspsNoMobCh[["CCS_[M+H]2+"]], suspsNoMobCh2[["CCS"]])
    
    # disable caching to get warnings (UNDONE?)
    withOpt(cache.mode = "none", {
        expect_warning(assignMobilities(suspsNoMob, susps, adducts = "[M-H]-"), "relevant") # susps only has M+H
        
        suspsNoMobMZ <- suspsNoMob[, -"mz"]
        expect_warning(assignMobilities(suspsNoMobMZ[, -"neutralMass"], susps, CCSParams = CCSParams,
                                        prepareChemProps = FALSE),
                       "no neutralMass data", fixed = TRUE) |>
            expect_warning("no mz data", fixed = TRUE) |>
            expect_warning("relevant")
        suspsNoMobNMNA <- copy(suspsNoMobMZ); suspsNoMobNMNA[2, neutralMass := NA]
        expect_warning(assignMobilities(suspsNoMobNMNA, susps, CCSParams = CCSParams, prepareChemProps = FALSE),
                       "rows with NA neutralMass values: 2", fixed = TRUE) |>
            expect_warning("no mz data", fixed = TRUE)
        expect_warning(assignMobilities(suspsNoMobMZ, suspsUnSpec, CCSParams = CCSParams, prepareChemProps = FALSE),
                       "no mz data")
        suspsNoMobMZNA <- copy(suspsNoMob); suspsNoMobMZNA[2, mz := NA]
        expect_warning(assignMobilities(suspsNoMobMZNA, suspsUnSpec, CCSParams = CCSParams, prepareChemProps = FALSE),
                       "rows with NA mz values: 2", fixed = TRUE)
        expect_warning(assignMobilities(suspsNoMobMZ[, -"SMILES"], "c3sdb", prepareChemProps = FALSE), "skipped")
        suspsNoMoSMNA <- copy(suspsNoMobMZ); suspsNoMoSMNA[2, SMILES := NA]
        expect_warning(assignMobilities(suspsNoMoSMNA, "c3sdb", prepareChemProps = FALSE), "C3SDB predictions: 2", fixed = TRUE)
    })    
})
