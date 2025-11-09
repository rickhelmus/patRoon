# SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

context("CCS")

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
    suspsNoMob <- susps[1:5, -"mobility_[M+H]+"]
    
    suspsPredC3SDB <- assignMobilities(suspsNoMob, "c3sdb")
    suspsPredPCL <- assignMobilities(suspsNoMob, "pubchemlite")
    suspsPredSelf <- assignMobilities(suspsNoMob, susps)
    
    expect_snapshot_value(list(
        suspsPredC3SDB,
        suspsPredPCL,
        suspsPredSelf
    ))
    
    verifyPred <- function(suspsPred)
    {
        checkmate::expect_data_table(suspsPred, nrows = nrow(suspsNoMob), ncols = ncol(suspsNoMob) + 1)
        checkmate::expect_names(names(suspsPred), must.include = "CCS_[M+H]+")
        checkmate::expect_numeric(suspsPred[["CCS_[M+H]+"]], any.missing = FALSE, lower = 0)
    }
    verifyPred(suspsPredC3SDB)
    verifyPred(suspsPredPCL)
    verifyPred(suspsPredSelf)
    
    expect_equal(assignMobilities(suspsNoMob, NULL), suspsNoMob)
    expect_equal(assignMobilities(suspsNoMob, susps, matchFromBy = "SMILES"), suspsPredSelf)
    expect_equal(assignMobilities(suspsPredC3SDB, susps, overwrite = FALSE), suspsPredC3SDB)
    epxect_false(isTRUE(all.equal(assignMobilities(suspsPredC3SDB, susps, overwrite = TRUE), suspsPredC3SDB)))
    expect_equal(assignMobilities(suspsPredC3SDB, susps, overwrite = TRUE), suspsPredSelf)
    expect_equal(assignMobilities(suspsPredC3SDB, susps[1:2], overwrite = TRUE)[-(1:2)], suspsPredC3SDB[-(1:2)])
})
