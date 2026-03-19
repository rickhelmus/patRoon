# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

initXCMS()

anaInfo <- getTestAnaInfo()
anaInfoOne <- getTestAnaInfo()[4, ]

ffOpenMS <- getTestFeatures(anaInfo)
# ffXCMS <- findFeatures(anaInfoOne, "xcms", noise = 3E4)
ffXCMS3 <- findFeatures(anaInfoOne, "xcms3", xcms::CentWaveParam(noise = 3E4))
# UNDONE: ignore warnings about clusters...
ffKPIC2 <- withCallingHandlers(findFeatures(anaInfoOne, "kpic2", level = 1E5),
                               warning = function(w) if (grepl("number of clusters", w, fixed = TRUE)) invokeRestart("muffleWarning"))
ffSIRIUS <- findFeatures(anaInfoOne, "sirius")

# generate mzXML files for enviPick
exDataFiles <- list.files(patRoonData::exampleDataPath(), "\\.mzML$", full.names = TRUE)
epAnaInfo <- makeMZXMLs(anaInfoOne)
ffEP <- findFeatures(epAnaInfo, "envipick", minpeak = 25)

getPiek <- \(genp, anaInfo, peakp = getDefPeakParams("chrom", "piek"), ...) findFeatures(anaInfo, "piek", genp, peakp, ...)
getPiekHRMS <- \(genp, ...) getPiek(genp, anaInfoOne, ..., suspects = patRoonData::suspectsPos, adduct = "[M+H]+")
ffPiekBins <- getPiekHRMS(getPiekEICParams(mzRange = c(200, 300)))
ffPiekSusp <- getPiekHRMS(getPiekEICParams(filter = "suspects"))
ffPiekMS2 <- getPiekHRMS(getPiekEICParams(filter = "ms2"))
ffPiekMS2OpenMS <- getPiekHRMS(getPiekEICParams(filter = "ms2"), getDefPeakParams("chrom", "openms"))
ffPiekMS2XCMS3 <- getPiekHRMS(getPiekEICParams(filter = "ms2"), getDefPeakParams("chrom", "xcms3"))
ffPiekMS2EP <- getPiekHRMS(getPiekEICParams(filter = "ms2"), getDefPeakParams("chrom", "envipick"))

anaInfoOneIMS <- getTestAnaInfoIMS()[4, ]
getPiekIMS <- \(..., args = list()) do.call(getPiek, c(list(getPiekEICParams(...), anaInfoOneIMS, IMS = TRUE,
                                                            suspects = patRoonDataIMS::suspectsPos, adduct = "[M+H]+"),
                                                       args))
getPiekSuspIMS <- \(...) getPiekIMS(filter = "suspects", filterIMS = "suspects", ...)
ffPiekBinsIMS <- getPiekIMS(mzRange = c(200, 300), mobRange = c(0.6, 0.8))
ffPiekSuspIMS <- getPiekSuspIMS()
ffPiekMS2IMS <- getPiekIMS(filter = "ms2", filterIMS = "ms2")

ffOpenMSQ <- calculatePeakQualities(ffOpenMS)

ffEmpty <- getTestFeatures(anaInfoOne, noiseThrInt = 1E9)
ffEmptyQ <- calculatePeakQualities(ffEmpty)

if (doDATests())
{
    ffDA <- findFeatures(getDAAnaInfo("std1")[1, ], "bruker")

    # NOTE: use 2nd analysis here so first can be re-used for MS peaklists/formulas...
    ffDAEmpty <- findFeatures(getDAAnaInfo("std2")[1, ], "bruker", endRange = 0.01, doFMF = "force")
}

# Remove ID column: not reproducible
OpenMSFTable <- function(ff) sapply(featureTable(ff), function(fts) fts[, -"ID"], simplify = FALSE)

test_that("verify feature finder output", {
    expect_known_val(OpenMSFTable(ffOpenMS), "ff-openms", tolerance = 1E-5) # increased tolerance value for win/lin deviations
    # expect_known_val(featureTable(ffXCMS), "ff-xcms")
    expect_known_val(featureTable(ffXCMS3), "ff-xcms3")
    expect_known_val(featureTable(ffEP), "ff-envipick")
    expect_known_val(featureTable(ffKPIC2), "ff-kpic2")
    expect_known_val(featureTable(ffSIRIUS), "ff-sirius")
    
    expect_known_val(OpenMSFTable(ffOpenMSQ), "ff-openms-qual")
    
    # extraOpts
    expect_equal(OpenMSFTable(ffOpenMS),
                 OpenMSFTable(getTestFeatures(anaInfo,
                                              extraOpts = list("-algorithm:common:noise_threshold_int" = 30000))))

    skip_if_not(doDATests())
    expect_known_val(featureTable(ffDA), "ff-DA")
})

test_that("verify show output", {
    expect_known_show(ffOpenMS, "ff-show-openms")
    # expect_known_show(ffXCMS, "ff-show-xcms")
    expect_known_show(ffXCMS3, "ff-show-xcms3")
    expect_known_show(ffEP, "ff-show-envipick")
    expect_known_show(ffKPIC2, "ff-show-kpic2")
    expect_known_show(ffSIRIUS, "ff-show-sirius")

    expect_known_show(ffOpenMSQ, "ff-show-openms-qual")
    
    skip_if_not(doDATests())
    expect_known_show(ffDA, "ff-DA")
})

test_that("piek", {
    expect_known_val(featureTable(ffPiekBins), "ff-piek-bins")
    expect_known_val(featureTable(ffPiekSusp), "ff-piek-suspects")
    expect_known_val(featureTable(ffPiekMS2), "ff-piek-ms2")
    expect_known_val(featureTable(ffPiekMS2OpenMS), "ff-piek-ms2-openms")
    expect_known_val(featureTable(ffPiekMS2XCMS3), "ff-piek-ms2-xcms3")
    expect_known_val(featureTable(ffPiekMS2EP), "ff-piek-ms2-envipick")
    expect_known_show(ffPiekBinsIMS, "ff-show-piek-bins-ims")
    expect_known_show(ffPiekSuspIMS, "ff-show-piek-suspects-ims")
    expect_known_show(ffPiekMS2IMS, "ff-show-piek-ms2-ims")
    
    expect_known_show(ffPiekBins, "ff-show-piek-bins")
    expect_known_show(ffPiekSusp, "ff-show-piek-suspects")
    expect_known_show(ffPiekMS2, "ff-show-piek-ms2")
    expect_known_show(ffPiekMS2OpenMS, "ff-show-piek-ms2-openms")
    expect_known_show(ffPiekMS2XCMS3, "ff-show-piek-ms2-xcms3")
    expect_known_show(ffPiekMS2EP, "ff-show-piek-ms2-envipick")
    expect_known_val(featureTable(ffPiekBinsIMS), "ff-piek-bins-ims")
    expect_known_val(featureTable(ffPiekSuspIMS), "ff-piek-suspects-ims")
    expect_known_val(featureTable(ffPiekMS2IMS), "ff-piek-ms2-ims")
    
    expect_range(as.data.table(ffPiekBins)$mz, c(200, 300+0.02)) # +0.02: bin width
    expect_lte(length(ffPiekSusp), nrow(patRoonData::suspectsPos))
    expect_gt(length(ffPiekMS2), length(getPiekHRMS(getPiekEICParams(filter = "ms2", minTIC = 1E5))))
    expect_equal(ffPiekBins, withOpt(cache.mode="none",
                                     getPiekHRMS(getPiekEICParams(mzRange = c(200, 300)), EICBatchSize = 5E3)))
    expect_false(isTRUE(all.equal(ffPiekMS2, getPiekHRMS(getPiekEICParams(filter = "ms2"), assignRTWindow = 5),
                                  check.attributes = FALSE)))
    
    expect_true(hasIMS(ffPiekBinsIMS))
    expect_false(hasIMS(ffPiekBins))
    expect_range(as.data.table(ffPiekBinsIMS)$mz, c(200, 300+0.02))
    expect_range(as.data.table(ffPiekBinsIMS)$mobility, c(0.6, 0.8+0.08)) # +0.08: bin width
    fgPiekSusp <- screenSuspects(groupFeatures(ffPiekSuspIMS, "greedy"), patRoonDataIMS::suspectsPos)
    expect_setequal(screenInfo(fgPiekSusp)$name, patRoonDataIMS::suspectsPos$name)
    expect_setequal(screenInfo(fgPiekSusp)$group, names(fgPiekSusp))
    expect_gt(length(ffPiekMS2IMS), length(getPiekIMS(filter = "ms2", filterIMS = "ms2", minTIC = 1E5)))
    
    expect_range(as.data.table(getPiekHRMS(getPiekEICParams(filter = "ms2", retRange = c(60, 120))))$ret, c(60, 120))
    expect_lt(length(getPiekHRMS(getPiekEICParams(filter = "ms2", minEICIntensity = 1E5))), length(ffPiekMS2))
    expect_lt(length(getPiekHRMS(getPiekEICParams(filter = "ms2", topMostEICMZ = 25))), length(ffPiekMS2))
    expect_lte(max(getPiekHRMS(getPiekEICParams(filter = "ms2"), getDefPeakParams("chrom", "piek", forcePeakWidth = c(0, 10)))[[1]][, retmax - retmin]), 10)
    
    expect_false(isTRUE(all.equal(ffPiekSuspIMS, getPiekSuspIMS(sumWindowMZ = 10))))
    expect_false(isTRUE(all.equal(ffPiekSuspIMS, getPiekSuspIMS(sumWindowMob = 10))))
    expect_false(isTRUE(all.equal(ffPiekSuspIMS, getPiekSuspIMS(smoothWindowMZ = 1))))
    expect_false(isTRUE(all.equal(ffPiekSuspIMS, getPiekSuspIMS(smoothWindowMob = 10))))

    ffPiekDup <- getPiekIMS(filter = "ms2", args = list(keepDups = TRUE))
    checkmate::expect_names(names(ffPiekDup[[1]]), must.include = c("mzCentered", "dup"))
    expect_gt(nrow(ffPiekDup[[1]]), nrow(ffPiekSusp[[1]]))
    expect_lt(length(getPiekIMS(filter = "ms2", filterIMS = "ms2", args = list(rtWindowDup = 300))), length(ffPiekMS2IMS))
    expect_lt(length(getPiekIMS(filter = "ms2", filterIMS = "ms2", args = list(mzWindowDup = 0.1))), length(ffPiekMS2IMS))
    expect_lt(length(getPiekIMS(filter = "ms2", filterIMS = "ms2", args = list(mobWindowDup = 0.2))), length(ffPiekMS2IMS))
    # NOTE: the effect of minPeakOverlapDup is very small it seems...
    expect_lt(length(getPiekIMS(filter = "ms2", filterIMS = "ms2", args = list(rtWindowDup = 300, minPeakOverlapDup = 0.1))),
              length(getPiekIMS(filter = "ms2", filterIMS = "ms2", args = list(rtWindowDup = 300, minPeakOverlapDup = 0.9))))

    ffPiekIMSProfs <- getPiekSuspIMS(saveMZProfiles = TRUE, saveEIMs = TRUE)
    expect_length(ffPiekIMSProfs@mzProfiles[[1]], length(ffPiekIMSProfs[1]))
    expect_length(ffPiekIMSProfs@EIMs[[1]], length(ffPiekIMSProfs[1]))
})

test_that("verify empty object can be generated", {
    expect_length(ffEmpty, 0)
    #expect_length(suppressWarnings(findFeatures(anaInfoOne, "xcms", noise = 1E9)), 0)
    expect_length(findFeatures(epAnaInfo, "envipick", minint = 1E8, maxint = 1E9, minpeak = 100), 0) # add minpeak to speed-up
    expect_length(ffEmptyQ, 0)

    skip_if_not(doDATests())
    expect_length(ffDAEmpty, 0)
})

test_that("basic subsetting", {
    expect_length(ffOpenMS["nope"], 0)
    expect_equal(analyses(ffOpenMS[1:2]), anaInfo$analysis[1:2])
    expect_equal(analyses(ffOpenMS[anaInfo$analysis[2:3]]), anaInfo$analysis[2:3])
    expect_equal(analyses(ffOpenMS[c(TRUE, FALSE)]), anaInfo$analysis[c(TRUE, FALSE)])
    expect_equal(length(ffOpenMS[FALSE]), 0)
    expect_length(ffEmpty[1:5], 0)

    expect_equal(ffOpenMS[[2]], featureTable(ffOpenMS)[[2]])
    expect_equal(ffOpenMS[[analyses(ffOpenMS)[2]]], featureTable(ffOpenMS)[[2]])
    expect_equal(callDollar(ffOpenMS, analyses(ffOpenMS)[2]), ffOpenMS[[2]])
    
    expect_equal(replicates(ffOpenMS[, ni = replicate == "standard-pos"]), "standard-pos")
})

qr <- list(ZigZag = c(0.2, 0.6), TPASRScore = c(0.5, 0.9))
ffOpenMSQF <- filter(ffOpenMSQ, qualityRange = qr)
ffOpenMSQFTab <- as.data.table(ffOpenMSQF)
ffOpenMSQFN <- filter(ffOpenMSQ, qualityRange = qr, negate = TRUE)
ffOpenMSQFNTab <- as.data.table(ffOpenMSQFN)

test_that("delete and filter", {
    checkmate::expect_names(analyses(delete(ffOpenMS, i = 1)), disjunct.from = analyses(ffOpenMS)[1])
    checkmate::expect_names(analyses(delete(ffOpenMS, i = analyses(ffOpenMS)[1])), disjunct.from = analyses(ffOpenMS)[1])
    expect_length(delete(ffOpenMS, i = analyses(ffOpenMS)), 0)
    expect_false(delete(ffOpenMS, j = 1)[[1]]$ID[1] == ffOpenMS[[1]]$ID[1])
    expect_length(delete(ffOpenMS, j = 3:4), length(ffOpenMS) - (length(analyses(ffOpenMS)) * 2))
    expect_false(delete(ffOpenMS, j = function(...) 1)[[1]]$ID[1] == ffOpenMS[[1]]$ID[1])
    expect_length(delete(ffOpenMS, j = function(...) 3:4), length(ffOpenMS) - (length(analyses(ffOpenMS)) * 2))
    expect_length(delete(ffOpenMS, j = function(...) TRUE), 0)
    expect_equal(delete(ffOpenMS, i = character()), ffOpenMS)
    expect_equal(delete(ffOpenMS, j = integer()), ffOpenMS)
    expect_length(delete(ffOpenMS), 0)
    expect_length(delete(ffEmpty), 0)
    
    expect_gte(min(filter(ffOpenMS, absMinIntensity = 500)[[1]]$intensity), 500)
    expect_gte(min(filter(ffOpenMS, relMinIntensity = 0.2)[[1]]$intensity), 0.2 * max(ffOpenMS[[1]]$intensity))

    expect_range(filter(ffOpenMS, retentionRange = c(120, 300))[[1]]$ret, c(120, 300))
    expect_equal(filter(ffOpenMS, retentionRange = c(0, Inf)), ffOpenMS)
    expect_range(filter(ffOpenMS, mzRange = c(200, 300))[[1]]$mz, c(200, 300))
    expect_equal(filter(ffOpenMS, mzRange = c(0, Inf)), ffOpenMS)
    expect_range(filter(ffOpenMS, mzDefectRange = c(0.1, 0.2))[[1]]$mz %% 1, c(0.1, 0.2))
    expect_equal(filter(ffOpenMS, mzDefectRange = c(0, 1)), ffOpenMS)
    expect_lt(length(filter(ffOpenMS, chromWidthRange = c(0, 30))), length(ffOpenMS))
    expect_equal(filter(ffOpenMS, chromWidthRange = c(0, Inf)), ffOpenMS)

    expect_range(ffOpenMSQFTab[[names(qr)[1]]], qr[[1]])
    expect_range(ffOpenMSQFTab[[names(qr)[2]]], qr[[2]])
    expect_true(all(ffOpenMSQFNTab[[names(qr)[1]]] < qr[[c(1, 1)]] | ffOpenMSQFNTab[[names(qr)[1]]] > qr[[c(1, 2)]]))
    expect_true(all(ffOpenMSQFNTab[[names(qr)[2]]] < qr[[c(2, 1)]] | ffOpenMSQFNTab[[names(qr)[2]]] > qr[[c(2, 2)]]))
    
    expect_known_show(filter(ffOpenMS, absMinIntensity = 500, retentionRange = c(120, Inf),
                             mzRange = c(100, 400)),
                      "ff-combi")
    expect_known_show(filter(ffOpenMS, absMinIntensity = 500, retentionRange = c(120, Inf),
                             mzRange = c(100, 400), negate = TRUE),
                      "ff-combi-neg")
    expect_length(filter(ffEmpty, absMinIntensity = 500, retentionRange = c(120, Inf),
                         mzRange = c(100, 400)), 0)
    
    expect_range(as.data.table(filter(ffPiekBinsIMS, IMSRangeParams = getIMSRangeParams("mobility", 0.65, 0.75)))$mobility,
                 c(0.65, 0.75))
    expect_range(as.data.table(filter(ffPiekBinsIMS, IMSRangeParams = getIMSRangeParams("mobility", 0.0025, 0.003, TRUE)))[, mobility / mz],
                 c(0.0025, 0.003))
    expect_equal(filter(ffPiekBinsIMS, IMSRangeParams = getIMSRangeParams("mobility", 0, 10)), ffPiekBinsIMS)
    expect_error(filter(ffPiekBinsIMS, IMSRangeParams = getIMSRangeParams("CCS", 0, 1000)), "CCS") # no CCS calculated
})

test_that("basic usage", {
    expect_equal(nrow(as.data.table(ffOpenMS)), length(ffOpenMS))
    
    expect_identical(getFeatureQualityNames(ffOpenMSQ), featureQualityNames(group = FALSE))
    checkmate::expect_names(names(as.data.table(ffOpenMSQ)),
                            must.include = c(getFeatureQualityNames(ffOpenMSQ),
                                             getFeatureQualityNames(ffOpenMSQ, scores = TRUE)))
})

# XCMSImpXCMS <- getXCMSSet(ffXCMS)
XCMSImpXCMS3 <- getXCMSSet(ffXCMS3, loadRawData = FALSE)
XCMSImpOpenMS <- doExportXCMS(ffOpenMS, loadRawData = FALSE)
XCMSImpEP <- getXCMSSet(ffEP, loadRawData = FALSE)
XCMSImpKPIC2 <- getXCMSSet(ffKPIC2, loadRawData = FALSE)
XCMSImpSIRIUS <- getXCMSSet(ffSIRIUS, loadRawData = FALSE)
featMZs <- function(f) lapply(featureTable(f), "[[", "mz")
test_that("XCMS conversion", {
    # expect_equal(nrow(xcms::peaks(XCMSImpXCMS)), length(ffXCMS))
    expect_equal(nrow(xcms::peaks(XCMSImpXCMS3)), length(ffXCMS3))
    expect_equal(nrow(xcms::peaks(XCMSImpOpenMS)), length(getExpFeats(ffOpenMS)))
    expect_equal(nrow(xcms::peaks(XCMSImpEP)), length(ffEP))
    expect_equal(nrow(xcms::peaks(XCMSImpKPIC2)), length(ffKPIC2))
    expect_equal(nrow(xcms::peaks(XCMSImpSIRIUS)), length(ffSIRIUS))
    
    # expect_known_val(xcms::peaks(XCMSImpXCMS), "ff-xcms_import_xcms")
    expect_known_val(xcms::peaks(XCMSImpXCMS3), "ff-xcms_import_xcms3")
    expect_known_val(xcms::peaks(XCMSImpOpenMS), "ff-xcms_import_openms")
    expect_known_val(xcms::peaks(XCMSImpEP), "ff-xcms_import_ep")
    expect_known_val(xcms::peaks(XCMSImpKPIC2), "ff-xcms_import_kpic2")
    expect_known_val(xcms::peaks(XCMSImpSIRIUS), "ff-xcms_import_sirius")
    
    # expect_equal(featMZs(importFeatures(XCMSImpXCMS, "xcms", anaInfoOne)), featMZs(ffXCMS))
    expect_equal(featMZs(importFeatures(XCMSImpXCMS3, "xcms", anaInfoOne)), featMZs(ffXCMS3))
    expect_equal(featMZs(importFeatures(XCMSImpEP, "xcms", epAnaInfo)), featMZs(ffEP))
    expect_equal(featMZs(importFeatures(XCMSImpKPIC2, "xcms", anaInfoOne)), featMZs(ffKPIC2))
    expect_equal(featMZs(importFeatures(XCMSImpSIRIUS, "xcms", anaInfoOne)), featMZs(ffSIRIUS))
    
    expect_equal(nrow(xcms::peaks(getXCMSSet(ffPiekBinsIMS, loadRawData = FALSE, IMS = FALSE))), 0)
    expect_equal(nrow(xcms::peaks(getXCMSSet(ffPiekBinsIMS, loadRawData = FALSE, IMS = TRUE))), length(ffPiekBinsIMS))
})

# XCMS3ImpXCMS <- getXCMSnExp(ffXCMS, loadRawData = FALSE)
XCMS3ImpXCMS3 <- getXCMSnExp(ffXCMS3)
XCMS3ImpOpenMS <- doExportXCMS3(ffOpenMS, loadRawData = FALSE)
XCMS3ImpEP <- getXCMSnExp(ffEP, loadRawData = FALSE)
XCMS3ImpKPIC2 <- getXCMSnExp(ffKPIC2, loadRawData = FALSE)
XCMS3ImpSIRIUS <- getXCMSnExp(ffSIRIUS, loadRawData = FALSE)
test_that("XCMS3 conversion", {
    # expect_equal(nrow(xcms::chromPeaks(XCMS3ImpXCMS)), length(ffXCMS))
    expect_equal(nrow(xcms::chromPeaks(XCMS3ImpXCMS3)), length(ffXCMS3))
    expect_equal(nrow(xcms::chromPeaks(XCMS3ImpOpenMS)), length(getExpFeats(ffOpenMS)))
    expect_equal(nrow(xcms::chromPeaks(XCMS3ImpEP)), length(ffEP))
    expect_equal(nrow(xcms::chromPeaks(XCMS3ImpKPIC2)), length(ffKPIC2))
    expect_equal(nrow(xcms::chromPeaks(XCMS3ImpSIRIUS)), length(ffSIRIUS))
    
    # expect_known_val(xcms::chromPeaks(XCMS3ImpXCMS), "ff-xcms3_import_xcms")
    expect_known_val(xcms::chromPeaks(XCMS3ImpXCMS3), "ff-xcms3_import_xcms3")
    expect_known_val(xcms::chromPeaks(XCMS3ImpOpenMS), "ff-xcms3_import_openms")
    expect_known_val(xcms::chromPeaks(XCMS3ImpEP), "ff-xcms3_import_ep")
    expect_known_val(xcms::chromPeaks(XCMS3ImpKPIC2), "ff-xcms3_import_kpic2")
    expect_known_val(xcms::chromPeaks(XCMS3ImpSIRIUS), "ff-xcms3_import_sirius")
    
    # expect_equal(featMZs(importFeatures(XCMS3ImpXCMS, "xcms3", anaInfoOne)), featMZs(ffXCMS))
    expect_equal(featMZs(importFeatures(XCMS3ImpXCMS3, "xcms3", anaInfoOne)), featMZs(ffXCMS3))
    expect_equal(featMZs(importFeatures(XCMS3ImpEP, "xcms3", epAnaInfo)), featMZs(ffEP))
    expect_equal(featMZs(importFeatures(XCMS3ImpKPIC2, "xcms3", anaInfoOne)), featMZs(ffKPIC2))
    expect_equal(featMZs(importFeatures(XCMS3ImpSIRIUS, "xcms3", anaInfoOne)), featMZs(ffSIRIUS))
    
    # BUG: xcms::chromPeaks() errors with zero features
    expect_false(xcms::hasChromPeaks(getXCMSnExp(ffPiekBinsIMS, loadRawData = FALSE, IMS = FALSE)))
    expect_equal(nrow(xcms::chromPeaks(getXCMSnExp(ffPiekBinsIMS, loadRawData = FALSE, IMS = TRUE))), length(ffPiekBinsIMS))
})

KPIC2ImpKPIC2 <- getPICSet(ffKPIC2)
KPIC2ImpOpenMS <- getPICSet(ffOpenMS, loadRawData = FALSE)
test_that("KPIC2 conversion", {
    expect_equal(KPIC2ImpKPIC2[[1]]$peakinfo[, "mz"], ffKPIC2[[1]]$mz)
    expect_equal(KPIC2ImpOpenMS[[1]]$peakinfo[, "mz"], ffOpenMS[[1]]$mz)
})

test_that("importFeaturesTable", {
    ffImp <- importFeaturesTable(as.data.table(ffOpenMS), analysisInfo(ffOpenMS))
    adtImp <- as.data.table(ffImp)
    adtRef <- as.data.table(ffOpenMS)
    expect_equal(as.data.table(adtImp), as.data.table(adtRef)[, names(adtImp), with = FALSE])
    checkmate::expect_names(names(as.data.table(adtImp)), disjunct.from = "isocount")
    checkmate::expect_names(names(as.data.table(importFeaturesTable(as.data.table(ffOpenMS), analysisInfo(ffOpenMS), addCols = "isocount"))),
                            must.include = "isocount")
    
    setsCols <- c("set", "adduct", "ion_mz")    
    mandCols <- c("analysis", "ret", "mz", "intensity")
    calcCols <- c("ID", "retmin", "retmax", "mzmin", "mzmax", "area")
    checkmate::expect_names(names(as.data.table(importFeaturesTable(as.data.table(ffOpenMS)[, c(mandCols, setsCols), with = FALSE], analysisInfo(ffOpenMS)))),
                            must.include = c(mandCols, calcCols, setsCols))
    
    mandColsIMS <- c(mandCols, "mobility")
    calcColsIMS <- c(calcCols, "mobmin", "mobmax", "mob_intensity", "mob_area", "ims_precursor_ID", "mob_assign_method", "mob_reintegr_method")
    checkmate::expect_names(names(as.data.table(importFeaturesTable(as.data.table(ffPiekSuspIMS)[, mandColsIMS, with = FALSE], analysisInfo(ffPiekSuspIMS)))),
                            must.include = c(mandColsIMS, calcColsIMS))
    
    expect_warning(importFeaturesTable(as.data.table(ffOpenMS)[, intensity_rel := 1], analysisInfo(ffOpenMS), addCols = "intensity_rel"),
                   "not be considered: intensity_rel")
    expect_warning(importFeaturesTable(as.data.table(ffOpenMS)[, -"set"], analysisInfo(ffOpenMS)),
                   "Assuming this is not a sets workflow")
    expect_warning(importFeaturesTable(as.data.table(ffPiekSuspIMS)[, -"mobility"], analysisInfo(ffPiekSuspIMS)),
                   "invalid in non-IMS")
    expect_warning(importFeaturesTable(as.data.table(ffOpenMS)[analysis != analysisInfo(ffOpenMS)$analysis[1]], analysisInfo(ffOpenMS)),
                   paste("in the input data:", analysisInfo(ffOpenMS)$analysis[1]))
    expect_warning(importFeaturesTable(as.data.table(ffOpenMS), analysisInfo(ffOpenMS)[-1]),
                   paste("in analysisInfo:", analysisInfo(ffOpenMS)$analysis[1]))
    
    suppressWarnings(expect_setequal(analyses(importFeaturesTable(as.data.table(ffOpenMS)[analysis != analysisInfo(ffOpenMS)$analysis[1]], analysisInfo(ffOpenMS)[-2])),
                                     analyses(ffOpenMS)[-(1:2)]))
})

getModifiedAI <- function(ai)
{
    ret <- ffOpenMS
    analysisInfo(ret) <- ai
    return(ret)
}
expectModAI <- function(cb)
{
    ret <- ffOpenMS
    ai <- cb(copy(analysisInfo(ret)))
    analysisInfo(ret) <- ai
    expect_equal(analysisInfo(ret), ai)
}

test_that("anaInfo modification", {
    expect_equal(getModifiedAI(analysisInfo(ffOpenMS)), ffOpenMS)
    expectModAI(\(x) { x$har <- 1; x})
    expectModAI(\(x) x[seq(nrow(x), 1)])
    expect_error(expectModAI(\(x) { x$analysis[1] <- "nope"; x }), regexp = "Cannot modify analysis")
    expect_error(expectModAI(\(x) x[, -"replicate"]), regexp = "missing")
    expect_error(expectModAI(\(x) rbind(x, x)), regexp = "Cannot add analyses")
    expect_error(expectModAI(\(x) x[-1]), regexp = "Cannot remove analyses")
    
    expect_equal(ffOpenMS[seq(length(analyses(ffOpenMS)), 1)], ffOpenMS)
    checkmate::expect_names(analyses(ffOpenMS[seq(length(analyses(ffOpenMS)), 1), reorder = TRUE]),
                                     identical.to = rev(analyses(ffOpenMS)))
    checkmate::expect_names(names(featureTable(ffOpenMS[seq(length(analyses(ffOpenMS)), 1), reorder = TRUE])),
                            identical.to = rev(analyses(ffOpenMS)))

    expect_equal(xcms::peaks(getXCMSSet(ffXCMS3[seq(length(analyses(ffXCMS3)), 1), reorder = TRUE], loadRawData = FALSE)),
                 xcms::peaks(XCMSImpXCMS3)[order(xcms::peaks(XCMSImpXCMS3)[, "sample"], decreasing = TRUE), ])
    
    expect_equal(xcms::chromPeaks(getXCMSnExp(ffXCMS3[seq(length(analyses(ffXCMS3)), 1), reorder = TRUE])),
                 xcms::chromPeaks(XCMS3ImpXCMS3)[order(xcms::chromPeaks(XCMS3ImpXCMS3)[, "sample"], decreasing = TRUE), ])
    expect_equal(xcms::chromPeaks(getXCMSnExp(ffEP[seq(length(analyses(ffEP)), 1), reorder = TRUE], loadRawData = FALSE)),
                 xcms::chromPeaks(XCMS3ImpEP)[order(xcms::chromPeaks(XCMS3ImpEP)[, "sample"], decreasing = TRUE), ])
    
    expect_equal(getPICSet(ffKPIC2[seq(length(analyses(ffKPIC2)), 1), reorder = TRUE]), KPIC2ImpKPIC2[seq(length(KPIC2ImpKPIC2), 1)])
    expect_equal(getPICSet(ffOpenMS[seq(length(analyses(ffOpenMS)), 1), reorder = TRUE], loadRawData = FALSE),
                 KPIC2ImpOpenMS[seq(length(KPIC2ImpOpenMS), 1)])
})

test_that("Sets functionality", {
    # proper (de)neutralization
    expect_equal(patRoon:::calculateMasses(unset(ffOpenMS, "positive")[[1]]$mz, as.adduct("[M+H]+"), "neutral"),
                 ffOpenMS[[1]]$mz)
    expect_equal(patRoon:::calculateMasses(ffOpenMS[[1]]$mz, as.adduct("[M+H]+"), "mz"), ffOpenMS[[1]]$ion_mz)
    expect_equal(analysisInfo(unset(ffOpenMS, "positive"), TRUE), getTestAnaInfoPos())
    expect_equal(analysisInfo(ffOpenMS[, sets = "positive"])[, -"set"], getTestAnaInfoPos(), ignore_attr = TRUE)
    expect_equal(unique(ffOpenMS[[1]]$adduct), "[M+H]+")
    expect_equal(sets(filter(ffOpenMS, sets = "positive", negate = TRUE)), "negative")
    expect_length(ffOpenMS[, sets = character()], 0)
    expect_equal(sets(ffOpenMS[, sets = "positive"]), "positive")
    expect_equal(sets(ffOpenMS[, sets = rev(sets(ffOpenMS)), reorder = TRUE]), rev(sets(ffOpenMS)))
    expect_equal(sets(ffOpenMS[, ni = set == "positive"]), "positive")
    expect_length(makeSet(ffXCMS3, ffXCMS3[FALSE], adducts = "[M+H]+"), length(ffXCMS3))
    expect_length(makeSet(ffXCMS3[FALSE], ffXCMS3[FALSE], adducts = "[M+H]+"), 0)
})

anaInfoNS <- getTestAnaInfoNS()
ffNS <- getTestFeaturesNS(anaInfoNS)
XCMSImpOpenMSNS <- doExportXCMSNS(ffNS, loadRawData = FALSE)
XCMS3ImpOpenMSNS <- doExportXCMS3NS(ffNS, loadRawData = FALSE)

test_that("set unsupported functionality", {
    expect_equal(featMZs(importFeatures(XCMSImpOpenMSNS, "xcms", anaInfoNS)), featMZs(ffNS))    
    expect_equal(featMZs(importFeatures(XCMS3ImpOpenMSNS, "xcms3", anaInfoNS)), featMZs(ffNS))
})
