context("MS library")

peakCount <- function(msl) sum(sapply(spectra(msl), nrow))

mslibraryMSP <- loadMSLibrary(file.path(getTestDataPathGeneric(), "MoNA-export-CASMI_2012.msp"), "msp")
mslibraryJSON <- loadMSLibrary(file.path(getTestDataPathGeneric(), "MoNA-export-CASMI_2016.json"), "json")
emptyFile <- tempfile(); file.create(emptyFile)

test_that("verify library loading", {
    expect_known_value(mslibraryMSP, testFile("mslibrary-msp"))
    expect_known_value(mslibraryJSON, testFile("mslibrary-json"))
    expect_known_show(mslibraryMSP, testFile("mslibrary-msp", text = TRUE))
    expect_known_show(mslibraryJSON, testFile("mslibrary-json", text = TRUE))

    expect_length(loadMSLibrary(emptyFile, "msp"), 0)
    expect_length(loadMSLibrary(emptyFile, "json"), 0)
})

mslibraryEmpty <- mslibraryMSP["none"]
test_that("basic usage",{
    expect_equal(names(mslibraryMSP[1:3]), names(mslibraryMSP)[1:3])
    expect_equal(names(mslibraryMSP[names(mslibraryMSP)[1:3]]), names(mslibraryMSP)[1:3])
    expect_equal(names(mslibraryMSP[c(TRUE, FALSE)]), names(mslibraryMSP)[c(TRUE, FALSE)])
    
    expect_length(mslibraryEmpty, 0)
    expect_length(mslibraryEmpty[1:5], 0)
    
    expect_equal(mslibraryMSP[[1]], spectra(mslibraryMSP)[[1]])
    expect_equal(mslibraryMSP[[names(mslibraryMSP)[1]]], spectra(mslibraryMSP)[[1]])
    expect_equal(callDollar(mslibraryMSP, names(mslibraryMSP)[2]), mslibraryMSP[[2]])

    expect_equal(nrow(records(mslibraryMSP)), length(mslibraryMSP))
    expect_equal(names(mslibraryMSP), records(mslibraryMSP)$DB_ID)
    expect_equal(nrow(as.data.table(mslibraryMSP)), peakCount(mslibraryMSP))
    expect_equal(nrow(as.data.table(mslibraryEmpty)), 0)
})

test_that("delete and filter", {
    checkmate::expect_names(names(delete(mslibraryMSP, i = 1)), disjunct.from = names(mslibraryMSP)[1])
    checkmate::expect_names(names(delete(mslibraryMSP, i = names(mslibraryMSP)[1])), disjunct.from = names(mslibraryMSP)[1])
    expect_length(delete(mslibraryMSP, i = names(mslibraryMSP)), 0)
    expect_length(delete(mslibraryMSP, i = 1:5), length(mslibraryMSP) - 5)
    expect_equal(records(delete(mslibraryMSP, i = 1:5)), records(mslibraryMSP)[-(1:5)])
    expect_length(delete(mslibraryMSP, i = function(...) 1:5), length(mslibraryMSP) - 5)
    expect_equal(delete(mslibraryMSP, j = 1)[[1]]$mz[1], mslibraryMSP[[1]]$mz[2])
    expect_equal(peakCount(delete(mslibraryMSP, j = 1:2)), peakCount(mslibraryMSP) - (length(mslibraryMSP) * 2))
    expect_false(delete(mslibraryMSP, j = function(...) 1)[[1]]$mz[1] == mslibraryMSP[[1]]$mz[1])
    expect_equal(peakCount(delete(mslibraryMSP, j = function(...) 1:2)), peakCount(mslibraryMSP) - (length(mslibraryMSP) * 2))
    expect_length(delete(mslibraryMSP, j = function(...) TRUE), 0)
    expect_equal(delete(mslibraryMSP, i = character()), mslibraryMSP)
    expect_equal(delete(mslibraryMSP, j = integer()), mslibraryMSP)
    expect_length(delete(mslibraryMSP), 0)
    expect_length(delete(mslibraryEmpty), 0)
    
    expect_true(all(records(filter(mslibraryMSP, properties = list(Instrument_type = "LC-ESI-QTOF")))$Instrument_type == "LC-ESI-QTOF"))
    expect_false(any(records(filter(mslibraryMSP,
                                    properties = list(Instrument_type = "LC-ESI-QTOF"),
                                    negate = TRUE))$Instrument_type == "LC-ESI-QTOF"))
    expect_range(records(filter(mslibraryMSP, massRange = c(100, 200)))$neutralMass, c(100, 200))
    expect_gt(min(records(filter(mslibraryMSP, massRange = c(0, 200), negate = TRUE))$neutralMass), 200)
    expect_range(as.data.table(filter(mslibraryMSP, mzRangeSpec = c(100, 200)))$mz, c(100, 200))
    expect_gt(min(as.data.table(filter(mslibraryMSP, mzRangeSpec = c(0, 200), negate = TRUE))$mz), 200)
    expect_gte(min(as.data.table(filter(mslibraryMSP, relMinIntensity = 0.1))$intensity), 10)
    expect_lt(max(as.data.table(filter(mslibraryMSP, relMinIntensity = 0.1, negate = TRUE))$intensity), 10)
    expect_lte(max(sapply(spectra(filter(mslibraryMSP, topMost = 3)), nrow)), 3)
    expect_lte(max(sapply(spectra(filter(mslibraryMSP, topMost = 3, negate = TRUE)), nrow)), 3)
    expect_length(filter(mslibraryMSP, onlyAnnotated = TRUE), 0)
    expect_equal(filter(mslibraryMSP, onlyAnnotated = TRUE, negate = TRUE), mslibraryMSP)
    expect_true(all(nzchar(as.data.table(filter(mslibraryJSON, onlyAnnotated = TRUE))$annotation)))
    # HACK: merge MSP/JSON to get library containing records with and without annotations
    expect_null(as.data.table(filter(merge(mslibraryMSP, mslibraryJSON), onlyAnnotated = TRUE, negate = TRUE))[["annotation"]])
})
