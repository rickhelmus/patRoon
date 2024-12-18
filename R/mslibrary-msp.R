# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include mslibrary.R
NULL

#' Load MS library data from MSP files
#'
#' This function loads, verifies and curates MS library data from MSP files.
#'
#' @templateVar algo an efficient \code{C++} MSP loader
#' @templateVar do load MS library data
#' @templateVar generic loadMSLibrary
#' @templateVar algoParam msp
#' @template algo_generator
#'
#' @templateVar format MSP
#' @template loadMSLibrary
#'
#' @templateVar whatCP MS library
#' @template chemPropCalc
#'
#' @details This function uses \code{C++} with \CRANpkg{Rcpp} to efficiently load and parse MSP files, and is mainly
#'   optimized for loading the \file{.msp} files from \href{https://massbank.eu/MassBank/}{MassBank EU} and
#'   \href{https://mona.fiehnlab.ucdavis.edu/}{MoNA}. Files from other sources may also work, any feedback on this is
#'   welcome!
#'
#' @param parseComments If \code{TRUE} then comments in the file are parsed to obtain additional fields, such as
#'   \acronym{SMILES}, \code{PubChemCID} and \code{Resolution}. Note that some records specify this data either in the
#'   comments or as a regular field, hence, to ensure that loaded data is most complete it is recommend to set
#'   \code{parseComments=TRUE}.
#'
#' @note The mass spectrum parser currently only supports space separated entries (MSP formerly also allows other
#'   formats).
#'
#' @export
loadMSLibraryMSP <- function(file, parseComments = TRUE, prefCalcChemProps = TRUE, neutralChemProps = FALSE,
                             potAdducts = TRUE, potAdductsLib = TRUE, absMzDev = 0.002, calcSPLASH = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFileExists(file, "r", add = ac)
    aapply(checkmate::assertFlag, . ~ parseComments + prefCalcChemProps + neutralChemProps + potAdductsLib + calcSPLASH,
           fixed = list(add = ac))
    checkmate::assert(checkmate::checkFlag(potAdducts),
                      checkmate::checkCharacter(potAdducts, any.missing = FALSE, min.chars = 1),
                      checkmate::checkList(potAdducts, types = c("adduct", "character"), any.missing = FALSE),
                      .var.name = "potAdducts")
    checkmate::assertNumber(absMzDev, lower = 0, finite = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    hash <- makeHash(makeFileHash(file), parseComments, prefCalcChemProps, neutralChemProps, potAdducts, potAdductsLib,
                     absMzDev, calcSPLASH)
    cd <- loadCacheData("MSLibraryMSP", hash)
    if (!is.null(cd))
        return(cd)
    
    lib <- readMSP(normalizePath(file), parseComments)
    lib <- sanitizeMSLibrary(lib, prefCalcChemProps, neutralChemProps, potAdducts, potAdductsLib, absMzDev, calcSPLASH)
    
    ret <- MSLibrary(records = lib$records[], spectra = lib$spectra, algorithm = "msp")
    
    saveCacheData("MSLibraryMSP", ret, hash)
    
    return(ret)
}
