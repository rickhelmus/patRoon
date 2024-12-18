# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include mslibrary.R
NULL

#' Load MS library data from MassBank of North America (\acronym{MoNA})
#'
#' This function loads, verifies and curates MS library data from \href{https://mona.fiehnlab.ucdavis.edu/}{MoNA}
#' \file{.json} files.
#'
#' @templateVar algo an efficient \code{C++} JSON loader
#' @templateVar do load MS library data
#' @templateVar generic loadMSLibrary
#' @templateVar algoParam json
#' @template algo_generator
#'
#' @templateVar format JSON
#' @template loadMSLibrary
#'
#' @templateVar whatCP MS library
#' @template chemPropCalc
#'
#' @details This function uses \code{C++} with \CRANpkg{Rcpp} and \CRANpkg{rapidjsonr} to efficiently load and parse
#'   JSON files from \href{https://mona.fiehnlab.ucdavis.edu/}{MoNA}. An advantage compared to
#'   \code{\link{loadMSLibraryMSP}} is that this function supports loading spectral annotations.
#'   
#'   The record field names are converted to those used in \file{.msp} files.
#'
#' @export
loadMSLibraryMoNAJSON <- function(file, prefCalcChemProps = TRUE, neutralChemProps = FALSE, potAdducts = TRUE,
                                  potAdductsLib = TRUE, absMzDev = 0.002, calcSPLASH = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFileExists(file, "r", add = ac)
    aapply(checkmate::assertFlag, . ~ prefCalcChemProps + neutralChemProps + potAdductsLib + calcSPLASH,
           fixed = list(add = ac))
    checkmate::assert(checkmate::checkFlag(potAdducts),
                      checkmate::checkCharacter(potAdducts, any.missing = FALSE, min.chars = 1),
                      checkmate::checkList(potAdducts, types = c("adduct", "character"), any.missing = FALSE),
                      .var.name = "potAdducts")
    checkmate::assertNumber(absMzDev, lower = 0, finite = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    hash <- makeHash(makeFileHash(file), prefCalcChemProps, neutralChemProps, potAdducts, potAdductsLib, absMzDev,
                     calcSPLASH)
    cd <- loadCacheData("MSLibraryJSON", hash)
    if (!is.null(cd))
        return(cd)
    
    lib <- readMoNAJSON(normalizePath(file))
    lib <- sanitizeMSLibrary(lib, prefCalcChemProps, neutralChemProps, potAdducts, potAdductsLib, absMzDev, calcSPLASH)
    
    ret <- MSLibrary(records = lib$records[], spectra = lib$spectra, algorithm = "json")
    
    saveCacheData("MSLibraryJSON", ret, hash)
    
    return(ret)
}
