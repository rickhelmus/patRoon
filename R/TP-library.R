# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
#' @include TP-structure.R
NULL

#' @rdname transformationProductsStructure-class
transformationProductsLibrary <- setClass("transformationProductsLibrary", contains = "transformationProductsStructure")

setMethod("initialize", "transformationProductsLibrary",
          function(.Object, ...) callNextMethod(.Object, algorithm = "library", ...))


#' Obtain transformation products (TPs) from a library
#'
#' Automatically obtains transformation products from a library.
#'
#' @templateVar algo a library
#' @templateVar do obtain transformation products
#' @templateVar generic generateTPs
#' @templateVar algoParam library
#' @template algo_generator
#'
#' @details By default, a library is used that is based on data from
#'   \href{https://doi.org/10.5281/zenodo.5644560}{PubChem}. However, it also possible to use your own library.
#'
#' @param TPLibrary If \code{NULL}, a default \href{https://doi.org/10.5281/zenodo.5644560}{PubChem} based library is
#'   used. Otherwise, \code{TPLibrary} should be a \code{data.frame}. See the details below.
#'
#' @templateVar id SMILES
#' @template tp_lib
#'
#' @templateVar parNULL TRUE
#' @template tp_gen-scr
#'
#' @template tp_gen-sim
#' @template fp-args
#'
#' @templateVar whatCP parent suspect list
#' @template chemPropCalc
#'
#' @return The TPs are stored in an object derived from the \code{\link{transformationProductsStructure}} class.
#'
#' @export
generateTPsLibrary <- function(parents = NULL, TPLibrary = NULL, generations = 1, skipInvalid = TRUE,
                               prefCalcChemProps = TRUE, neutralChemProps = FALSE, neutralizeTPs = FALSE,
                               matchParentsBy = "InChIKey", matchGenerationsBy = "InChIKey", calcSims = FALSE,
                               fpType = "extended", fpSimMethod = "tanimoto")
{
    # UNDONE: default match by IK or IK1?
    
    checkmate::assert(
        checkmate::checkNull(parents),
        checkmate::checkClass(parents, "data.frame"),
        checkmate::checkClass(parents, "compounds"),
        checkmate::checkClass(parents, "featureGroupsScreening"),
        checkmate::checkClass(parents, "featureGroupsScreeningSet"),
        .var.name = "parents"
    )
    
    ac <- checkmate::makeAssertCollection()
    if (!is.null(TPLibrary))
    {
        checkmate::assertDataFrame(TPLibrary, any.missing = FALSE, min.rows = 1, add = ac)
        libCols <- c("parent_name", "parent_SMILES", "TP_name", "TP_SMILES")
        assertHasNames(TPLibrary, libCols, add = ac)
        for (cl in libCols)
            assertListVal(TPLibrary, cl, checkmate::assertCharacter, min.chars = 1, any.missing = FALSE, add = ac)
    }
        
    if (is.data.frame(parents))
        assertSuspectList(parents, needsAdduct = FALSE, skipInvalid = TRUE, add = ac)
    checkmate::assertCount(generations, positive = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ skipInvalid + prefCalcChemProps + neutralChemProps + neutralizeTPs + calcSims,
           fixed = list(add = ac))
    aapply(checkmate::assertChoice, . ~ matchParentsBy + matchGenerationsBy, null.ok = FALSE,
           fixed = list(choices = c("InChIKey", "InChIKey1", "InChI", "SMILES", "formula", "name"), add = ac))
    aapply(checkmate::assertString, . ~ fpType + fpSimMethod, min.chars = 1, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    hash <- makeHash(parents, TPLibrary, generations, skipInvalid, prefCalcChemProps, neutralChemProps, neutralizeTPs,
                     matchParentsBy, matchGenerationsBy, calcSims, fpType, fpSimMethod)
    cd <- loadCacheData("TPsLib", hash)
    if (!is.null(cd))
        return(cd)
    
    if (is.null(TPLibrary))
        TPLibrary <- copy(PubChemTransformations) # default to embedded PC transformations

    if (!is.null(parents))
        parents <- getTPParents(parents, skipInvalid, prefCalcChemProps, neutralChemProps)
    
    prep <- prepareDataForTPLibrary(parents, TPLibrary, generations, matchParentsBy, matchGenerationsBy, "InChIKey",
                                    neutralizeTPs)

    ret <- transformationProductsLibrary(calcSims = calcSims, fpType = fpType, fpSimMethod = fpSimMethod,
                                         parents = prep$parents, products = prep$products)
    saveCacheData("TPsLib", ret, hash)
    return(ret)
}
