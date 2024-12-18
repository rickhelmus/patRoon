# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
#' @include TP-formula.R
NULL

#' @rdname transformationProductsFormula-class
transformationProductsLibraryFormula <- setClass("transformationProductsLibraryFormula",
                                                 contains = "transformationProductsFormula")

setMethod("initialize", "transformationProductsLibraryFormula",
          function(.Object, ...) callNextMethod(.Object, algorithm = "library_formula", ...))


#' Obtain transformation products (TPs) from a library with formula data
#'
#' Automatically obtains transformation products from a library with formula data.
#'
#' @templateVar algo a library
#' @templateVar do obtain transformation products
#' @templateVar generic generateTPs
#' @templateVar algoParam library_formula
#' @template algo_generator
#'
#' @details This function is similar to \code{\link{generateTPsLibrary}}, however, it only require formula information
#'   of the parent and TPs.
#'
#' @param parents The parents for which transformation products should be obtained. This should be either a suspect list
#'   (see \link[=suspect-screening]{suspect screening} for more information) or the resulting output of
#'   \code{\link{screenSuspects}}. The suspect (hits) are used as parents. If \code{NULL} then TPs for all parents in
#'   the library are obtained.
#' @param TPLibrary A \code{data.frame}. See the details below.
#'
#' @templateVar id formula
#' @template tp_lib
#'
#' @templateVar whatCP parent suspect list
#' @template chemPropCalc
#'
#' @inheritParams generateTPsLibrary
#'
#' @return The TPs are stored in an object derived from the \code{\link{transformationProductsFormula}} class.
#'
#' @note Unlike \code{\link{generateTPsLibrary}}, this function defaults the \code{matchParentsBy} and
#'   \code{matchGenerationsBy} arguments to \code{"name"}. While matching by \code{formula} is also possible, it is
#'   likely that duplicate parent formulae (\emph{i.e.} isomers) are present in \code{parents} and/or \code{TPLibrary},
#'   making matching by formula unsuitable. However, if you are sure that no duplicate formulae are present, it may be
#'   better to set the matching method to \code{"formula"}.
#'
#' @seealso \code{\link{generateTPsLibrary}} to generate TPs from a library that contains structural information.
#' @seealso \code{\link{genFormulaTPLibrary}} to automatically generate formula TP libraries.
#'
#' @export
generateTPsLibraryFormula <- function(parents = NULL, TPLibrary, generations = 1, skipInvalid = TRUE,
                                      prefCalcChemProps = TRUE, neutralChemProps = FALSE, matchParentsBy = "name",
                                      matchGenerationsBy = "name")
{
    # NOTE: this is mainly a simplified version of generateTPsLibrary
    
    checkmate::assert(
        checkmate::checkNull(parents),
        checkmate::checkClass(parents, "data.frame"),
        checkmate::checkClass(parents, "featureGroupsScreening"),
        checkmate::checkClass(parents, "featureGroupsScreeningSet"),
        .var.name = "parents"
    )
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertDataFrame(TPLibrary, any.missing = FALSE, min.rows = 1, add = ac)
    libCols <- c("parent_name", "parent_formula", "TP_name", "TP_formula")
    assertHasNames(TPLibrary, libCols, add = ac)
    for (cl in libCols)
        assertListVal(TPLibrary, cl, checkmate::assertCharacter, min.chars = 1, any.missing = FALSE, add = ac)
    
    if (is.data.frame(parents))
        assertSuspectList(parents, needsAdduct = FALSE, skipInvalid = TRUE, add = ac)
    checkmate::assertCount(generations, positive = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ skipInvalid + prefCalcChemProps + neutralChemProps, fixed = list(add = ac))
    aapply(checkmate::assertChoice, . ~ matchParentsBy + matchGenerationsBy, null.ok = FALSE,
           fixed = list(choices = c("InChIKey", "InChIKey1", "InChI", "SMILES", "formula", "name"), add = ac))
    checkmate::reportAssertions(ac)
    
    hash <- makeHash(parents, TPLibrary, generations, skipInvalid, prefCalcChemProps, neutralChemProps, matchParentsBy,
                     matchGenerationsBy)
    cd <- loadCacheData("TPsLibFormula", hash)
    if (!is.null(cd))
        return(cd)
    
    if (!is.null(parents))
        parents <- getTPParents(parents, skipInvalid, prefCalcChemProps, neutralChemProps, checkWhat = "formula")
    
    prep <- prepareDataForTPLibrary(parents, TPLibrary, generations, matchParentsBy, matchGenerationsBy, "formula",
                                    FALSE)
    
    ret <- transformationProductsLibraryFormula(parents = prep$parents, products = prep$products)
    saveCacheData("TPsLibFormula", ret, hash)
    return(ret)
}

#' Automatically generate a transformation product library with formula data.
#'
#' Functionality to automatically generate a TP library with formula data from a set of transformation rules, which can
#' be used with \code{\link{generateTPsLibraryFormula}}. TP calculation will be skipped if the transformation involves
#' subtraction of elements not present in the parent.
#'
#' @param parents The parents to which the given transformation rules should be used to generate the TP library. Should
#'   be either a suspect list (see \link[=suspect-screening]{suspect screening} for more information) or the resulting
#'   output of \code{\link{screenSuspects}}.
#' @param minMass The minimum mass for a TP to be kept.
#' @param generations An \code{integer} that specifies the number of transformation generations that should be
#'   calculated. If \code{generations>1} then TPs are calculated by applying the transformation rules to the TPs
#'   generated in the previous generation.
#' @param skipInvalid Set to \code{TRUE} to skip parents without formula information. Otherwise an error is thrown.
#'
#' @template tp_trans
#'
#' @templateVar whatCP parent suspect list
#' @template chemPropCalc
#'
#' @return A \code{data.table} that is suitable for the \code{TPLibrary} argument to
#'   \code{\link{generateTPsLibraryFormula}}.
#'
#' @seealso \code{\link{generateTPsLibraryFormula}} and \code{\link{generateTPsLogic}}
#'
#' @export
genFormulaTPLibrary <- function(parents, transformations = NULL, minMass = 40, generations = 1, skipInvalid = TRUE,
                                prefCalcChemProps = TRUE, neutralChemProps = FALSE)
{
    checkmate::assert(
        checkmate::checkClass(parents, "data.frame"),
        checkmate::checkClass(parents, "featureGroupsScreening"),
        checkmate::checkClass(parents, "featureGroupsScreeningSet"),
        .var.name = "parents"
    )
    
    ac <- checkmate::makeAssertCollection()
    if (is.data.frame(parents))
        assertSuspectList(parents, needsAdduct = FALSE, skipInvalid = TRUE, add = ac)
    assertLogicTransformations(transformations, null.ok = TRUE, add = ac)
    checkmate::assertNumber(minMass, lower = 0, finite = TRUE, add = ac)
    checkmate::assertCount(generations, positive = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ skipInvalid + prefCalcChemProps + neutralChemProps, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    parents <- getTPParents(parents, skipInvalid, prefCalcChemProps, neutralChemProps, checkWhat = "formula")
    transformations <- getTPLogicTransformations(transformations)
    
    genLibItems <- function(parNames, parFormulas, parNeutralMasses, gen)
    {
        rbindlist(Map(parNames, parFormulas, parNeutralMasses, f = function(parName, parForm, parMass)
        {
            tab <- data.table(parent_name = parName, parent_formula = parForm, parent_neutralMass = parMass,
                              TP_name = paste0(parName, "-", transformations$transformation),
                              TP_formula = mapply(transformations$add, transformations$sub, FUN = function(a, s)
                              {
                                  return(subtractFormula(addFormula(parForm, a), s))
                              }),
                              TP_neutralMass = parMass + transformations$deltaMZ,
                              generation = gen,
                              retDir = transformations$retDir)
            return(tab[TP_neutralMass >= minMass & sapply(TP_formula, function(f) min(splitFormulaToList(f)) > 0)])
        }))
    }
    
    ret <- genLibItems(parents$name, parents$formula, parents$neutralMass, 1)
    
    if (generations > 1)
    {
        for (gen in seq(2, generations))
        {
            p <- ret[generation == (gen-1)]
            ret <- rbind(ret, genLibItems(p$TP_name, p$TP_formula, p$TP_neutralMass, gen))
        }
    }
    
    return(ret)
}
