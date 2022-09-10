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
#'   (see \link[=suspect-screening]{suspect screening} for more information) or the resulting output
#'   \code{\link{screenSuspects}}. The suspect (hits) are used as parents. If \code{NULL} then TPs for all parents in
#'   the library are obtained.
#' @param TPLibrary A \code{data.frame}. See the details below.

#' @inheritParams generateTPsLibrary
#'
#' @templateVar whatCP parent suspect list
#' @template chemPropCalc
#'
#' @return The TPs are stored in an object derived from the \code{\link{transformationProductsFormula}} class.
#'
#' @section TP libraries: The \code{TPLibrary} argument is used to specify a custom TP library. This should be a
#'   \code{data.frame} where each row specifies a TP for a parent, with the following columns: \itemize{
#'
#'   \item \code{parent_name} and \code{TP_name}: The name of the parent/TP.
#'
#'   \item \code{parent_formula} and \code{TP_formula} The chemical formula of the parent/TP.
#'
#'   }
#'
#'   Other columns are allowed, and will be included in the final object. Multiple TPs for a single parent are specified
#'   by repeating the value within \code{parent_} columns.
#'
#' @seealso \code{\link{generateTPsLibrary}} to generate TPs from a library that contains structural information.
#'
#' @export
generateTPsLibraryFormula <- function(parents = NULL, TPLibrary, generations = 1, skipInvalid = TRUE,
                                      prefCalcChemProps = TRUE)
{
    # NOTE: this is mainly a simplified version of generateTPsLibrary
    
    checkmate::assert(
        checkmate::checkNull(parents),
        checkmate::checkClass(parents, "data.frame"),
        checkmate::checkClass(parents, "compounds"),
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
    aapply(checkmate::assertFlag, . ~ skipInvalid + prefCalcChemProps, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    hash <- makeHash(parents, TPLibrary, generations, skipInvalid, prefCalcChemProps)
    cd <- loadCacheData("TPsLibFormula", hash)
    if (!is.null(cd))
        return(cd)
    
    TPLibrary <- copy(as.data.table(TPLibrary))
    
    # add chem infos where necessary
    for (wh in c("parent", "TP"))
    {
        whmcol <- paste0(wh, "_neutralMass")
        if (is.null(TPLibrary[[whmcol]]))
            TPLibrary[, (whmcol) := sapply(get(paste0(wh, "_formula")), getFormulaMass)]
    }

    if (!is.null(parents))
    {
        parents <- getTPParents(parents, skipInvalid, prefCalcChemProps)

        # only take data in both
        dataInBoth <- intersect(TPLibrary$parent_name, parents$name)
        TPLibrary <- TPLibrary[parent_name %chin% dataInBoth]
        parents <- parents[name %chin% dataInBoth]
    }
    else
    {
        parents <- unique(TPLibrary[, grepl("^parent_", names(TPLibrary)), with = FALSE], by = "parent_name")
        setnames(parents, sub("^parent_", "", names(parents)))
    }
    
    results <- getProductsFromLib(TPLibrary, generations)
    parents <- parents[name %in% names(results)]
    results <- results[match(parents$name, names(results))] # sync order
    
    ret <- transformationProductsLibraryFormula(parents
                                                = parents, products = results)
    saveCacheData("TPsLibFormula", ret, hash)
    return(ret)
}
