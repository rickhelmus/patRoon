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
#' @param generations An \code{integer} that specifies the number of transformation generations. TPs for subsequent
#'   iterations obtained by repeating the library search where the TPs from the previous generation are considered
#'   parents.
#' @param matchParentsBy A \code{character} that specifies how the input parents are matched with the data from the TP
#'   library. Valid options are: \code{"InChIKey"}, \code{"InChIKey1"}, \code{"InChI"}, \code{"SMILES"},
#'   \code{"formula"}, \code{"name"}. If the parent from the TP library is matched with multiple input parents then only
#'   the first is considered.
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
#' @section Custom TP libraries: The \code{TPLibrary} argument is used to specify a custom TP library. This should be a
#'   \code{data.frame} where each row specifies a TP for a parent, with the following columns: \itemize{
#'
#'   \item \code{parent_name} and \code{TP_name}: The name of the parent/TP.
#'
#'   \item \code{parent_SMILES} and \code{TP_SMILES} The \acronym{SMILES} of the parent/TP structure.
#'
#'   \item \code{parent_LogP} and \code{TP_LogP} The \code{log P} values for the parent/TP. (\strong{optional})
#'
#'   \item \code{LogPDiff} The difference between parent and TP \code{Log P} values. Ignored if \emph{both}
#'   \code{parent_LogP} and \code{TP_LogP} are specified. (\strong{optional})
#'
#'   }
#'
#'   Other columns are allowed, and will be included in the final object. Multiple TPs for a single parent are specified
#'   by repeating the value within \code{parent_} columns.
#'
#' @export
generateTPsLibrary <- function(parents = NULL, TPLibrary = NULL, generations = 1, skipInvalid = TRUE,
                               prefCalcChemProps = TRUE, matchParentsBy = "InChIKey", calcSims = FALSE,
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
    aapply(checkmate::assertFlag, . ~ skipInvalid + prefCalcChemProps + calcSims, fixed = list(add = ac))
    checkmate::assertChoice(matchParentsBy, c("InChIKey", "InChIKey1", "InChI", "SMILES", "formula", "name"),
                            null.ok = FALSE, add = ac)
    aapply(checkmate::assertString, . ~ fpType + fpSimMethod, min.chars = 1, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    hash <- makeHash(parents, TPLibrary, generations, skipInvalid, prefCalcChemProps, matchParentsBy, calcSims,
                     fpType, fpSimMethod)
    cd <- loadCacheData("TPsLib", hash)
    if (!is.null(cd))
        return(cd)
    
    if (is.null(TPLibrary))
        TPLibrary <- copy(PubChemTransformations) # default to embedded PC transformations
    else
    {
        TPLibrary <- copy(as.data.table(TPLibrary))
        
        # add chem infos where necessary
        for (wh in c("parent", "TP"))
        {
            for (col in c("formula", "InChI", "InChIKey"))
            {
                whcol <- paste0(wh, "_", col)
                if (is.null(TPLibrary[[whcol]]))
                {
                    whSMI <- paste0(wh, "_SMILES")
                    TPLibrary[, (whcol) := switch(col,
                                                  formula = babelConvert(get(whSMI), "smi", "formula"),
                                                  InChI = babelConvert(get(whSMI), "smi", "inchi"),
                                                  InChIKey = babelConvert(get(whSMI), "smi", "inchikey"))]
                }
            }
            whmcol <- paste0(wh, "_neutralMass")
            if (is.null(TPLibrary[[whmcol]]))
                TPLibrary[, (whmcol) := sapply(get(paste0(wh, "_formula")), getFormulaMass)]
        }
    }

    if (!is.null(parents))
        parents <- getTPParents(parents, skipInvalid, prefCalcChemProps)
    
    prep <- prepareParentsForLib(parents, TPLibrary, matchParentsBy)
    parents <- prep$parents; TPLibrary <- prep$TPLibrary

    results <- getProductsFromLib(TPLibrary, generations)
    parents <- parents[name %in% names(results)]
    results <- results[match(parents$name, names(results))] # sync order
    
    ret <- transformationProductsLibrary(calcSims = calcSims, fpType = fpType, fpSimMethod = fpSimMethod,
                                         parents = parents, products = results)
    saveCacheData("TPsLib", ret, hash)
    return(ret)
}
