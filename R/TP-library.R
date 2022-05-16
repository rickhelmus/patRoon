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
#'   library. Valid options are: \code{"InChIKey"}, \code{"InChIKey1"}, \code{"InChI"}, \code{"SMILES"}.
#'
#' @templateVar parNULL TRUE
#' @template tp_gen-scr
#'
#' @template tp_gen-sim
#' @template fp-args
#'
#' @return The TPs are stored in an object derived from the \code{\link{transformationProductsStructure}} class.
#'
#' @section Custom TP libraries: The \code{TPLibrary} argument is used to specify a custom TP library. This should be a
#'   \code{data.frame} where each row specifies a TP for a parent, with the following columns: \itemize{
#'
#'   \item \code{parent_name} and \code{TP_name}: The name of the parent/TP.
#'
#'   \item \code{parent_SMILES} and \code{TP_SMILES} The \acronym{SMILES} of the parent structure.
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
    checkmate::assertChoice(matchParentsBy, c("InChIKey", "InChIKey1", "InChI", "SMILES"), null.ok = FALSE, add = ac)
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
    {
        parents <- getTPParents(parents, skipInvalid, prefCalcChemProps)
        
        # match with library
        if (matchParentsBy == "InChIKey1")
        {
            dataLib <- getIKBlock1(TPLibrary$parent_InChIKey)
            dataSusp <- getIKBlock1(parents$InChIKey)
        }
        else
        {
            dataLib <- TPLibrary[[paste0("parent_", matchParentsBy)]]
            dataSusp <- parents[[matchParentsBy]]
        }
        
        # rename from suspect list
        TPLibrary[, parent_name_lib := parent_name] # store original
        TPLibrary[, parent_name := parents[match(dataLib, dataSusp)]$name]
        
        # only take data in both
        dataInBoth <- intersect(dataLib, dataSusp)
        TPLibrary <- TPLibrary[dataLib %chin% dataInBoth]
        parents <- parents[dataSusp %chin% dataInBoth]
    }
    else
    {
        parents <- unique(TPLibrary[, grepl("^parent_", names(TPLibrary)), with = FALSE], by = "parent_name")
        setnames(parents, sub("^parent_", "", names(parents)))
    }

    results <- split(TPLibrary, by = "parent_name")
    
    curTPIDs <- setNames(vector("integer", length = length(results)), names(results))
    prepTPs <- function(r, pn, pid, gen, prvLogPDiff)
    {
        # remove parent columns
        set(r, j = grep("^parent_", names(r), value = TRUE), value = NULL)
        
        # remove TP_ prefix
        cols <- grep("^TP_", names(r), value = TRUE)
        setnames(r, cols, sub("^TP_", "", cols))
        setnames(r, "name", "name_lib")
        
        r[, retDir := 0] # may be changed below
        r[, generation := gen]
        
        r[, ID := curTPIDs[pn] + seq_len(nrow(r))]
        curTPIDs[pn] <<- curTPIDs[pn] + nrow(r)
        r[, parent_ID := pid]
        
        # make it additive so LogPDiff corresponds to the original parent
        if (!is.null(r[["LogPDiff"]]) && !is.null(prvLogPDiff))
            r[, LogPDiff := LogPDiff + prvLogPDiff]
        
        return(r)
    }
    results <- Map(results, names(results), f = prepTPs, MoreArgs = list(pid = NA_integer_, gen = 1, prvLogPDiff = NULL))
    
    if (generations > 1)
    {
        for (gen in seq(2, generations))
        {
            results <- Map(results, names(results), f = function(r, pn)
            {
                tps <- r[generation == (gen-1)]
                nexttps <- rbindlist(lapply(split(tps, seq_len(nrow(tps))), function(tpRow)
                {
                    nt <- copy(TPLibrary[parent_InChIKey == tpRow$InChIKey])
                    return(prepTPs(nt, pn, tpRow$ID, gen, tpRow$LogPDiff))
                }))
                return(rbind(r, nexttps))
            })
        }
    }
    
    # fill in chem IDs and names now that we sorted out all TPs
    results <- Map(results, names(results), f = function(r, pn)
    {
        set(r, j = "chem_ID", value = match(r$InChIKey, unique(r$InChIKey)))
        set(r, j = "name", value = paste0(pn, "-TP", r$chem_ID))
    })
    
    if (!is.null(TPLibrary[["parent_LogP"]]) && !is.null(TPLibrary[["TP_LogP"]]))
    {
        results <- Map(results, TPLibrary[match(names(results), parent_name)]$parent_LogP,
                       f = function(r, pLogP) set(r, j = "retDir", value = fifelse(r$LogP < pLogP, -1, 1)))
    }
    else if (!is.null(TPLibrary[["LogPDiff"]]))
    {
        results <- lapply(results, function(x) set(x, j = "retDir", value = fcase(x$LogPDiff < 0, -1,
                                                                                  x$LogPDiff > 0, 1,
                                                                                  default = 0)))
    }
    
    results <- pruneList(results, checkZeroRows = TRUE)
    parents <- parents[name %in% names(results)]
    results <- results[match(parents$name, names(results))] # sync order
    
    ret <- transformationProductsLibrary(calcSims = calcSims, fpType = fpType, fpSimMethod = fpSimMethod,
                                         parents = parents, products = results)
    saveCacheData("TPsLib", ret, hash)
    return(ret)
}
