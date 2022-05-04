#' @include main.R
#' @include TP-structure.R
NULL

#' @rdname transformationProductsStructure-class
transformationProductsCTS <- setClass("transformationProductsCTS", contains = "transformationProductsStructure")

setMethod("initialize", "transformationProductsCTS",
          function(.Object, ...) callNextMethod(.Object, algorithm = "cts", ...))

# NOTE: this function is called by a withProg() block, so handles progression updates
runCTS <- function(parentRow, transLibrary, generations, errorRetries, calcXLogP)
{
    CTSURL <- "https://qed.epa.gov/cts/rest/metabolizer/run"
    # CTSURL <- "https://qed.epacdx.net/cts/rest/metabolizer/run" # older version??
    
    resp <- httr::RETRY("POST", url = CTSURL, encode = "json", times = errorRetries,
                        body = list(structure = parentRow$SMILES, generationLimit = generations,
                                    transformationLibraries = list(transLibrary)))
    
    httr::warn_for_status(resp)
    
    # HACK: do here to simplify NULL returns below
    doProgress()
    
    if (httr::status_code(resp) != 200)
        return(NULL)
    
    cont <- httr::content(resp, "parsed")
    if (is.null(cont) || is.null(cont[["data"]]) || cont$data$total_products == 0)
        return(NULL)

    curTPID <- 0L
    processChilds <- function(chi, parent_ID)
    {
        curTPID <<- curTPID + 1L
        res <- c(list(ID = curTPID, parent_ID = parent_ID), chi$data)
        if (length(chi$children) > 0)
        {
            sub <- lapply(chi$children, processChilds, parent_ID = res$ID) # NOTE: this will alter curTPID, so don't use it for parent_ID!
            return(rbindlist(c(list(res), sub)))
        }
        return(as.data.table(res))
    }
    
    ret <- rbindlist(lapply(cont$data$data$children, processChilds, parent_ID = NA_integer_))
    
    setnames(ret, c("smiles", "routes"), c("SMILES", "transformation"))
    
    ret[, c("InChI", "InChIKey", "formula") := .(babelConvert(SMILES, "smi", "inchi", mustWork = TRUE),
                                                 babelConvert(SMILES, "smi", "inchikey", mustWork = TRUE),
                                                 convertToFormulaBabel(SMILES, "smi", mustWork = TRUE))]
    ret[, neutralMass := sapply(formula, getFormulaMass)]
    
    if (calcXLogP)
    {
        ret[, XLogP := calculateXLogP(SMILES, mustWork = FALSE)]
        ret[, retDir := fifelse(XLogP < parentRow$XLogP, -1, 1)]
    }
    else
        ret[, retDir := 0]
    
    # convert row IDs to unique IDs
    ret[, chem_ID := match(InChIKey, unique(InChIKey))]
    
    # Assign some unique identifier
    ret[, name := paste0(parentRow$name, "-TP", chem_ID)]
    
    setcolorder(ret, c("name", "ID", "parent_ID", "chem_ID", "SMILES", "InChI", "InChIKey", "formula", "neutralMass"))
    
    return(ret)
}


#' Obtain transformation products (TPs) with Chemical Transformation Simulator (CTS)
#'
#' Uses \href{https://qed.epa.gov/cts/}{Chemical Transformation Simulator (CTS)} to predict TPs.
#'
#' @templateVar algo CTS
#' @templateVar do obtain transformation products
#' @templateVar generic generateTPs
#' @templateVar algoParam cts
#' @template algo_generator
#'
#' @details This function uses the \CRANpkg{httr} package to access the Web API of CTS for automatic TP prediction.
#'   Hence, an Internet connection is mandatory. Please take care to not 'abuse' the CTS servers, \emph{e.g.} by running
#'   very large batch calculations in parallel, as this may result in rejected connections.
#'
#' @param transLibrary A \code{character} specifying which transformation library should be used. Currently supported
#'   are: \code{"hydrolysis"}, \code{"abiotic_reduction"}, \code{"photolysis_unranked"}, \code{"photolysis_ranked"},
#'   \code{"mammalian_metabolism"}, \code{"combined_abioticreduction_hydrolysis"},
#'   \code{"combined_photolysis_abiotic_hydrolysis"}.
#' @param generations An \code{integer} that specifies the number of transformation generations to predict.
#' @param errorRetries The maximum number of connection retries. Sets the \code{times} argument to the
#'   \code{\link[httr:RETRY]{http::RETRY}} function.
#' @param calcXLogP If \code{TRUE} then \code{XLogP} values will be calculated of parent and TPs to predict their
#'   retention order (\code{retDir}). The calculations are done with \CRANpkg{rcdk}.
#'
#' @return The TPs are stored in an object derived from the \code{\link{transformationProductsStructure}} class.
#'
#' @template parallel-arg
#' @template tp_gen-scr
#' @template tp_gen-sim
#' @template fp-args
#'
#' @seealso The website: \url{https://qed.epa.gov/cts/} and the
#'   \href{https://www.epa.gov/chemical-research/users-guide-chemical-transformation-simulator-cts}{CTS User guide}.
#'
#' @references \insertRef{Wolfe2016}{patRoon} \cr\cr \insertRef{TebesStevens2017}{patRoon} \cr\cr
#'   \insertRef{Yuan2020}{patRoon} \cr\cr \insertRef{Yuan2021}{patRoon}
#'
#' @export
generateTPsCTS <- function(parents, transLibrary, generations = 1, errorRetries = 3, skipInvalid = TRUE,
                           calcXLogP = TRUE, calcSims = FALSE, fpType = "extended", fpSimMethod = "tanimoto",
                           parallel = TRUE)
{
    checkmate::assert(
        checkmate::checkClass(parents, "data.frame"),
        checkmate::checkClass(parents, "compounds"),
        checkmate::checkClass(parents, "featureGroupsScreening"),
        checkmate::checkClass(parents, "featureGroupsScreeningSet"),
        .var.name = "parents"
    )
    
    ac <- checkmate::makeAssertCollection()
    if (is.data.frame(parents))
        assertSuspectList(parents, needsAdduct = FALSE, skipInvalid = TRUE, add = ac)
    checkmate::assertChoice(transLibrary, c("hydrolysis", "abiotic_reduction", "photolysis_unranked", 
                                            "photolysis_ranked", "mammalian_metabolism",
                                            "combined_abioticreduction_hydrolysis", 
                                            "combined_photolysis_abiotic_hydrolysis"), add = ac)
    aapply(checkmate::assertCount, . ~ generations + errorRetries, positive = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertFlag, . ~ skipInvalid + calcXLogP + calcSims + parallel, fixed = list(add = ac))
    aapply(checkmate::assertString, . ~ fpType + fpSimMethod, min.chars = 1, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    parents <- getTPParents(parents, skipInvalid)
    
    if (nrow(parents) == 0)
        results <- list()
    else
    {
        cacheDB <- openCacheDBScope()
        
        if (calcXLogP)
        {
            ph <- makeHash(parents)
            cd <- loadCacheData("TPsCTSXLogP", ph, cacheDB)
            if (is.null(cd))
            {
                printf("Calculating parent XLogP values... ")
                parents[, XLogP := calculateXLogP(SMILES, mustWork = FALSE)]
                printf("Done!\n")
                saveCacheData("TPsCTSXLogP", parents, ph, cacheDB)
            }
            else
                parents <- cd
        }
        
        parsSplit <- split(parents, seq_len(nrow(parents)))
        names(parsSplit) <- parents$name
        
        baseHash <- makeHash(transLibrary, generations, errorRetries, skipInvalid, fpType, fpSimMethod)
        setHash <- makeHash(parents, baseHash)
        cachedSet <- loadCacheSet("TPsCTS", setHash, cacheDB)
        hashes <- sapply(parsSplit, function(par) makeHash(baseHash, par[, c("name", "SMILES")], with = FALSE))
        cachedResults <- pruneList(sapply(hashes, function(h)
        {
            result <- NULL
            if (!is.null(cachedSet))
                result <- cachedSet[[h]]
            if (is.null(result))
                result <- loadCacheData("TPsCTS", h, cacheDB)
            return(result)
        }, simplify = FALSE))
        
        parsTBD <- setdiff(parents$name, names(cachedResults))
        newResults <- list()
        if (length(parsTBD) > 0)
        {
            lapfunc <- if (parallel) future.apply::future_sapply else sapply
            newResults <- withProg(length(parsTBD), parallel,
                                   do.call(lapfunc, list(parsSplit[parsTBD], patRoon:::runCTS,
                                                         transLibrary, generations, errorRetries, calcXLogP,
                                                         simplify = FALSE)))

            for (pn in names(newResults))
            {
                if (!is.null(newResults[[pn]]))
                    saveCacheData("TPsCTS", newResults[[pn]], hashes[[pn]], cacheDB)
            }
                        
            newResults <- pruneList(newResults, checkZeroRows = TRUE)
        }
        
        if (is.null(cachedSet))
            saveCacheSet("TPsCTS", hashes, setHash, cacheDB)
        
        results <- c(cachedResults, newResults)
        results <- results[intersect(parents$name, names(results))]
        parents <- parents[name %in% names(results)]
    }
    
    return(transformationProductsCTS(calcSims = calcSims, fpType = fpType, fpSimMethod = fpSimMethod,
                                     parents = parents, products = results))
}
