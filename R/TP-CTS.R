#' @include main.R
#' @include TP-structure.R
NULL

#' @export
transformationProductsCTS <- setClass("transformationProductsCTS", contains = "transformationProductsStructure")

setMethod("initialize", "transformationProductsCTS",
          function(.Object, ...) callNextMethod(.Object, algorithm = "cts", ...))

# NOTE: this function is called by a withProg() block, so handles progression updates
runCTS <- function(parentRow, transLibrary, generationLimit, errorRetries, calcXLogP)
{
    CTSURL <- "https://qed.epa.gov/cts/rest/metabolizer/run"
    # CTSURL <- "https://qed.epacdx.net/cts/rest/metabolizer/run" # older version??
    
    resp <- httr::RETRY("POST", url = CTSURL, encode = "json", times = errorRetries,
                        body = list(structure = parentRow$SMILES, generationLimit = generationLimit,
                                    transformationLibraries = list(transLibrary)))
    
    httr::warn_for_status(resp)
    
    # HACK: do here to simplify NULL returns below
    doProgress()
    
    if (httr::status_code(resp) != 200)
        return(NULL)
    
    cont <- httr::content(resp, "parsed")
    if (is.null(cont) || is.null(cont[["data"]]) || cont$data$total_products == 0)
        return(NULL)

    curTPID <- 0
    processChilds <- function(chi, parentID = 1)
    {
        curTPID <<- curTPID + 1
        res <- c(list(ID = curTPID, parentID = parentID), chi$data)
        if (length(chi$children) > 0)
        {
            sub <- lapply(chi$children, processChilds, parentID = curTPID-1) # NOTE: -1 since we just incremented
            return(rbindlist(c(list(res), sub)))
        }
        return(as.data.table(res))
    }
    
    ret <- rbindlist(lapply(cont$data$data$children, processChilds))
    setnames(ret, "smiles", "SMILES")
    
    # Assign some unique identifier
    ret[, name := paste0(parentRow$name, "-TP", seq_len(nrow(ret)))]
    
    if (calcXLogP)
    {
        ret[, XLogP := calculateXLogP(SMILES, mustWork = FALSE)]
        ret[, retDir := fifelse(XLogP < parentRow$XLogP, -1, 1)]
    }
    else
        ret[, retDir := 0]
    
    return(ret)
}


#' @export
generateTPsCTS <- function(parents, transLibrary, generationLimit = 1, errorRetries = 3, skipInvalid = TRUE,
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
    aapply(checkmate::assertCount, . ~ generationLimit + errorRetries, positive = TRUE, fixed = list(add = ac))
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
        
        baseHash <- makeHash(transLibrary, generationLimit, errorRetries, skipInvalid, fpType, fpSimMethod)
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
                                                         transLibrary, generationLimit, errorRetries, calcXLogP,
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
