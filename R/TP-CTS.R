#' @include main.R
#' @include TP.R
NULL

#' @export
transformationProductsCTS <- setClass("transformationProductsCTS", contains = "transformationProducts")

setMethod("initialize", "transformationProductsCTS",
          function(.Object, ...) callNextMethod(.Object, algorithm = "cts", ...))

# NOTE: this function is called by a withProg() block, so handles progression updates
runCTS <- function(SMILES, transLibrary, generationLimit)
{
    CTSURL <- "https://qed.epa.gov/cts/rest/metabolizer/run"
    # CTSURL <- "https://qed.epacdx.net/cts/rest/metabolizer/run" # older version??
    
    # UNDONE: args for parameters and retries
    resp <- httr::RETRY("POST", url = CTSURL, encode = "json",
                        body = list(structure = SMILES, generationLimit = generationLimit,
                                    transformationLibraries = list(transLibrary)))
    
    # HACK: do here to simplify NULL returns below
    doProgress()
    
    if (httr::status_code(resp) != 200)
        return(NULL)
    
    cont <- httr::content(resp, "parsed")
    if (cont$data$total_products == 0)
        return(NULL)

    curTPID <- 0
    processChilds <- function(chi, parentID = 0)
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
    return(ret)
}


#' @export
generateTPsCTS <- function(parents, transLibrary, generationLimit = 1, skipInvalid = TRUE, fpType = "extended",
                           fpSimMethod = "tanimoto", parallel = TRUE)
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
    checkmate::assertCount(generationLimit, positive = TRUE, add = ac)
    checkmate::assertFlag(skipInvalid, add = ac)
    aapply(checkmate::assertString, . ~ fpType + fpSimMethod, min.chars = 1, fixed = list(add = ac))
    checkmate::assertFlag(parallel, add = ac)
    checkmate::reportAssertions(ac)
    
    parents <- getTPParents(parents, skipInvalid)
    
    baseHash <- makeHash(skipInvalid, fpType, fpSimMethod)
    setHash <- makeHash(parents, baseHash)

    # UNDONE: patRoon:::
    lapfunc <- if (parallel) future.apply::future_lapply else lapply
    results <- withProg(nrow(parents), parallel, do.call(lapfunc, list(parents$SMILES, runCTS, transLibrary,
                                                                       generationLimit)))
    names(results) <- parents$name

    results <- pruneList(results, checkZeroRows = TRUE)
    parents <- parents[name %in% names(results)]
    
    return(transformationProductsBT(parents = parents, products = results))
}
