#' @include main.R
#' @include components.R
NULL

#' @export
componentsTPs <- setClass("componentsTPs", contains = "components")

setMethod("initialize", "componentsTPs",
          function(.Object, ...) callNextMethod(.Object, ..., algorithm = "tp"))


#' @export
generateComponentsTPs <- function(fGroups, pred, adduct, mzWindow = 0.005, rGroupsIn = NULL, rGroupsEff = NULL,
                                  inThreshold = 0, effThreshold = 0)
{
    checkmate::assertClass(fGroups, "featureGroups")

    rGroups <- replicateGroups(fGroups)
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(pred, "TPPredictions")
    
    aapply(checkmate::assertNumber, . ~ mzWindow + inThreshold, effThreshold, lower = 0, finite = TRUE, fixed = list(add = ac))
    
    if (!is.null(rGroupsIn))
        checkmate::assertSubset(rGroupsIn, empty.ok = FALSE, choices = rGroups, add = ac)
    if (!is.null(rGroupsEff))
        checkmate::assertSubset(rGroupsEff, empty.ok = FALSE, choices = rGroups, add = ac)
    
    checkmate::reportAssertions(ac)
    
    adduct <- checkAndToAdduct(adduct)
    
    if (length(fGroups) == 0)
        return(componentsTPs(componentInfo = data.table(), components = list()))
    
    hash <- makeHash(fGroups, pred, adduct, mzWindow, rGroupsIn, rGroupsEff, inThreshold, effThreshold)
    cd <- loadCacheData("componentsTPs", hash)
    if (!is.null(cd))
        return(cd)
    

    # for every precursor:
    #   - check for any matching fGroups (based on mass)
    #   - if none, skip
    #   - similarly, check which precursors are present (only mz)
    #   - for any fGroup that matches the precursor:
    #       - filter TPs (retention, intensity, ...)
    
    
    suspList <- convertToSuspects(pred, adduct, includePrec = FALSE)

    printf("Screening TPs...\n")
    
    screeningTPs <- screenSuspects(fGroups, suspList, mzWindow = mzWindow)
    screeningTPs <- screeningTPs[!is.na(group)] # remove non-hits
    susps <- suspects(pred)
    precGroups <- linkPrecursorsToFGroups(pred, fGroups, adduct, mzWindow)
    
    if (!is.null(rGroupsIn))
    {
        # UNDONE: or isolate rGroup?
        fg <- replicateGroupSubtract(fGroups, rGroupsEff, inThreshold)
        precGroups <- precGroups[group %in% names(fg)]
    }
    if (!is.null(rGroupsEff))
    {
        fg <- replicateGroupSubtract(fGroups, rGroupsIn, effThreshold)
        screeningTPs <- screeningTPs[group %in% names(fg)]
    }
    
    compTab <- rbindlist(mapply(susps$name, predictions(pred), SIMPLIFY = FALSE, FUN = function(pname, preds)
    {
        scrP <- precGroups[precGroups$name == pname]
        scrTP <- screeningTPs[name %in% preds$Identifier]
        
        if (nrow(scrP) == 0 || nrow(scrTP) == 0)
            return(NULL)
        
        comps <- rbindlist(lapply(split(scrP, by = "group"), function(scrRow)
        {
            # UNDONE: do more checks etc
            # UNDONE: omit precursor columns? Or move to compInfo?
            # UNDONE: make column names more clear (e.g. from suspects)
            
            ret <- merge(scrTP, preds, by.x = "name", by.y = "Identifier")
            setnames(ret, c("name", "group"), c("TP_name", "TP_group"))
            return(ret)
        }), idcol = "precursor_group")
    }), idcol = "precursor_susp_name")

    
    compTab[, links := list(list(unique(name))), by = c("TP_name", "TP_group")]
    
    browser()
    
    compList <- split(compTab, by = c("precursor_susp_name", "precursor_group"), keep.by = FALSE)
    names(compList) <- paste0("CMP", seq_along(compList))

    compInfo <- unique(compTab[, c("precursor_susp_name", "precursor_group")])
    compInfo[, name := names(compList)]
    compInfo[, size := sapply(compList, nrow)]
    setcolorder(compInfo, "name")

    ret <- componentsTPs(componentInfo = compInfo, components = compList)    
    saveCacheData("componentsTPs", ret, hash)
    
    return(ret)
}
