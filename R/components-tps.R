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
    
    hash <- makeHash(fGroups, rGroupsIn, rGroupsEff)
    cd <- loadCacheData("componentsTPs", hash)
    if (!is.null(cd))
        return(cd)
    

    # for every precursor:
    #   - check for any matching fGroups (based on mass)
    #   - if none, skip
    #   - similarly, check which precursors are present (only mz)
    #   - for any fGroup that matches the precursor:
    #       - filter TPs (retention, intensity, ...)
    
    
    suspList <- convertToSuspects(pred, adduct)
    
    printf("Screening precursors and TPs...\n")
    screening <- screenSuspects(fGroups, suspList, mzWindow = mzWindow)
    screening <- screening[!is.na(group)] # remove non-hits
    susps <- suspects(pred)
    screeningPrec <- screening[name %in% susps$name]
    screeningTPs <- screening[!name %in% susps$name]

    if (!is.null(rGroupsIn))
    {
        fg <- replicateGroupSubtract(fGroups, rGroupsEff, inThreshold)
        screeningPrec <- screeningPrec[group %in% names(fg)]
    }
    if (!is.null(rGroupsEff))
    {
        fg <- replicateGroupSubtract(fGroups, rGroupsIn, effThreshold)
        screeningTPs <- screeningTPs[group %in% names(fg)]
    }
    
    compTab <- rbindlist(mapply(susps$name, predictions(pred), SIMPLIFY = FALSE, FUN = function(pname, preds)
    {
        scrP <- screeningPrec[screeningPrec$name == pname]
        scrTP <- screeningTPs[name %in% preds$Identifier]
        
        if (nrow(scrP) == 0 || nrow(scrTP) == 0)
            return(NULL)
        
        comps <- rbindlist(lapply(split(scrP, by = "group"), function(scrRow)
        {
            # UNDONE: do more checks etc
            
            ret <- merge(scrTP, preds, by.x = "name", by.y = "Identifier")
            setnames(ret, c("name", "group"), c("TP_name", "TP_group"))
            return(ret)
        }), idcol = "precursor_group")
    }), idcol = "precursor_susp_name")

    compList <- split(compTab, by = c("precursor_susp_name", "precursor_group"), keep.by = FALSE)
    names(compList) <- paste0("CMP", seq_along(compList))

    compInfo <- unique(compTab[, c("precursor_susp_name", "precursor_group")])
    compInfo[, name := names(compList)]
    compInfo[, size := sapply(compList, nrow)]

    ret <- componentsTPs(componentInfo = compInfo, components = compList)    
    saveCacheData("componentsTPs", ret, hash)
    
    return(ret)
}
