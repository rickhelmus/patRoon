#' @include main.R
#' @include components.R
#' @include utils-components.R
#' @include feature_groups-screening.R
#' @include feature_groups-screening-set.R
NULL

genTPSpecSimilarities <- function(pl1, pl2, method, precDiff, removePrecursor,
                                  mzWeight, intWeight, absMzDev, relMinIntensity)
{
    getSim <- function(shift)
    {
        if (!is.null(pl2))
        {
            return(specSimilarity(pl1, pl2, method = method, shift = shift,
                                  precDiff = precDiff, removePrecursor = removePrecursor, mzWeight = mzWeight,
                                  intWeight = intWeight, absMzDev = absMzDev,
                                  relMinIntensity = relMinIntensity))
        }
        return(NA)
    }
    
    return(list(specSimilarity = getSim("none"), specSimilarityPrec = getSim("precursor"),
                specSimilarityBoth = getSim("both")))
}

doGenComponentsTPs <- function(fGroups, fGroupsTPs, pred, MSPeakLists, minRTDiff,
                               simMethod, removePrecursor, mzWeight, intWeight, absMzDev, relMinIntensity)
{
    if (length(fGroups) == 0)
        return(componentsTPs(componentInfo = data.table(), components = list()))
    
    hash <- makeHash(fGroups, fGroupsTPs, pred, MSPeakLists, minRTDiff, simMethod, removePrecursor,
                     mzWeight, intWeight, absMzDev, relMinIntensity)
    cd <- loadCacheData("componentsTPs", hash)
    if (!is.null(cd))
        return(cd)
    
    # for every precursor:
    #   - check for any matching fGroups (based on mass)
    #   - if none, skip
    #   - similarly, check which precursors are present (only mz)
    #   - for any fGroup that matches the precursor:
    #       - filter TPs (retention, intensity, ...)
    
    precFGMapping <- linkPrecursorsToFGroups(pred, fGroups)
    TPFGMapping <- linkTPsToFGroups(pred, fGroupsTPs)
    gInfoPrec <- groupInfo(fGroups); gInfoTPs <- groupInfo(fGroupsTPs)
    
    cat("Linking precursors and TPs ...\n")
    precCount <- length(names(pred))
    prog <- openProgBar(0, precCount)
    
    compTab <- rbindlist(Map(names(pred), predictions(pred), seq_len(precCount), f = function(pname, preds, i)
    {
        precFGs <- precFGMapping[name == pname][["group"]]
        TPs <- TPFGMapping[TP_name %in% preds$name]
        
        if (length(precFGs) == 0 || nrow(TPs) == 0)
            return(NULL)
        
        # limit columns a bit to not bloat components too much
        # UNDONE: column selection OK?
        predCols <- c("name", "InChIKey", "formula", "mass", "RTDir",
                      "reaction_add", "reaction_sub", "deltaMZ")
        preds <- preds[, intersect(names(preds), predCols), with = FALSE]
        
        comps <- rbindlist(sapply(precFGs, function(precFG)
        {
            # UNDONE: do more checks etc
            
            ret <- merge(TPs, preds, by.x = "TP_name", by.y = "name")
            
            # merge rows with duplicate fGroups, for instance, caused by different TPs with equal mass
            ret[, TP_name := paste0(TP_name, collapse = ","), by = "group"]
            ret <- unique(ret, by = "group")
            
            # dummy intensity value so e.g. plotSpec works            
            ret[, intensity := 1]
            
            ret[, c("ret", "mz") := gInfoTPs[group, c("rts", "mzs")]]
            ret[, mzDiff := gInfoPrec[precFG, "mzs"] - mz]
            
            if (minRTDiff > 0)
            {
                rtDiffs <- gInfoTPs[ret$group, "rts"] - gInfoPrec[precFG, "rts"]
                ret <- ret[RTDir == 0 | abs(rtDiffs) <= minRTDiff | (rtDiffs < 0 & RTDir < 0) | (rtDiffs > 0 & RTDir > 0)]
            }
            
            precMSMS <- MSPeakLists[[precFG]][["MSMS"]]
            if (nrow(ret) > 0 && !is.null(precMSMS))
            {
                sims <- rbindlist(Map(ret$group, ret$mzDiff, f = function(g, mzd)
                {
                    genTPSpecSimilarities(precMSMS, MSPeakLists[[g]][["MSMS"]], method = simMethod,
                                          precDiff = -mzd, removePrecursor = removePrecursor,
                                          mzWeight = mzWeight, intWeight = intWeight, absMzDev = absMzDev,
                                          relMinIntensity = relMinIntensity)
                }))
                ret <- cbind(ret, sims)
            }
            else
                ret[, c("specSimilarity", "specSimilarityPrec", "specSimilarityBoth") := NA]
            
            return(ret)
        }, simplify = FALSE), idcol = "precursor_group")
        
        setTxtProgressBar(prog, i)
        
        return(comps)
    }), idcol = "precursor_name")
    
    setTxtProgressBar(prog, precCount)
    close(prog)
    
    if (nrow(compTab) > 0)
    {
        compTab[, name := paste0("CMP", .GRP), by = c("precursor_name", "precursor_group")]
        compTab[, links := list(list(unique(name))), by = c("TP_name", "group")] # link to other components having this TP
        compTab[, links := mapply(links, name, FUN = setdiff)] # remove self-links
        
        compList <- split(compTab[, -"name"], by = c("precursor_name", "precursor_group"), keep.by = FALSE)
        
        compInfo <- unique(compTab[, c("name", "precursor_name", "precursor_group")])
        compInfo[, size := sapply(compList, nrow)]
        compInfo[, links := lapply(compList, function(cmp) unique(unlist(cmp$links)))] # overal links
        setcolorder(compInfo, "name")
        
        names(compList) <- compInfo$name
    }
    else
    {
        compInfo <- data.table()
        compList <- list()
    }
    
    ret <- componentsTPs(componentInfo = compInfo[], components = compList)
    saveCacheData("componentsTPs", ret, hash)
    
    return(ret)
}

#' @template components_noint
#' @export
componentsTPs <- setClass("componentsTPs", contains = "components")

setMethod("initialize", "componentsTPs",
          function(.Object, ...) callNextMethod(.Object, ..., algorithm = "tp"))


#' @describeIn componentsTPs Returns all component data in a table.
#' @export
setMethod("as.data.table", "componentsTPs", function(x)
{
    # as regular method, but get rid of double links column when merging
    # components / componentInfo
    
    x@componentInfo <- x@componentInfo[, -"links"]
    ret <- callNextMethod(x)
    
    return(ret)
})

#' @export
setMethod("filter", "componentsTPs", function(obj, ..., formulas = NULL, negate = FALSE)
{
    # UNDONE: if formulas is set, also remove fGroups without assignment? Otherwise document!
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(formulas, "formulas", null.ok = TRUE, add = ac)
    checkmate::assertFlag(negate, add = ac)
    checkmate::reportAssertions(ac)
    
    if (length(obj) == 0)
        return(obj)
    
    if (!is.null(formulas))
    {
        if (is.null(obj[[1]][["reaction_add"]]))
            stop("formula filter is only available for logic TP predictions")
        
        oldn <- length(obj)
        
        obj@components <- Map(componentTable(obj), componentInfo(obj)$precursor_group, f = function(cmp, gName)
        {
            # check if subtracting is possible, ie by checking if
            # subtraction doesn't lead to negative element counts
            canSub <- function(f, fg)
            {
                if (is.null(formulas[[fg]]) || length(f) == 0 || !nzchar(f))
                    return(TRUE) # UNDONE?
                candidateForms <- unique(formulas[[fg]]$neutral_formula)
                for (cf in candidateForms)
                {
                    fl <- splitFormulaToList(subtractFormula(cf, f))
                    if (all(fl >= 0))
                        return(TRUE)
                }
                return(FALSE)
            }
            if (negate)
                canSub <- Negate(canSub)
            
            # filter results where subtraction of any of the precursor formulas is impossible
            cmp <- cmp[!nzchar(reaction_sub) | sapply(reaction_sub, canSub, gName)]
            
            # filter results where addition is not part of TP candidate formulas
            if (nrow(cmp) > 0)
                cmp <- cmp[!nzchar(reaction_add) | mapply(reaction_add, group, FUN = canSub)]

            return(cmp)
        })
        
        obj@components <- pruneList(obj@components, checkZeroRows = TRUE)
        
        newn <- length(obj)
        printf("Done! Filtered %d (%.2f%%) components. Remaining: %d\n", oldn - newn, if (oldn == 0) 0 else (1-(newn/oldn))*100, newn)
    }
    
    if (...length() > 0)
        return(callNextMethod(obj, ..., negate = negate))
    
    return(obj)
})

#' @describeIn componentsTPs Plots an interactive network graph for linked
#'   components. Components are linked when (partial) overlap occurs of their
#'   containing transformation products. The graph is constructed with the
#'   \pkg{\link{igraph}} package and rendered with \pkg{\link{visNetwork}}.
#'
#' @param obj The \code{componentsTPs} object to plot.
#'
#' @template plotGraph
#'
#' @export
setMethod("plotGraph", "componentsTPs", function(obj, onlyLinked)
{
    checkmate::assertFlag(onlyLinked)

    cInfo <- componentInfo(obj)
    cTable <- componentTable(obj)
    
    titles <- sprintf("<b>%s</b> (%s - %s)<br>TPs: <i>%s (%s)</i>",
                      names(obj), cInfo$precursor_name, cInfo$precursor_group,
                      sapply(cTable, function(cmp) paste0(cmp$TP_name, collapse = ", ")),
                      sapply(cTable, function(cmp) paste0(unique(cmp$group), collapse = ", ")))
    makeGraph(obj, onlyLinked, titles)
})

#' @export
setMethod("generateComponentsTPs", "featureGroups", function(fGroups, fGroupsTPs = fGroups, pred,
                                                             MSPeakLists, minRTDiff = 20, simMethod,
                                                             removePrecursor = FALSE, mzWeight = 0,
                                                             intWeight = 1, absMzDev = 0.005,
                                                             relMinIntensity = 0.1)
{
    # UNDONE: optionally remove TPs with equal formula as parent (how likely are these?)
    # UNDONE: optional MSPeakLists? and MSPeakListsTPs?
    # UNDONE: doc when fGroups/fGroupsTPs needs to be fGroupsScreening
    
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertClass, . ~ fGroupsTPs + pred + MSPeakLists,
           c("featureGroups", "TPPredictions", "MSPeakLists"), fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ minRTDiff + mzWeight + intWeight + absMzDev + relMinIntensity,
           lower = 0, finite = TRUE, fixed = list(add = ac))
    checkmate::assertChoice(simMethod, c("cosine", "jaccard"), add = ac)
    checkmate::assertFlag(removePrecursor, add = ac)
    checkmate::reportAssertions(ac)
    
    if (needsScreening(pred) &&
        (!inherits(fGroups, "featureGroupsScreening") || !inherits(fGroupsTPs, "featureGroupsScreening")))
        stop("Input feature groups need to be screened for (TP) suspects!")

    return(doGenComponentsTPs(fGroups, fGroupsTPs, pred, MSPeakLists, minRTDiff, simMethod, removePrecursor,
                              mzWeight, intWeight, absMzDev, relMinIntensity))
})

#' @export
setMethod("generateComponentsTPs", "featureGroupsSet", function(fGroups, fGroupsTPs = fGroups, pred,
                                                                MSPeakLists, minRTDiff = 20,
                                                                simMethod, removePrecursor = FALSE, mzWeight = 0,
                                                                intWeight = 1, absMzDev = 0.005,
                                                                relMinIntensity = 0.1)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertClass, . ~ fGroupsTPs + pred + MSPeakLists,
           c("featureGroupsSet", "TPPredictions", "MSPeakListsSet"), fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ minRTDiff + mzWeight + intWeight + absMzDev + relMinIntensity,
           lower = 0, finite = TRUE, fixed = list(add = ac))
    checkmate::assertChoice(simMethod, c("cosine", "jaccard"), add = ac)
    checkmate::assertFlag(removePrecursor, add = ac)
    checkmate::reportAssertions(ac)

    if (needsScreening(pred) &&
        (!inherits(fGroups, "featureGroupsScreeningSet") || !inherits(fGroupsTPs, "featureGroupsScreeningSet")))
        stop("Input feature groups need to be screened for (TP) suspects!")
    ret <- doGenComponentsTPs(fGroups, fGroupsTPs, pred, MSPeakLists, minRTDiff, simMethod, removePrecursor,
                              mzWeight, intWeight, absMzDev, relMinIntensity)
    
    # UNDONE: more efficient method to get set specific fGroups?
    gNamesTPsSets <- sapply(sets(fGroupsTPs), function(s) names(fGroupsTPs[, sets = s]), simplify = FALSE)
    ionizedMSPeaksLists <- sapply(sets(MSPeakLists), ionize, obj = MSPeakLists, simplify = FALSE)
    ret@components <- Map(ret@components, ret@componentInfo$precursor_group, f = function(cmp, precFG)
    {
        for (s in sets(fGroupsTPs))
        {
            # mark set presence
            set(cmp, j = s, value = cmp$group %in% gNamesTPsSets[[s]])
            
            # calculate per set spectrum similarities
            simColNames <- paste0(c("specSimilarity", "specSimilarityPrec", "specSimilarityBoth"), "-", s)
            # if (any(!is.na(cmp$specSimilarity))) browser()
            precMSMS <- ionizedMSPeaksLists[[s]][[precFG]][["MSMS"]]
            if (!is.null(precMSMS))
            {
                sims <- rbindlist(Map(cmp$group, cmp$mzDiff, f = function(g, mzd)
                {
                    genTPSpecSimilarities(precMSMS, ionizedMSPeaksLists[[s]][[g]][["MSMS"]], method = simMethod,
                                          precDiff = -mzd, removePrecursor = removePrecursor,
                                          mzWeight = mzWeight, intWeight = intWeight, absMzDev = absMzDev,
                                          relMinIntensity = relMinIntensity)
                }))
                cmp <- cbind(cmp, setnames(sims, simColNames))
            }
            else
                cmp[, (simColNames) := NA]
        }
        
        # move spec similarity columns to end
        setcolorder(cmp, setdiff(names(cmp), grep("^specSimilarity", names(cmp), value = TRUE)))
        
        return(cmp)
    })
    
    return(ret)
})
