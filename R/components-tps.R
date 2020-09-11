#' @include main.R
#' @include components.R
#' @include utils-components.R
NULL

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
                      names(obj), cInfo$precursor_susp_name, cInfo$precursor_group,
                      sapply(cTable, function(cmp) paste0(cmp$TP_name, collapse = ", ")),
                      sapply(cTable, function(cmp) paste0(unique(cmp$group), collapse = ", ")))
    makeGraph(obj, onlyLinked, titles)
})

#' @export
generateComponentsTPs <- function(fGroups, fGroupsTPs = NULL, pred, MSPeakLists, mzWindow = 0.005, minRTDiff = 20)
{
    # UNDONE: make method
    # UNDONE: optionally remove TPs with equal formula as parent (how likely are these?)
    # UNDONE: optional MSPeakLists? and MSPeakListsTPs?
    
    
    
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertClass, . ~ fGroups + fGroupsTPs + pred + MSPeakLists,
           c("featureGroupsScreening", "featureGroupsScreening", "TPPredictions", "MSPeakLists"),
           null.ok = c(FALSE, TRUE, FALSE, FALSE), fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ mzWindow + minRTDiff,
           lower = 0, finite = TRUE, fixed = list(add = ac))
    checkmate::reportAssertions(ac)

    if (length(fGroups) == 0)
        return(componentsTPs(componentInfo = data.table(), components = list()))

    if (is.null(fGroupsTPs))
        fGroupsTPs <- fGroups
    
    hash <- makeHash(fGroups, fGroupsTPs, pred, MSPeakLists, mzWindow, minRTDiff)
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
    gInfoPrec <- groupInfo(fGroups); gInfoTPs <- groupInfo(fGroupsTPs)
    screeningTPs <- screenInfo(fGroupsTPs)
    
    compTab <- rbindlist(mapply(names(pred), predictions(pred), SIMPLIFY = FALSE, FUN = function(pname, preds)
    {
        precFGs <- precFGMapping[name == pname][["group"]]
        scrTP <- screeningTPs[name %in% preds$name]
        
        if (length(precFGs) == 0 || nrow(scrTP) == 0)
            return(NULL)
        
        rnCols <- c("name", "mz")
        scrTP <- scrTP[, c("group", rnCols), with = FALSE]
        setnames(scrTP, rnCols, c("TP_name", "TP_mz"))
        
        # limit columns a bit to not bloat components too much
        # UNDONE: column selection OK?
        predCols <- c("name", "InChIKey", "formula", "mass", "RTDir")
        preds <- preds[, intersect(names(preds), predCols), with = FALSE]
        
        return(rbindlist(lapply(precFGs, function(precFG)
        {
            # UNDONE: do more checks etc
            
            ret <- merge(scrTP, preds, by.x = "TP_name", by.y = "name")
            
            # merge rows with duplicate fGroups, for instance, caused by different TPs with equal mass
            ret[, TP_name := paste0(TP_name, collapse = ","), by = "group"]
            ret <- unique(ret, by = "group")
            
            # dummy intensity value so e.g. plotSpec works            
            ret[, intensity := 1]
            
            if (minRTDiff > 0)
            {
                rtDiffs <- gInfoTPs[ret$group, "rts"] - gInfoPrec[precFG, "rts"]
                ret <- ret[RTDir == 0 | abs(rtDiffs) <= minRTDiff | (rtDiffs < 0 & RTDir < 0) | (rtDiffs > 0 & RTDir > 0)]
            }
            
            precMSMS <- MSPeakLists[[precFG]][["MSMS"]]
            if (!is.null(precMSMS))
            {
                getSim <- function(g, shift)
                {
                    TPMSMS <- MSPeakLists[[g]][["MSMS"]]
                    if (!is.null(TPMSMS))
                    {
                        # UNDONE: make some arguments configurable (eg omit precursor)
                        
                        if (shift)
                            shift <- diff(gInfoPrec[c(precFG, g), "mzs"])
                        else
                            shift <- 0
                        return(spectrumSimilarity(MSPeakLists, precFG, g, MSLevel = 2, doPlot = FALSE, mzShift = shift))
                    }
                    
                    return(NA)
                }
                
                ret[, specSimilarity := sapply(group, getSim, shift = FALSE)]
                ret[, specSimilarityShift := sapply(group, getSim, shift = TRUE)]
            }
            else
                ret[, c("specSimilarity", "specSimilarityShift") := NA]
            
            return(ret)
        }), idcol = "precursor_group"))
    }), idcol = "precursor_susp_name")
    
    if (nrow(compTab) > 0)
    {
        compTab[, name := paste0("CMP", .GRP), by = c("precursor_susp_name", "precursor_group")]
        compTab[, links := list(list(unique(name))), by = c("TP_name", "group")] # link to other components having this TP
        compTab[, links := mapply(links, name, FUN = setdiff)] # remove self-links
        
        compList <- split(compTab[, -"name"], by = c("precursor_susp_name", "precursor_group"), keep.by = FALSE)
        
        compInfo <- unique(compTab[, c("name", "precursor_susp_name", "precursor_group")])
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
