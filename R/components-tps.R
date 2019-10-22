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
generateComponentsTPs <- function(fGroups, pred, MSPeakLists, adduct, mzWindow = 0.005,
                                  fGroupsPrec = fGroups, fGroupsTPs = fGroups, minRTDiff = 20)
{
    # UNDONE: optionally remove TPs with equal formula as parent (how likely are these?)
    
    aapply(checkmate::assertClass, . ~ fGroups + fGroupsPrec + fGroupsTPs, "featureGroups")

    rGroups <- replicateGroups(fGroups)
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(pred, "TPPredictions")
    checkmate::assertClass(MSPeakLists, "MSPeakLists")
    
    aapply(checkmate::assertNumber, . ~ mzWindow + minRTDiff,
           lower = 0, finite = TRUE, fixed = list(add = ac))
    
    checkmate::reportAssertions(ac)
    
    adduct <- checkAndToAdduct(adduct)
    
    if (length(fGroups) == 0)
        return(componentsTPs(componentInfo = data.table(), components = list()))
    
    hash <- makeHash(fGroups, pred, adduct, mzWindow, fGroupsPrec, fGroupsTPs)
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
    
    screeningTPs <- screenSuspects(fGroupsTPs, suspList, mzWindow = mzWindow)
    screeningTPs <- screeningTPs[!is.na(group)] # remove non-hits
    susps <- suspects(pred)
    precGroups <- linkPrecursorsToFGroups(pred, fGroupsPrec, adduct, mzWindow)
    
    gInfo <- groupInfo(fGroups)
    
    compTab <- rbindlist(mapply(susps$name, predictions(pred), SIMPLIFY = FALSE, FUN = function(pname, preds)
    {
        scrP <- precGroups[precGroups$name == pname]
        scrTP <- screeningTPs[name %in% preds$name]
        
        if (nrow(scrP) == 0 || nrow(scrTP) == 0)
            return(NULL)

        rnCols <- c("name", "mz", "exp_rt", "exp_mz")
        scrTP <- scrTP[, c("group", rnCols), with = FALSE]
        setnames(scrTP, rnCols, c("TP_name", "TP_mz", "rt", "mz"))
        
        # limit columns a bit to not bloat components too much
        # UNDONE: column selection OK?
        prCols <- c("name", "InChIKey", "formula", "mass", "RTDir")
        preds <- preds[, intersect(names(preds), prCols), with = FALSE]

        comps <- rbindlist(lapply(split(scrP, by = "group"), function(scrRow)
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
                rtDiffs <- ret$rt - gInfo[scrRow$group, "rts"]
                ret <- ret[RTDir == 0 | abs(rtDiffs) <= minRTDiff | (rtDiffs < 0 & RTDir < 0) | (rtDiffs > 0 & RTDir > 0)]
            }

            precMSMS <- MSPeakLists[[scrRow$group]][["MSMS"]]
            if (!is.null(precMSMS))
            {
                getSim <- function(g, shift)
                {
                    TPMSMS <- MSPeakLists[[g]][["MSMS"]]
                    if (!is.null(TPMSMS)) # UNDONE: make some arguments configurable
                    {
                        if (shift)
                            shift <- diff(gInfo[c(scrRow$group, g), "mzs"])
                        else
                            shift <- 0
                        
                        return(spectrumSimilarity(MSPeakLists, scrRow$group, g, MSLevel = 2, doPlot = FALSE, mzShift = shift))
                    }
                        
                    return(NA)
                }

                ret[, specSimilarity := sapply(group, getSim, shift = FALSE)]
                ret[, specSimilarityShift := sapply(group, getSim, shift = TRUE)]
            }
            else
                ret[, c("specSimilarity", "specSimilarityShift") := NA]
            
            return(ret)
        }), idcol = "precursor_group")
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
