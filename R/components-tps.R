#' @include main.R
#' @include components.R
#' @include utils-components.R
#' @include feature_groups-screening.R
#' @include feature_groups-screening-set.R
NULL

genTPSpecSimilarities <- function(obj, groupName1, groupName2, ...)
{
    gn1 <- intersect(groupName1, groupNames(obj))
    gn2 <- intersect(groupName2, groupNames(obj))
    otherGN2 <- setdiff(groupName2, groupNames(obj))
    
    getSim <- function(shift)
    {
        if (length(gn1) == 0)
            return(NA_real_)
        
        ret <- numeric()
        if (length(gn2) > 0)
        {
            ret <- drop(do.call(spectrumSimilarity, list(obj, gn1, gn2, shift = shift, MSLevel = 2, drop = FALSE, ...)))
            names(ret) <- gn2
        }
        if (length(otherGN2) > 0)
        {
            ret <- c(ret, setNames(rep(NA_real_, length(otherGN2)), otherGN2))
            ret <- ret[groupName2]
        }
        
        return(ret)
    }
    
    return(list(specSimilarity = getSim("none"), specSimilarityPrec = getSim("precursor"),
                specSimilarityBoth = getSim("both")))
}

genTPAnnSimilarities <- function(precFG, TPFGs, MSPeakLists, formulas, compounds)
{
    getAllAnnPLs <- function(grp)
    {
        annPLsForm <- annPLsComp <- NULL
        
        if (!is.null(formulas) && !is.null(formulas[[grp]]))
        {
            precs <- unique(formulas[[grp]]$neutral_formula)
            annPLsForm <- rbindlist(pruneList(lapply(precs, function(p)
            {
                ret <- annotatedPeakList(formulas, precursor = p, groupName = grp, MSPeakLists = MSPeakLists,
                                         onlyAnnotated = TRUE)
            })))
        }
        
        if (!is.null(compounds) && !is.null(compounds[[grp]]))
        {
            annPLsComp <- rbindlist(pruneList(lapply(seq_len(nrow(compounds[[grp]])), function(i)
            {
                # NOTE: we don't set formulas as we handle them separately in order to get all their results as well
                annotatedPeakList(compounds, index = i, groupName = grp, MSPeakLists = MSPeakLists, formulas = NULL,
                                  onlyAnnotated = TRUE)
            })), fill = TRUE)
        }
        
        # NOTE: use rbindlist to force DT method, which has the fill argument we need
        # NOTE: rbind/rbindlist deals with NULLs
        return(rbindlist(list(annPLsForm, annPLsComp), fill = TRUE)) 
    }
    
    annPLPrec <- getAllAnnPLs(precFG)
    
    if (is.null(annPLPrec))
        return(list(fragmentMatches = rep(NA_integer_, length(TPFGs)),
                    neutralLossMatches = rep(NA_integer_, length(TPFGs))))
    
    precFrags <- unique(annPLPrec$formula); precNLs <- unique(annPLPrec$neutral_loss)
    annPLTPs <- lapply(TPFGs, getAllAnnPLs)
    fragMatches <- lapply(annPLTPs, function(ann)
    {
        if (is.null(ann))
            return(NA_integer_)
        return(sum(unique(ann$formula) %chin% precFrags))
    })
    NLMatches <- lapply(annPLTPs, function(ann)
    {
        if (is.null(ann))
            return(NA_integer_)
        return(sum(unique(ann$neutral_loss) %chin% precNLs))
    })
    
    return(list(fragmentMatches = fragMatches,
                neutralLossMatches = NLMatches))
}

doGenComponentsTPs <- function(fGroups, fGroupsTPs, ignorePrecursors, pred, MSPeakLists, formulas, compounds, minRTDiff,
                               specSimArgs)
{
    if (length(fGroups) == 0)
        return(componentsTPs(componentInfo = data.table(), components = list()))
    
    hash <- makeHash(fGroups, fGroupsTPs, pred, MSPeakLists, formulas, compounds, minRTDiff, specSimArgs)
    cd <- loadCacheData("componentsTPs", hash)
    if (!is.null(cd))
        return(cd)
    
    if (ignorePrecursors)
        fGroupsTPs <- delete(fGroupsTPs, j = names(fGroups))
    
    if (length(fGroupsTPs) == 0)
    {
        msg <- "fGroupsTPs is doesn't contain any feature groups"
        if (ignorePrecursors)
            msg <- paste(msg, "or doesn't contain any (unique) feature groups")
        warning(msg, call. = FALSE)
        
        return(componentsTPs(componentInfo = data.table(), components = list()))
    }
    
    gInfoPrec <- groupInfo(fGroups); gInfoTPs <- groupInfo(fGroupsTPs)
    
    cat("Linking precursors and TPs ...\n")
    
    prepareComponent <- function(cmp, precFG)
    {
        # UNDONE: do more checks etc
        
        # dummy intensity value so e.g. plotSpec works            
        cmp[, intensity := 1]
        
        cmp[, c("ret", "mz") := gInfoTPs[group, c("rts", "mzs")]]
        cmp[, retDiff := gInfoPrec[precFG, "rts"] - ret]
        cmp[, mzDiff := gInfoPrec[precFG, "mzs"] - mz]
        
        if (minRTDiff > 0)
        {
            rtDiffs <- gInfoTPs[cmp$group, "rts"] - gInfoPrec[precFG, "rts"]
            cmp <- cmp[RTDir == 0 | abs(rtDiffs) <= minRTDiff | (rtDiffs < 0 & RTDir < 0) | (rtDiffs > 0 & RTDir > 0)]
        }
        
        if (!is.null(MSPeakLists))
        {
            if (nrow(cmp) > 0 && !is.null(MSPeakLists[[precFG]][["MSMS"]]))
            {
                sims <- do.call(genTPSpecSimilarities, c(list(MSPeakLists, precFG, cmp$group), specSimArgs))
                cmp[, (names(sims)) := sims]
            }
            else
                cmp[, c("specSimilarity", "specSimilarityPrec", "specSimilarityBoth") := NA_real_]

            if (!is.null(formulas) || !is.null(compounds))
            {
                sims <- genTPAnnSimilarities(precFG, cmp$group, MSPeakLists, formulas, compounds)
                cmp[, (names(sims)) := sims]
            }
        }
        
        return(cmp)
    }
    
    compTab <- NULL
    if (is.null(pred))
    {
        # simply link each given parent with all given TPs, while relying on the filtering of prepareComponent() to get
        # sensible components
        
        precCount <- length(fGroups)
        prog <- openProgBar(0, precCount)
        
        compTab <- rbindlist(Map(names(fGroups), seq_len(precCount), f = function(grp, i)
        {
            grpsTPs <- names(fGroupsTPs)
            
            comp <- data.table(group = grpsTPs, RTDir = 0, precursor_group = grp)
            comp <- prepareComponent(comp, grp)
            
            # NOTE: name afterwards as the component may have been filtered
            comp[, TP_name := paste0(grp, "-TP", seq_len(nrow(comp)))]
            
            setTxtProgressBar(prog, i)
            
            return(comp)
        }), idcol = "precursor_name")
        
        close(prog)
    }
    else
    {
        # for every precursor:
        #   - check for any matching fGroups (based on mass)
        #   - if none, skip
        #   - similarly, check which precursors are present (only mz)
        #   - for any fGroup that matches the precursor:
        #       - filter TPs (retention, intensity, ...)
        
        precFGMapping <- linkPrecursorsToFGroups(pred, fGroups)
        TPFGMapping <- linkTPsToFGroups(pred, fGroupsTPs)
        susps <- suspects(pred)
        
        precCount <- length(names(pred))
        prog <- openProgBar(0, precCount)
        
        compTab <- rbindlist(Map(names(pred), predictions(pred), seq_len(precCount), f = function(pname, preds, i)
        {
            precFGs <- precFGMapping[name == pname][["group"]]
            TPs <- TPFGMapping[TP_name %in% preds$name]
            
            if (length(precFGs) == 0 || nrow(TPs) == 0)
                comps <- NULL
            else
            {
                # limit columns a bit to not bloat components too much
                # UNDONE: column selection OK?
                predCols <- c("name", "InChIKey", "formula", "mass", "RTDir",
                              "reaction_add", "reaction_sub", "deltaMZ")
                preds <- preds[, intersect(names(preds), predCols), with = FALSE]
                
                comps <- rbindlist(sapply(precFGs, function(precFG)
                {
                    ret <- merge(TPs, preds, by.x = "TP_name", by.y = "name")
                    
                    # merge rows with duplicate fGroups, for instance, caused by different TPs with equal mass
                    ret[, TP_name := paste0(TP_name, collapse = ","), by = "group"]
                    ret <- unique(ret, by = "group")
                    
                    ret <- prepareComponent(ret, precFG)
                    
                    if (!is.null(ret[["formula"]])) # eg TRUE for BT
                    {
                        precForm <- susps[name == pname]$formula
                        ret[, formulaDiff := sapply(formula, subtractFormula, formula1 = precForm)]
                    }
                    
                    return(ret)
                }, simplify = FALSE), idcol = "precursor_group")
            }
            
            setTxtProgressBar(prog, i)
            
            return(comps)
        }), idcol = "precursor_name")
        
        close(prog)
    }
    
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
    
    return(ret[])
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
setMethod("filter", "componentsTPs", function(obj, ..., minSpecSim = NULL, minSpecSimPrec = NULL, minSpecSimBoth = NULL,
                                              minFragMatches = NULL, minNLMatches = NULL, formulas = NULL,
                                              verbose = TRUE, negate = FALSE)
{
    # UNDONE: if formulas is set, also remove fGroups without assignment? Otherwise document!
    # UNDONE: also filter set separate similarities?
    
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertNumber, . ~ minSpecSim + minSpecSimPrec + minSpecSimBoth + minFragMatches + minNLMatches,
           lower = 0, finite = TRUE, null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertClass(formulas, "formulas", null.ok = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ verbose + negate, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    if (length(obj) == 0)
        return(obj)
    
    if (!is.null(formulas) && is.null(obj[[1]][["reaction_add"]]))
        stop("formula filter is only available for logic TP predictions")
    
    old <- obj
    
    minColFilter <- function(ct, col, val)
    {
        if (!is.null(val) && val > 0)
        {
            if (!is.null(ct[[col]]))
            {
                if (negate)
                    pred <- function(x, v) is.na(x) | x < v
                else
                    pred <- function(x, v) !is.na(x) & numGTE(x, v)
                ct[keep == TRUE, keep := pred(get(col), val)]
            }
        }
        return(ct)
    }
    
    anyTPFilters <- !is.null(minSpecSim) || !is.null(minSpecSimPrec) || !is.null(minSpecSimBoth) ||
        !is.null(minFragMatches) || !is.null(minNLMatches) || !is.null(formulas)
    
    if (anyTPFilters)
    {
        obj <- delete(obj, j = function(ct, cmp, ...)
        {
            ct <- copy(ct)
            ct[, keep := TRUE]
            
            ct <- minColFilter(ct, "specSimilarity", minSpecSim)
            ct <- minColFilter(ct, "specSimilarityPrec", minSpecSimPrec)
            ct <- minColFilter(ct, "minSpecSimBoth", minSpecSimBoth)
            ct <- minColFilter(ct, "fragMatches", minFragMatches)
            ct <- minColFilter(ct, "neutralLossMatches", minNLMatches)
            
            if (!is.null(formulas))
            {
                # check if subtracting is possible, ie by checking if subtraction doesn't lead to negative element
                # counts
                canSub <- function(f, ft)
                {
                    if (is.null(ft) || length(f) == 0 || !nzchar(f))
                        return(TRUE) # UNDONE?
                    candidateForms <- unique(ft$neutral_formula)
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
                
                precFG <- componentInfo(obj)[name == cmp]$precursor_group
                if (!is.null(formulas[[precFG]]))
                {
                    # filter results where subtraction of any of the precursor formulas is impossible
                    ct[keep == TRUE & nzchar(reaction_sub), keep := sapply(reaction_sub, canSub, formulas[[precFG]])]
                }
                
                # filter results where addition is not part of TP candidate formulas
                ct[keep == TRUE & nzchar(reaction_add), keep := mapply(reaction_add, formulaTable(formulas)[group],
                                                                       FUN = canSub)]
            }
            
            return(!ct$keep)
        })
        
        if (verbose)
            printComponentsFiltered(old, obj)
    }
    
    
    if (...length() > 0)
        return(callNextMethod(obj, ..., verbose = verbose, negate = negate))
    
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
setMethod("generateComponentsTPs", "featureGroups", function(fGroups, fGroupsTPs = fGroups, ignorePrecursors = FALSE,
                                                             pred = NULL, MSPeakLists = NULL, formulas = NULL,
                                                             compounds = NULL, minRTDiff = 20,
                                                             simMethod, removePrecursor = FALSE, mzWeight = 0,
                                                             intWeight = 1, absMzDev = 0.005,
                                                             relMinIntensity = 0.05, minSimMSMSPeaks = 0)
{
    # UNDONE: optionally remove TPs with equal formula as parent (how likely are these?)
    # UNDONE: doc when fGroups/fGroupsTPs needs to be fGroupsScreening
    
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertClass, . ~ fGroupsTPs + pred + MSPeakLists + formulas + compounds,
           c("featureGroups", "TPPredictions", "MSPeakLists", "formulas", "compounds"),
           null.ok = c(FALSE, TRUE, TRUE, TRUE, TRUE), fixed = list(add = ac))
    aapply(checkmate::assertFlag, . ~ ignorePrecursors + removePrecursor, fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ minRTDiff + mzWeight + intWeight + absMzDev + relMinIntensity,
           lower = 0, finite = TRUE, fixed = list(add = ac))
    checkmate::assertChoice(simMethod, c("cosine", "jaccard"), add = ac)
    checkmate::assertCount(minSimMSMSPeaks, add = ac)
    checkmate::reportAssertions(ac)
    
    if (!is.null(pred) && needsScreening(pred) &&
        (!inherits(fGroups, "featureGroupsScreening") || !inherits(fGroupsTPs, "featureGroupsScreening")))
        stop("Input feature groups need to be screened for (TP) suspects!")

    return(doGenComponentsTPs(fGroups, fGroupsTPs, ignorePrecursors, pred, MSPeakLists, formulas, compounds, minRTDiff,
                              specSimArgs = list(method = simMethod, removePrecursor = removePrecursor,
                              mzWeight = mzWeight, intWeight = intWeight, absMzDev = absMzDev,
                              relMinIntensity = relMinIntensity, minPeaks = minSimMSMSPeaks)))
})

#' @export
setMethod("generateComponentsTPs", "featureGroupsSet", function(fGroups, fGroupsTPs = fGroups, ignorePrecursors = FALSE,
                                                                pred = NULL, MSPeakLists = NULL, formulas = NULL,
                                                                compounds = NULL, minRTDiff = 20,
                                                                simMethod, removePrecursor = FALSE, mzWeight = 0,
                                                                intWeight = 1, absMzDev = 0.005,
                                                                relMinIntensity = 0.05, minSimMSMSPeaks = 0)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertClass, . ~ fGroupsTPs + pred + MSPeakLists + formulas + compounds,
           c("featureGroupsSet", "TPPredictions", "MSPeakListsSet", "formulasSet", "compoundsSet"),
           null.ok = c(FALSE, TRUE, TRUE, TRUE, TRUE), fixed = list(add = ac))
    aapply(checkmate::assertFlag, . ~ ignorePrecursors + removePrecursor, fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ minRTDiff + mzWeight + intWeight + absMzDev + relMinIntensity, lower = 0,
           finite = TRUE, fixed = list(add = ac))
    checkmate::assertChoice(simMethod, c("cosine", "jaccard"), add = ac)
    aapply(checkmate::assertFlag, . ~ ignorePrecursors + removePrecursor, fixed = list(add = ac))
    checkmate::assertCount(minSimMSMSPeaks, add = ac)
    checkmate::reportAssertions(ac)

    if (!is.null(pred) && needsScreening(pred) &&
        (!inherits(fGroups, "featureGroupsScreeningSet") || !inherits(fGroupsTPs, "featureGroupsScreeningSet")))
        stop("Input feature groups need to be screened for (TP) suspects!")
    
    specSimArgs <- list(method = simMethod, removePrecursor = removePrecursor, mzWeight = mzWeight,
                        intWeight = intWeight, absMzDev = absMzDev, relMinIntensity = relMinIntensity,
                        minPeaks = minSimMSMSPeaks)
    
    ret <- doGenComponentsTPs(fGroups, fGroupsTPs, ignorePrecursors, pred, MSPeakLists, formulas, compounds,
                              minRTDiff, specSimArgs = specSimArgs)
    
    # UNDONE: more efficient method to get set specific fGroups?
    gNamesTPsSets <- sapply(sets(fGroupsTPs), function(s) names(fGroupsTPs[, sets = s]), simplify = FALSE)
    
    unsetORNULL <- function(x) if (!is.null(x)) sapply(sets(x), unset, obj = x, simplify = FALSE) else NULL
    
    unsetMSPeakLists <- unsetORNULL(MSPeakLists)
    unsetFormulas <- unsetORNULL(formulas)
    unsetCompounds <- unsetORNULL(compounds)

    cat("Adding sets related data...\n")
    ret@components <- withProg(length(ret), FALSE,
                               Map(ret@components, ret@componentInfo$precursor_group, f = function(cmp, precFG)
    {
        for (s in sets(fGroupsTPs))
        {
            # mark set presence
            set(cmp, j = s, value = cmp$group %in% gNamesTPsSets[[s]])
            
            if (!is.null(unsetMSPeakLists))
            {
                # calculate per set spectrum similarities
                simColNames <- paste0(c("specSimilarity", "specSimilarityPrec", "specSimilarityBoth"), "-", s)
                grpsInSet <- intersect(cmp$group, groupNames(unsetMSPeakLists[[s]]))
                
                if (length(grpsInSet) > 0 && !is.null(unsetMSPeakLists[[s]][[precFG]][["MSMS"]]))
                {
                    sims <- do.call(genTPSpecSimilarities, c(list(unsetMSPeakLists[[s]], precFG, grpsInSet), specSimArgs))
                    cmp[match(grpsInSet, group), (simColNames) := sims]
                }
                else
                    cmp[, (simColNames) := NA_real_]
                
                if (!is.null(unsetFormulas) || !is.null(unsetCompounds))
                {
                    # calculate per set spectrum similarities
                    sims <- genTPAnnSimilarities(precFG, cmp$group, unsetMSPeakLists[[s]], unsetFormulas[[s]],
                                                 unsetCompounds[[s]])
                    annColNames <- paste0(c("fragmentMatches", "neutralLossMatches"), "-", s)
                    cmp[, (annColNames) := sims]
                }
            }
            
            # move spec similarity columns to end
            setcolorder(cmp, setdiff(names(cmp), grep("^specSimilarity", names(cmp), value = TRUE)))
            
            # ... and annotation sim columns after that
            setcolorder(cmp, setdiff(names(cmp), grep("Matches", names(cmp), fixed = TRUE, value = TRUE)))
            
            doProgress()
        }
        
        return(cmp)
    }))
        
    return(ret)
})
