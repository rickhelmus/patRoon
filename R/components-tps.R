#' @include main.R
#' @include components.R
#' @include utils-components.R
#' @include feature_groups-screening.R
#' @include feature_groups-screening-set.R
NULL

genTPSpecSimilarities <- function(obj, groupName1, groupName2, specSimParams, ...)
{
    # NOTE: groupName2 may have duplicate group names, which will be removed by intersect below. We don't need to repeat
    # their calculation. Just re-add them below.
    
    gn1 <- intersect(groupName1, groupNames(obj))
    gn2 <- intersect(groupName2, groupNames(obj))
    otherGN2 <- setdiff(groupName2, groupNames(obj))
    
    getSim <- function(shift)
    {
        sp <- specSimParams
        sp$shift <- shift
        
        if (length(gn1) == 0)
            return(NA_real_)
        
        ret <- numeric()
        if (length(gn2) > 0)
        {
            ret <- drop(spectrumSimilarity(obj, gn1, gn2, MSLevel = 2, drop = FALSE, specSimParams = sp, ...))
            names(ret) <- gn2
        }
        if (length(otherGN2) > 0)
            ret <- c(ret, setNames(rep(NA_real_, length(otherGN2)), otherGN2))
        
        ret <- ret[groupName2] # re-order and re-add duplicate group columns if needed
        
        return(ret)
    }
    
    return(list(specSimilarity = getSim("none"), specSimilarityPrec = getSim("precursor"),
                specSimilarityBoth = getSim("both")))
}

genTPAnnSimilarities <- function(parentFG, TPFGs, MSPeakLists, formulas, compounds)
{
    getAllAnnPLs <- function(grp)
    {
        annPLsForm <- annPLsComp <- NULL
        
        if (!is.null(formulas) && !is.null(formulas[[grp]]))
        {
            annPLsForm <- rbindlist(pruneList(lapply(seq_len(nrow(formulas[[grp]])), function(i)
            {
                annotatedPeakList(formulas, index = i, groupName = grp, MSPeakLists = MSPeakLists,
                                  onlyAnnotated = TRUE)
            })), fill = TRUE)
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
    
    annPLParent <- getAllAnnPLs(parentFG)
    
    if (is.null(annPLParent))
        return(list(fragmentMatches = rep(NA_integer_, length(TPFGs)),
                    neutralLossMatches = rep(NA_integer_, length(TPFGs))))
    
    parentFrags <- unique(annPLParent$ion_formula); parentNLs <- unique(annPLParent$neutral_loss)
    annPLTPs <- lapply(TPFGs, getAllAnnPLs)
    fragMatches <- sapply(annPLTPs, function(ann)
    {
        if (is.null(ann))
            return(NA_integer_)
        return(sum(unique(ann$ion_formula) %chin% parentFrags))
    })
    NLMatches <- sapply(annPLTPs, function(ann)
    {
        if (is.null(ann))
            return(NA_integer_)
        return(sum(unique(ann$neutral_loss) %chin% parentNLs))
    })
    
    return(list(fragmentMatches = fragMatches,
                neutralLossMatches = NLMatches))
}

doGenComponentsTPs <- function(fGroups, fGroupsTPs, ignoreParents, TPs, MSPeakLists, formulas, compounds, minRTDiff,
                               specSimParams)
{
    if (length(fGroups) == 0 || (!is.null(TPs) && length(TPs) == 0))
        return(componentsTPs(componentInfo = data.table(), components = list()))
    
    hash <- makeHash(fGroups, fGroupsTPs, ignoreParents, TPs, MSPeakLists, formulas, compounds, minRTDiff,
                     specSimParams)
    cd <- loadCacheData("componentsTPs", hash)
    if (!is.null(cd))
        return(cd)
    
    if (ignoreParents)
        fGroupsTPs <- delete(fGroupsTPs, j = names(fGroups))
    
    if (length(fGroupsTPs) == 0)
    {
        msg <- "fGroupsTPs is doesn't contain any feature groups"
        if (ignoreParents)
            msg <- paste(msg, "or doesn't contain any unique feature groups")
        warning(msg, call. = FALSE)
        
        return(componentsTPs(componentInfo = data.table(), components = list()))
    }
    
    gInfoParents <- groupInfo(fGroups); gInfoTPs <- groupInfo(fGroupsTPs)
    
    cat("Linking parents and TPs ...\n")
    
    prepareComponent <- function(cmp, parentFG)
    {
        # UNDONE: do more checks etc
        
        # dummy intensity value so e.g. plotSpectrum works            
        cmp[, intensity := 1]
        
        cmp[, c("ret", "mz") := gInfoTPs[group, c("rts", "mzs")]]
        cmp[, retDiff := ret - gInfoParents[parentFG, "rts"]]
        cmp[, mzDiff := mz - gInfoParents[parentFG, "mzs"]]
        
        cmp[, retDir := fcase((retDiff + minRTDiff) < 0, -1,
                             (retDiff - minRTDiff) > 0, 1,
                             default = 0)]
        
        if (!is.null(MSPeakLists))
        {
            if (nrow(cmp) > 0 && !is.null(MSPeakLists[[parentFG]][["MSMS"]]))
            {
                sims <- genTPSpecSimilarities(MSPeakLists, parentFG, cmp$group, specSimParams = specSimParams)
                cmp[, (names(sims)) := sims]
            }
            else
                cmp[, c("specSimilarity", "specSimilarityPrec", "specSimilarityBoth") := NA_real_]

            if (!is.null(formulas) || !is.null(compounds))
            {
                sims <- genTPAnnSimilarities(parentFG, cmp$group, MSPeakLists, formulas, compounds)
                cmp[, (names(sims)) := sims]
            }
        }
        
        return(cmp)
    }
    
    compTab <- NULL
    if (is.null(TPs))
    {
        # simply link each given parent with all given TPs, while relying on filter() to get sensible components
        
        parentCount <- length(fGroups)
        prog <- openProgBar(0, parentCount)
        
        compTab <- rbindlist(Map(names(fGroups), seq_len(parentCount), f = function(grp, i)
        {
            grpsTPs <- names(fGroupsTPs)
            
            comp <- data.table(group = grpsTPs, parent_group = grp)
            comp <- prepareComponent(comp, grp)
            
            # NOTE: name afterwards as the component may have been filtered
            comp[, TP_name := paste0(grp, "-TP", seq_len(nrow(comp)))]
            
            setTxtProgressBar(prog, i)
            
            return(comp)
        }), idcol = "parent_name")
        
        close(prog)
    }
    else
    {
        # for every parent:
        #   - check for any matching fGroups (based on mass)
        #   - if none, skip
        #   - similarly, check which parents are present (only mz)
        #   - for any fGroup that matches the parent:
        #       - filter TPs (retention, intensity, ...)
        
        parentFGMapping <- linkParentsToFGroups(TPs, fGroups)
        TPFGMapping <- linkTPsToFGroups(TPs, fGroupsTPs)
        pars <- parents(TPs)
        
        parentCount <- length(names(TPs))
        prog <- openProgBar(0, parentCount)
        
        compTab <- rbindlist(Map(names(TPs), products(TPs), seq_len(parentCount), f = function(pname, prods, i)
        {
            parentFGs <- parentFGMapping[name == pname][["group"]]
            TPs <- TPFGMapping[TP_name %in% prods$name]
            
            if (length(parentFGs) == 0 || nrow(TPs) == 0)
                comps <- NULL
            else
            {
                # limit columns a bit to not bloat components too much
                # UNDONE: column selection OK?
                prodCols <- c("name", "SMILES", "InChI", "InChIKey", "formula", "CID", "mass", "retDir", "trans_add",
                              "trans_sub", "deltaMZ")
                prods <- prods[, intersect(names(prods), prodCols), with = FALSE]
                
                comps <- rbindlist(sapply(parentFGs, function(parentFG)
                {
                    ret <- merge(TPs, prods, by.x = "TP_name", by.y = "name")
                    setnames(ret, "retDir", "TP_retDir")
                    
                    ret <- prepareComponent(ret, parentFG)
                    
                    if (!is.null(ret[["formula"]])) # eg TRUE for BT
                    {
                        parentForm <- pars[name == pname]$formula
                        ret[, formulaDiff := sapply(formula, subtractFormula, formula2 = parentForm)]
                    }
                    
                    return(ret)
                }, simplify = FALSE), idcol = "parent_group")
            }
            
            setTxtProgressBar(prog, i)
            
            return(comps)
        }), idcol = "parent_name")
        
        close(prog)
    }
    
    if (nrow(compTab) > 0)
    {
        compTab[, name := paste0("CMP", .GRP), by = c("parent_name", "parent_group")]
        compTab[, links := list(list(unique(name))), by = c("TP_name", "group")] # link to other components having this TP
        compTab[, links := mapply(links, name, FUN = setdiff)] # remove self-links
        
        compList <- split(compTab[, -"name"], by = c("parent_name", "parent_group"), keep.by = FALSE)
        
        compInfo <- unique(compTab[, c("name", "parent_name", "parent_group")])
        
        if (!is.null(TPs))
        {
            pars <- parents(TPs)
            cols <- c("formula", "SMILES", "InChI", "InChIKey", "CID", "neutralMass", "rt", "mz") # more/less?
            cols <- intersect(names(pars), cols)
            cols <- cols[sapply(cols, function(cl) any(!is.na(pars[[cl]])))]
            targetCols <- paste0("parent_", cols)
            compInfo[, (targetCols) := pars[match(parent_name, pars$name), cols, with = FALSE]]
        }
        
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

setMethod("collapseComponents", "componentsTPs", function(obj)
{
    obj@components <- lapply(obj@components, function(cmp)
    {
        keep <- intersect(c("TP_name", "group", "set"), names(cmp))
        cmp <- cmp[, keep, with = FALSE]
        cmp[, TP_name := paste0(unique(TP_name), collapse = ","), by = "group"]
        if (!is.null(cmp[["set"]]))
            cmp[, set := paste0(unique(unlist(strsplit(set, ","))), collapse = ","), by = "group"]
        return(unique(cmp, by = "group"))
    })
    return(obj)
})

#' @describeIn componentsTPs Returns all component data in a table.
#' @export
setMethod("as.data.table", "componentsTPs", function(x)
{
    # as regular method, but get rid of double links column when merging
    # components / componentInfo
    
    if (length(x) > 0)
        x@componentInfo <- x@componentInfo[, -"links"]
    
    return(callNextMethod(x))
})

#' @export
setMethod("filter", "componentsTPs", function(obj, ..., retDirMatch = FALSE,
                                              minSpecSim = NULL, minSpecSimPrec = NULL, minSpecSimBoth = NULL,
                                              minFragMatches = NULL, minNLMatches = NULL, formulas = NULL,
                                              verbose = TRUE, negate = FALSE)
{
    # UNDONE: optionally remove TPs with equal formula as parent (how likely are these?)
    # UNDONE: if formulas is set, also remove fGroups without assignment? Otherwise document!
    # UNDONE: also filter set separate similarities?
    
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertNumber, . ~ minSpecSim + minSpecSimPrec + minSpecSimBoth + minFragMatches + minNLMatches,
           lower = 0, finite = TRUE, null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertClass(formulas, "formulas", null.ok = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ retDirMatch + verbose + negate, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    if (length(obj) == 0)
        return(obj)
    
    if (!is.null(formulas) && is.null(obj[[1]][["trans_add"]]))
        stop("formula filter is only available for logic TP products")
    
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
    
    anyTPFilters <- retDirMatch || !is.null(minSpecSim) || !is.null(minSpecSimPrec) || !is.null(minSpecSimBoth) ||
        !is.null(minFragMatches) || !is.null(minNLMatches) || !is.null(formulas)
    
    if (anyTPFilters)
    {
        obj <- delete(obj, j = function(ct, cmp, ...)
        {
            ct <- copy(ct)
            ct[, keep := TRUE]
            
            if (retDirMatch)
            {
                if (negate)
                    ct[TP_retDir != 0 & retDir != 0, keep := TP_retDir != retDir]
                else
                    ct[TP_retDir != 0 & retDir != 0, keep := TP_retDir == retDir]
            }
            
            ct <- minColFilter(ct, "specSimilarity", minSpecSim)
            ct <- minColFilter(ct, "specSimilarityPrec", minSpecSimPrec)
            ct <- minColFilter(ct, "specSimilarityBoth", minSpecSimBoth)
            ct <- minColFilter(ct, "fragmentMatches", minFragMatches)
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
                
                parentFG <- componentInfo(obj)[name == cmp]$parent_group
                if (!is.null(formulas[[parentFG]]))
                {
                    # filter results where subtraction of any of the parent formulas is impossible
                    ct[keep == TRUE & nzchar(trans_sub), keep := sapply(trans_sub, canSub, formulas[[parentFG]])]
                }
                
                # filter results where addition is not part of TP candidate formulas
                ct[keep == TRUE & nzchar(trans_add), keep := mapply(trans_add, annotations(formulas)[group],
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
                      names(obj), cInfo$parent_name, cInfo$parent_group,
                      sapply(cTable, function(cmp) paste0(cmp$TP_name, collapse = ", ")),
                      sapply(cTable, function(cmp) paste0(unique(cmp$group), collapse = ", ")))
    makeGraph(obj, onlyLinked, titles)
})

#' @export
setMethod("generateComponentsTPs", "featureGroups", function(fGroups, fGroupsTPs = fGroups, ignoreParents = FALSE,
                                                             TPs = NULL, MSPeakLists = NULL, formulas = NULL,
                                                             compounds = NULL, minRTDiff = 20,
                                                             specSimParams = getDefSpecSimParams())
{
    # UNDONE: doc when fGroups/fGroupsTPs needs to be fGroupsScreening
    
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertClass, . ~ fGroupsTPs + TPs + MSPeakLists + formulas + compounds,
           c("featureGroups", "transformationProducts", "MSPeakLists", "formulas", "compounds"),
           null.ok = c(FALSE, TRUE, TRUE, TRUE, TRUE), fixed = list(add = ac))
    checkmate::assertFlag(ignoreParents, add = ac)
    checkmate::assertNumber(minRTDiff, lower = 0, finite = TRUE, add = ac)
    assertSpecSimParams(specSimParams, add = ac)
    checkmate::reportAssertions(ac)
    
    if (!is.null(TPs) && needsScreening(TPs) &&
        (!inherits(fGroups, "featureGroupsScreening") || !inherits(fGroupsTPs, "featureGroupsScreening")))
        stop("Input feature groups need to be screened for parents/TPs!")

    return(doGenComponentsTPs(fGroups, fGroupsTPs, ignoreParents, TPs, MSPeakLists, formulas, compounds, minRTDiff,
                              specSimParams = specSimParams))
})

#' @export
setMethod("generateComponentsTPs", "featureGroupsSet", function(fGroups, fGroupsTPs = fGroups, ignoreParents = FALSE,
                                                                TPs = NULL, MSPeakLists = NULL, formulas = NULL,
                                                                compounds = NULL, minRTDiff = 20,
                                                                specSimParams = getDefSpecSimParams())
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertClass, . ~ fGroupsTPs + TPs + MSPeakLists + formulas + compounds,
           c("featureGroupsSet", "transformationProducts", "MSPeakListsSet", "formulasSet", "compoundsSet"),
           null.ok = c(FALSE, TRUE, TRUE, TRUE, TRUE), fixed = list(add = ac))
    checkmate::assertFlag(ignoreParents, add = ac)
    checkmate::assertNumber(minRTDiff, lower = 0, finite = TRUE, add = ac)
    assertSpecSimParams(specSimParams, add = ac)
    checkmate::reportAssertions(ac)

    if (!is.null(TPs) && needsScreening(TPs) &&
        (!inherits(fGroups, "featureGroupsScreeningSet") || !inherits(fGroupsTPs, "featureGroupsScreeningSet")))
        stop("Input feature groups need to be screened for parents/TPs!")
    
    ret <- doGenComponentsTPs(fGroups, fGroupsTPs, ignoreParents, TPs, MSPeakLists, formulas, compounds,
                              minRTDiff, specSimParams = specSimParams)
    
    # UNDONE: more efficient method to get set specific fGroups?
    gNamesTPsSets <- sapply(sets(fGroupsTPs), function(s) names(fGroupsTPs[, sets = s]), simplify = FALSE)
    
    unsetMSPeakLists <- checkAndUnSetOther(sets(fGroupsTPs), MSPeakLists, "MSPeakLists", TRUE)
    unsetFormulas <- checkAndUnSetOther(sets(fGroupsTPs), formulas, "formulas", TRUE)
    unsetCompounds <- checkAndUnSetOther(sets(fGroupsTPs), compounds, "compounds", TRUE)

    cat("Adding sets related data...\n")
    if (length(ret) > 0)
    {
        ret@components <- withProg(length(ret), FALSE, Map(ret@components, ret@componentInfo$parent_group, f = function(cmp, parentFG)
        {
            # mark set presence
            cmp[, set := sapply(group, function(g)
            {
                paste0(names(which(sapply(gNamesTPsSets, function(n) g %chin% n))), collapse = ",")
            })]
            
            for (s in sets(fGroupsTPs))
            {
                if (!is.null(unsetMSPeakLists[[s]]))
                {
                    # calculate per set spectrum similarities
                    simColNames <- paste0(c("specSimilarity", "specSimilarityPrec", "specSimilarityBoth"), "-", s)
                    grpsInSet <- intersect(cmp$group, groupNames(unsetMSPeakLists[[s]]))
                    
                    if (length(grpsInSet) > 0 && !is.null(unsetMSPeakLists[[s]][[parentFG]][["MSMS"]]))
                    {
                        sims <- genTPSpecSimilarities(unsetMSPeakLists[[s]], parentFG, grpsInSet,
                                                      specSimParams = specSimParams)
                        cmp[match(grpsInSet, group), (simColNames) := sims]
                    }
                    else
                        cmp[, (simColNames) := NA_real_]
                    
                    if (!is.null(unsetFormulas[[s]]) || !is.null(unsetCompounds[[s]]))
                    {
                        # calculate per set spectrum similarities
                        sims <- genTPAnnSimilarities(parentFG, cmp$group, unsetMSPeakLists[[s]], unsetFormulas[[s]],
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
    }
    
    return(ret)
})
