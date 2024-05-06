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
        return(list(totalFragmentMatches = rep(NA_integer_, length(TPFGs)),
                    totalNeutralLossMatches = rep(NA_integer_, length(TPFGs))))
    
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
    
    return(list(totalFragmentMatches = fragMatches,
                totalNeutralLossMatches = NLMatches))
}

getTPComponCandidatesScr <- function(prods, TPGroups, MSPeakLists, formulas, compounds, parFG, parFormula, parIK)
{
    if (nrow(prods) == 0)
        return(NULL)
    
    # limit columns a bit to not bloat components too much
    # UNDONE: column selection OK?
    tab <- copy(prods)
    tab <- subsetDTColumnsIfPresent(tab, c("name", "compoundName", "SMILES", "InChI", "InChIKey", "formula",
                                           "molNeutralized", "CID", "mass", "retDir", "trans_add",
                                           "trans_sub", "deltaMZ", "similarity", "fitFormula", "fitCompound", "simSusp",
                                            "simSuspSMILES", "TPScore", "mergedBy", "coverage"))
    tab <- unique(tab, by = "name") # omit duplicates from different routes
    setnames(tab, "retDir", "TP_retDir", skip_absent = TRUE)
    if (!is.null(tab[["formula"]]) && !is.null(parFormula))
        tab[, formulaDiff := sapply(formula, getFormulaDiffText, form2 = parFormula)]

    getAnnPL <- function(featAnn, grp, UID, ...)
    {
        if (is.null(UID) || is.null(featAnn) || is.null(MSPeakLists) || !grp %chin% groupNames(featAnn))
            return(NULL)
        wh <- which(UID == featAnn[[grp]]$UID)
        if (length(wh) == 0)
            return(NULL)
        return(annotatedPeakList(featAnn, groupName = grp, index = wh, MSPeakLists = MSPeakLists, onlyAnnotated = TRUE,
                                 ...))
    }
    getFragNLs <- function(grp, InChIKey, formula)
    {
        fAPL <- getAnnPL(formulas, grp, formula)
        cAPL <- getAnnPL(compounds, grp, getIKBlock1(InChIKey))
        APL <- rbindlist(list(fAPL, cAPL), fill = TRUE)
        if (nrow(APL) == 0)
            return(NULL)
        return(unique(APL[, c("ion_formula", "neutral_loss"), with = FALSE]))
    }

    hasForm <- !is.null(tab[["formula"]]); hasIK <- !is.null(tab[["InChIKey"]])
    if ((hasForm && !is.null(formulas)) || (hasIK && !is.null(compounds)))
    {
        parFragNLs <- getFragNLs(parFG, parIK, parFormula)    
        
        if (!is.null(parFragNLs) && nrow(parFragNLs) > 0)
        {
            tab[, c("fragmentMatches", "neutralLossMatches") := {
                fragNLs <- getFragNLs(TPGroups[.I], if (hasIK) NAToNULL(InChIKey) else NULL,
                                      if (hasForm) NAToNULL(formula) else NULL)
                if (!is.null(fragNLs) && nrow(fragNLs) > 0)
                    list(sum(fragNLs$ion_formula %chin% parFragNLs$ion_formula),
                         sum(fragNLs$neutral_loss %chin% parFragNLs$neutral_loss))
                else
                    list(NA_integer_, NA_integer_)
            }, by = seq_len(nrow(tab))]
        }
    }        
    return(tab)
}

mergeTPComponCandidatesTab <- function(compTab)
{
    if (is.null(compTab[["candidates"]]))
        return(compTab)
        
    compTab[, candidates := Map(name, group, candidates, f = function(cmpName, grp, tab)
    {
        tab <- copy(tab)
        tab[, c("cmpName", "group") := .(cmpName, grp)]
        setnames(tab, "name", "candidate_name")
        return(tab)
    })]
    
    allc <- rbindlist(compTab$candidates, fill = TRUE)
    compTab <- merge(compTab, allc, by.x = c("name", "group"), by.y = c("cmpName", "group"), all.x = TRUE, sort = FALSE)
    
    return(compTab)
}

doGenComponentsTPs <- function(fGroups, fGroupsTPs, TPs, MSPeakLists, formulas, compounds, ignoreParents, minRTDiff,
                               specSimParams, parallel)
{
    fromTPs <- !is.null(TPs)
    parFromScr <- fromTPs && parentsFromScreening(TPs)
    TPsFromScr <- fromTPs && TPsFromScreening(TPs)
    
    emptyRet <- componentsTPs(componentInfo = data.table(), components = list(), fromTPs = fromTPs,
                              parentsFromScreening = parFromScr, TPsFromScreening = TPsFromScr)
    
    if (length(fGroups) == 0 || (!is.null(TPs) && length(TPs) == 0))
        return(emptyRet)

    hash <- makeHash(fGroups, fGroupsTPs, TPs, MSPeakLists, formulas, compounds, ignoreParents, minRTDiff, specSimParams)
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
        
        return(emptyRet)
    }
    
    gInfoParents <- groupInfo(fGroups); gInfoTPs <- groupInfo(fGroupsTPs)
    haveFormulas <- !is.null(formulas) && length(formulas) > 0
    haveCompounds <- !is.null(compounds) && length(compounds) > 0

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
    
    # figure out components
    #   - for scr+unknowns: each component is a unique parent_name/parent_group pair --> parentFGMapping
    #   - for unknowns: each component is a unique susp name/group pair (ie scr table with these 2 cols)
    #   - for full unknowns: each component is a group
    # for each component
    #   - for scr: get products from parent_name, get component fGroups from intersection TPFGMapping/products, make candidates from TP data
    #   - for unknowns with TPs: consider all fGroups with annotation results, calculate their TP scores etc, keep those within thresholds and use calculated values for candidate tables
    #   - for unknowns with scr: same as above
    #   - for full unknowns: simply get component from prepareComponent(), no candidate table

    compInfo <- TPFGMapping <- NULL
    if (fromTPs)
    {
        compInfo <- linkParentsToFGroups(TPs, fGroups)
        setnames(compInfo, paste0("parent_", names(compInfo)))
        pars <- parents(TPs)
        cols <- c("formula", "SMILES", "InChI", "InChIKey", "CID", "neutralMass", "molNeutralized")
        cols <- intersect(names(pars), cols)
        cols <- cols[sapply(cols, function(cl) any(!is.na(pars[[cl]])))]
        targetCols <- paste0("parent_", cols)
        compInfo[, (targetCols) := pars[match(parent_name, pars$name), cols, with = FALSE]]
        compInfo[, c("parent_rt", "parent_mz") := gInfoParents[parent_group, c("rts", "mzs")]]
        setcolorder(compInfo, c("parent_name", "parent_group", "parent_rt", "parent_mz"))
        TPFGMapping <- linkTPsToFGroups(TPs, fGroupsTPs)
    }
    else
        compInfo <- data.table(parent_group = names(fGroups), parent_name = names(fGroups))

    compList <- doApply("Map", FALSE, compInfo$parent_name, compInfo$parent_group, f = function(parn, parfg)
    {
        cmpTab <- NULL
        if (fromTPs)
        {
            parRow <- parents(TPs)[name == parn][1]
            parForm <- NAToNULL(parRow[["formula"]]); parIK <- NAToNULL(parRow[["InChIKey"]])
            allFGCandidates <- unique(TPFGMapping$group)
            allFGCandidates <- setdiff(allFGCandidates, parfg)
            allFGCandidates <- intersect(allFGCandidates, names(fGroupsTPs))
            candidates <- pruneList(sapply(allFGCandidates, function(candFG)
            {
                getTPComponCandidatesScr(products(TPs)[[parn]][name %chin% TPFGMapping[group == candFG]$TP_name],
                                         candFG, MSPeakLists, formulas, compounds, parfg, parForm, parIK)
            }, simplify = FALSE), checkZeroRows = TRUE)
            cmpTab <- data.table(group = names(candidates), candidates = candidates)
        }
        else
        {
            cmpTab <- data.table(group = names(fGroupsTPs))
            cmpTab[, TP_name := paste0(parfg, "-TP", seq_len(nrow(cmpTab)))]
        }
        cmpTab <- prepareComponent(cmpTab, parfg)
        doProgress()
        return(cmpTab)
    })
    
    if (nrow(compInfo) > 0)
    {
        compInfo[, size := sapply(compList, nrow)]
        whNotEmpty <- compInfo[size != 0, which = TRUE]
        compList <- compList[whNotEmpty]
        compInfo <- compInfo[whNotEmpty]
        
        compInfo[, name := paste0("CMP", .I)]
        setcolorder(compInfo, "name")
        names(compList) <- compInfo$name

        if (length(compList) == 0)
        {
            compList <- Map(names(compList), compList, f = function(cmpName, cmpTab)
            {
                cmpTab <- copy(cmpTab)
                cmpTab[, links := list(list(character()))]
                return(cmpTab)
            })
        }
        else if (fromTPs)
        {
            # generate candidate links
            
            compListTab <- rbindlist(compList, idcol = "name", fill = TRUE)
            compListTab <- mergeTPComponCandidatesTab(compListTab)
            # UNDONE: change candidate_name to something more unique when we support multiple TPs objects
            compListTab[, links := list(list(unique(name))), by = c("candidate_name", "group")]
            compListTab[, links := Map(links, name, f = setdiff)] # remove self-links
        
            compList <- Map(names(compList), compList, f = function(cmpName, cmpTab)
            {
                cmpTab <- copy(cmpTab)
                cmpTab[, candidates := Map(group, candidates, f = function(grp, ct)
                {
                    ct <- copy(ct)
                    # UNDONE: as above
                    ct[compListTab[name == cmpName & group == grp], links := i.links, on = c(name = "candidate_name")]
                    return(ct)
                })]
                return(cmpTab)
            })
        }
        else
        {
            # generate fGroup links
            compListTab <- rbindlist(compList, idcol = "name", fill = TRUE)
            compListTab[, links := list(list(unique(name))), by = "group"]
            compListTab[, links := Map(links, name, f = setdiff)] # remove self-links

            compList <- Map(names(compList), compList, f = function(cmpName, cmpTab)
            {
                cmpTab <- copy(cmpTab)
                cmpTab[compListTab[name == cmpName], links := i.links, on = "group"]
                return(cmpTab)
            })
        }
        
        # overall links
        compInfo[, links := lapply(name, function(cn) unique(unlist(compListTab[name == cn]$links)))]
    }
    
    # UNDONE: update counting for candidates
    printf("Linked %d parents with %d TPs.\n", nrow(compInfo),
           if (length(compList) > 0) sum(sapply(compList, nrow)) else 0)
    
    ret <- componentsTPs(componentInfo = compInfo[], components = compList, fromTPs = fromTPs,
                         parentsFromScreening = parFromScr, TPsFromScreening = TPsFromScr)
    saveCacheData("componentsTPs", ret, hash)
    
    return(ret)
}

#' Components based on parent and transformation product (TP) linkage.
#'
#' This class is derived from \code{\link{components}} and is used to store components that result from linking feature
#' groups that are (predicted to be) parents with feature groups that (are predicted to be) transformation products. For
#' more details, see \code{\link{generateComponentsTPs}}.
#'
#' @param x,obj A \code{componentsTPs} object.
#'
#' @slot fromTPs A \code{logical} that is \code{TRUE} when the componentization was performed with
#'   \code{\link{transformationProducts}} data.
#'
#' @seealso \code{\link{components}} for other relevant methods and \code{\link{generateComponents}}
#'
#' @templateVar class componentsTPs
#' @template class-hierarchy
#'
#' @template components_noint
#' @export
componentsTPs <- setClass("componentsTPs", contains = "components", slots = c(fromTPs = "logical",
                                                                              parentsFromScreening = "logical",
                                                                              TPsFromScreening = "logical"))

setMethod("initialize", "componentsTPs",
          function(.Object, ...) callNextMethod(.Object, ..., algorithm = "tp"))

setMethod("groupNamesResults", "componentsTPs", function(obj)
{
    if (length(obj) == 0)
        return(character())
    return(union(componentInfo(obj)$parent_group, groupNames(obj)))
})

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

setMethod("parentsFromScreening", "componentsTPs", function(obj) obj@parentsFromScreening)
setMethod("TPsFromScreening", "componentsTPs", function(obj) obj@TPsFromScreening)

#' @describeIn componentsTPs Returns all component data as a \code{\link{data.table}}.
#' @export
setMethod("as.data.table", "componentsTPs", function(x, candidates = FALSE)
{
    checkmate::assertFlag(candidates)
    
    # get rid of double links column when merging components / componentInfo
    if (length(x) > 0)
        x@componentInfo <- x@componentInfo[, -"links"]
    
    ret <- callNextMethod(x)
    
    if (candidates)
        ret <- mergeTPComponCandidatesTab(ret)
    
    ret <- removeDTColumnsIfPresent(ret, c("candidates", "candidatesForm", "candidatesComp"))
    
    return(ret)
})

#' @describeIn componentsTPs Provides various rule based filtering options to clean and prioritize TP data.
#'
#' @param \dots,verbose Further arguments passed to the base \code{\link[=filter,components-method]{filter method}}.
#' @param retDirMatch If set to \code{TRUE}, only keep TPs for which the retention time direction (\code{retDir}, see
#'   Details in \link{componentsTPs}) matches with the observed direction. TPs will never be removed if the
#'   expected/observed direction is \samp{0} (\emph{i.e.} unknown or not significantly different than the parent).
#' @param minSpecSim,minSpecSimPrec,minSpecSimBoth The minimum spectral similarity of a TP compared to its parent
#'   (\samp{0-1}). The \code{minSpecSimPrec} and \code{minSpecSimBoth} apply to binned data that is shifted with the
#'   \code{"precursor"} and \code{"both"} method, respectively (see \link[=specSimParams]{MS spectral similarity
#'   parameters} for more details). Set to \code{NULL} to ignore.
#' @param minFragMatches,minNLMatches Minimum number of parent/TP fragment and neutral loss matches, respectively. Set
#'   to \code{NULL} to ignore. See the \verb{Linking parents and transformation products} section in
#'   \code{\link{generateComponentsTPs}} for more details.
#' @param formulas A \code{\link{formulas}} object. The formula annotation data in this object is to verify if elemental
#'   additions/subtractions from metabolic logic reactions are possible (hence, it only works with data from
#'   \code{\link{generateTPsLogic}}). To verify elemental additions, only TPs with at least one candidate formula that
#'   has these elements are kept. Similarly, for elemental subtractions, any of the parent candidate formulae must
#'   contain the subtraction elements. Note that TPs are currently not filtered if either the parent or the TP has no
#'   formula annotations. Set to \code{NULL} to ignore.
#' @param negate If \code{TRUE} then filters are applied in opposite manner.
#'
#' @return \code{filter} returns a filtered \code{componentsTPs} object.
#'
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
            ct <- minColFilter(ct, "totalFragmentMatches", minFragMatches)
            ct <- minColFilter(ct, "totalNeutralLossMatches", minNLMatches)
            
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

#' @describeIn componentsTPs Plots an interactive network graph for linked components. Components are linked with each
#'   other if one or more transformation products overlap. The graph is constructed with the \pkg{\link{igraph}} package
#'   and rendered with \pkg{\link{visNetwork}}.
#'
#' @inheritParams plotGraph,componentsNT-method
#' 
#' @template plotGraph
#'
#' @references \addCitations{igraph}{2}
#' 
#' @export
setMethod("plotGraph", "componentsTPs", function(obj, onlyLinked = TRUE, width = NULL, height = NULL)
{
    checkmate::assertFlag(onlyLinked)

    cInfo <- componentInfo(obj)
    cTable <- componentTable(obj)
    
    titles <- sprintf("<b>%s</b> (%s - %s)<br>TPs: <i>%s (%s)</i>",
                      names(obj), cInfo$parent_name, cInfo$parent_group,
                      sapply(cTable, function(cmp) paste0(cmp$TP_name, collapse = ", ")),
                      sapply(cTable, function(cmp) paste0(unique(cmp$group), collapse = ", ")))
    makeGraph(obj, onlyLinked, titles, width, height)
})

#' Generate components of transformation products
#'
#' Generates components by linking feature groups of transformation products and their parents.
#'
#' @templateVar algo transformation product screening
#' @templateVar do generate components
#' @templateVar generic generateComponents
#' @templateVar algoParam tp
#' @template algo_generator
#'
#' @details This method typically employs data from \link[=generateTPs]{generated transformation products} to find
#'   parents and their TPs. However, this data is not necessary, and components can also be made based on MS/MS
#'   similarity and/or other annotation similarities between the parent and its TPs. For more details see the
#'   \verb{Linking parents and transformation products} section below.
#'
#' @param fGroups The input \code{\link{featureGroups}} for componentization. See \code{fGroupsTPs}.
#' @param fGroupsTPs A \code{\link{featureGroups}} object containing the feature groups that are expected to be
#'   transformation products. If a distinction between parents and TPs is not yet known, \code{fGroupsTPs} should equal
#'   the \code{fGroups} argument. Otherwise, \code{fGroups} should only contain the parent feature groups, and both
#'   \code{fGroups} and \code{fGroupsTPs} \emph{must} be a subset of the same \code{\link{featureGroups}} object.
#' @param ignoreParents If \code{TRUE} then feature groups present in both \code{fGroups} and \code{fGroupsTPs} are not
#'   considered as TPs.
#' @param TPs A \code{\link{transformationProducts}} object. Set to \code{NULL} to perform linking without this data.
#' @param MSPeakLists,formulas,compounds A \code{\link{MSPeakLists}}/\code{\link{formulas}}/\code{\link{compounds}}
#'   object to calculate MS/MS or annotation similarities between parents and TPs. If \code{NULL} then this data is not
#'   calculated. For more details see the \verb{Linking parents and transformation products} section below.
#' @param minRTDiff Minimum retention time (in seconds) difference between the parent and a TP to determine whether a TP
#'   elutes prior/after the parent (to calculate \code{retDir} values, see Details in \link{componentsTPs}))
#'
#' @template specSimParams-arg
#'
#' @note The \code{shift} parameter of \code{specSimParams} is ignored by \code{generateComponentsTPs}, since it always
#'   calculates similarities with all supported options.
#'
#' @return The components are stored in objects derived from \code{\link{componentsTPs}}.
#'
#' @section Linking parents and transformation products: Each component consists of feature groups that are considered
#'   to be transformation products for one parent (the parent that 'belongs' to the component can be retrieved with the
#'   \code{\link{componentInfo}} method). The parent feature groups are taken from the \code{fGroups} parameter, while
#'   the feature groups for TPs are taken from \code{fGroupsTPs}. If a feature group occurs in both variables, it may
#'   therefore be considered as both a parent or TP.
#'
#'   If transformation product data is given, \emph{i.e.} the \code{TPs} argument is set, then a suspect screening of
#'   the TPs must be performed in advance (see \code{\link{screenSuspects}} and \code{\link{convertToSuspects}} to
#'   create the suspect list). Furthermore, if TPs were generated with \code{\link{generateTPsBioTransformer}} or
#'   \code{\link{generateTPsLibrary}} then the suspect screening must also include the parents (\emph{e.g.} by setting
#'   \code{includeParents=TRUE} when calling \code{convertToSuspects} or by amending results by setting
#'   \code{amend=TRUE} to \code{screenSuspects}). The suspect screening is necessary for the componentization algorithm
#'   to map the feature groups of the parent or TP. If the the suspect screening yields multiple TP hits, all will be
#'   reported. Similarly, if the suspect screening contains multiple hits for a parent, a component is made for each of
#'   the parent hits.
#'
#'   In case no transformation product data is provided (\code{TPs=NULL}), the componentization algorithm simply assumes
#'   that each feature group from \code{fGroupsTPs} is a potential TP for every parent feature group in \code{fGroups}.
#'   For this reason, it is highly recommended to specify which feature groups are parents/TPs (see the
#'   \code{fGroupsTPs} argument description above) and \emph{crucial} that the data is post-processed, for instance by
#'   only retaining TPs that have high annotation similarity with their parents (see the
#'   \code{\link[=filter,componentsTPs-method]{filter}} method for \code{\link{componentsTPs}}).
#'
#'   A typical way to distinguish which feature groups are parents or TPs from two different (groups of) samples is by
#'   calculating Fold Changes (see the \code{\link[=as.data.table,featureGroups-method]{as.data.table}} method for
#'   feature groups and \code{\link{plotVolcano}}). Of course, other statistical techniques from \R are also suitable.
#'
#'   During componentization, several characteristics are calculated which may be useful for post-processing: \itemize{
#'
#'   \item \code{specSimilarity}: the MS/MS spectral similarity between the feature groups of the TP and its parent
#'   (\samp{0-1}).
#'
#'   \item \code{specSimilarityPrec},\code{specSimilarityBoth}: as \code{specSimilarity}, but calculated with binned
#'   data using the \code{"precursor"} and \code{"both"} method, respectively (see \link[=specSimParams]{MS spectral
#'   similarity parameters} for more details).
#'
#'   \item \code{totalFragmentMatches} The number of MS/MS fragment formula annotations that overlap between the TP and
#'   parent. If both the \code{formulas} and \code{compounds} arguments are specified then the annotation data is pooled
#'   prior to calculation. Note that only unique matches are counted. Furthermore, note that annotations from \emph{all}
#'   candidates are considered, even if the formula/structure of the parent/TP is known. Hence, \code{totalFragmentMatches}
#'   is mainly useful when little or no chemical information is known on the parents/TPs, \emph{i.e.}, when
#'   \code{TPs=NULL} or originates from \code{\link{generateTPsLogic}}. Since annotations for all candidates are used,
#'   it is highly recommended that the annotation objects are first processed with the \code{\link{filter}} method, for
#'   instance, to select only the top ranked candidates.
#'
#'   \item \code{totalNeutralLossMatches} As \code{totalFragmentMatches}, but counting overlapping neutral loss formulae.
#'
#'   \item \code{retDir} The retention time direction of the TP relative to its parent. See Details in
#'   \link{componentsTPs}. If TP data was specified, the expected direction is stored in \code{TP_retDir}.
#'
#'   \item \code{retDiff},\code{mzDiff},\code{formulaDiff} The retention time, \emph{m/z} and formula difference between
#'   the parent and TP (latter only available if data TP formula is available).
#'
#'   }
#'
#' @section Sets workflows: In a \link[=sets-workflow]{sets workflow} the component tables are amended with extra
#'   information such as overall/specific set spectrum similarities. As sets data is mixed, transformation products are
#'   able to be linked with a parent, even if they were not measured in the same set.
#'
#' @templateVar what generateComponentsTPs
#' @template main-rd-method
#' @export
setMethod("generateComponentsTPs", "featureGroups", function(fGroups, fGroupsTPs = fGroups, TPs = NULL,
                                                             MSPeakLists = NULL, formulas = NULL, compounds = NULL,
                                                             ignoreParents = FALSE, minRTDiff = 20,
                                                             specSimParams = getDefSpecSimParams(), parallel = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertClass, . ~ fGroupsTPs + TPs + MSPeakLists + formulas + compounds,
           c("featureGroups", "transformationProducts", "MSPeakLists", "formulas", "compounds"),
           null.ok = c(FALSE, TRUE, TRUE, TRUE, TRUE), fixed = list(add = ac))
    checkmate::assertFlag(ignoreParents, add = add)
    checkmate::assertNumber(minRTDiff, lower = 0, finite = TRUE, add = add)
    assertSpecSimParams(specSimParams, add = ac)
    checkmate::assertFlag(parallel, add = ac)
    checkmate::reportAssertions(ac)
    
    if (!is.null(TPs) && (parentsFromScreening(TPs) || TPsFromScreening(TPs)) &&
        (!isScreening(fGroups) || !isScreening(fGroupsTPs)))
        stop("Input feature groups need to be screened for parents/TPs!")

    return(doGenComponentsTPs(fGroups, fGroupsTPs, TPs, MSPeakLists, formulas, compounds, ignoreParents, minRTDiff,
                              specSimParams, parallel))
})

#' @rdname generateComponentsTPs
#' @export
setMethod("generateComponentsTPs", "featureGroupsSet", function(fGroups, fGroupsTPs = fGroups, TPs = NULL,
                                                                MSPeakLists = NULL, formulas = NULL, compounds = NULL,
                                                                ignoreParents = FALSE, minRTDiff = 20,
                                                                specSimParams = getDefSpecSimParams(), parallel = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertClass, . ~ fGroupsTPs + TPs + MSPeakLists + formulas + compounds,
           c("featureGroupsSet", "transformationProducts", "MSPeakListsSet", "formulasSet", "compoundsSet"),
           null.ok = c(FALSE, TRUE, TRUE, TRUE, TRUE), fixed = list(add = ac))
    checkmate::assertFlag(ignoreParents, add = add)
    checkmate::assertNumber(minRTDiff, lower = 0, finite = TRUE, add = add)
    assertSpecSimParams(specSimParams, add = ac)
    checkmate::assertFlag(parallel, add = ac)
    checkmate::reportAssertions(ac)

    if (!is.null(TPs) && (parentsFromScreening(TPs) || TPsFromScreening(TPs)) &&
        (!isScreening(fGroups) || !isScreening(fGroupsTPs)))
        stop("Input feature groups need to be screened for parents/TPs!")
    
    ret <- doGenComponentsTPs(fGroups, fGroupsTPs, TPs, MSPeakLists, formulas, compounds, ignoreParents, minRTDiff,
                              specSimParams, parallel)
    
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
                    # NOTE: do _not_ use intersect below, since we want to keep duplicate group names (ie if multiple
                    # TPs are assigned to a same feat group) so that similarities are are returned for all.
                    grpsInSet <- cmp$group[cmp$group %in% groupNames(unsetMSPeakLists[[s]])]
                    
                    if (length(grpsInSet) > 0 && !is.null(unsetMSPeakLists[[s]][[parentFG]][["MSMS"]]))
                    {
                        sims <- genTPSpecSimilarities(unsetMSPeakLists[[s]], parentFG, grpsInSet,
                                                      specSimParams = specSimParams)
                        cmp[group %in% grpsInSet, (simColNames) := sims]
                    }
                    else
                        cmp[, (simColNames) := NA_real_]
                    
                    if (!is.null(unsetFormulas[[s]]) || !is.null(unsetCompounds[[s]]))
                    {
                        # calculate per set spectrum similarities
                        sims <- genTPAnnSimilarities(parentFG, cmp$group, unsetMSPeakLists[[s]], unsetFormulas[[s]],
                                                     unsetCompounds[[s]])
                        annColNames <- paste0(c("totalFragmentMatches", "totalNeutralLossMatches"), "-", s)
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
