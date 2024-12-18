# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
#' @include TP.R
NULL

#' Base transformation products (TP) class with structure information
#'
#' Holds information for all TPs for a set of parents, including structural information.
#'
#' This (virtual) class is derived from the \code{\link{transformationProducts}} base class, please see its
#' documentation for more details. Objects from this class are returned by \link[=generateTPs]{TP generators}. More
#' specifically, algorithms that works with chemical structures (\emph{e.g.} \code{biotransformer}), uses this class to
#' store their results. The methods defined for this class extend the functionality for the base
#' \code{\link{transformationProducts}} class.
#'
#' @param obj,TPs \code{transformationProductsStructure} derived object to be accessed
#' @param commonParents Only consider TPs from parents that are common to all compared objects.
#' @param \dots For \code{filter}: Further argument passed to the base
#'   \code{\link[=filter,transformationProducts-method]{filter method}}.
#'
#'   For \code{plotVenn}, \code{plotUpSet} and \code{consensus}: further (unique) \code{transformationProductsStructure}
#'   objects.
#'
#' @section Comparison between objects: The methods that compare different objects (\emph{e.g.} \code{plotVenn} and
#'   \code{consensus}) use the \acronym{InChIKey} to match TPs between objects. Moreover, the parents between objects
#'   are matched by their name. Hence, it is \emph{crucial} that the input parents to \code{\link{generateTPs}}
#'   (\emph{i.e.} the \code{parents} argument) are named equally.
#'
#' @seealso The base class \code{\link{transformationProducts}} for more relevant methods and \code{\link{generateTPs}}
#'
#' @templateVar class transformationProductsStructure
#' @template class-hierarchy
#'
#' @export
transformationProductsStructure <- setClass("transformationProductsStructure",
                                            contains = c("VIRTUAL", "transformationProducts"))

setMethod("initialize", "transformationProductsStructure", function(.Object, calcSims, fpType, fpSimMethod, ...)
{
    .Object <- callNextMethod(.Object, ...)
    
    if (length(.Object) == 0)
        return(.Object)
    
    # remove neutralized TPs that became duplicates
    rmNeutTPs <- 0
    .Object@products <- lapply(products(.Object), function(pr)
    {
        if (!is.null(pr[["molNeutralized"]]))
        {
            prNN <- pr[molNeutralized == FALSE]
            wh <- pr$molNeutralized & pr$InChIKey %chin% prNN$InChIKey
            pr <- pr[!wh]
            rmNeutTPs <<- rmNeutTPs + sum(wh)
        }
        return(pr)
    })
    
    if (rmNeutTPs > 0)
        printf(sprintf("Removed %d TPs that were duplicates after neutralization.\n", rmNeutTPs))
    
    if (calcSims)
    {
        hash <- makeHash(.Object)
        cd <- loadCacheData("TPsParentSims", hash)
        if (!is.null(cd))
            .Object <- cd
        else
        {
            printf("Calculating parent/TP structural similarities...\n")
            pspl <- split(parents(.Object), seq_len(nrow(parents(.Object))))
            .Object@products <- withProg(nrow(parents(.Object)), FALSE, Map(products(.Object), pspl, f = function(pr, pa)
            {
                pr <- copy(pr)
                pr[, similarity := sapply(SMILES, distSMILES, SMI1 = pa$SMILES, fpType = fpType, fpSimMethod = fpSimMethod)]
                doProgress()
                return(pr)
            }))
            saveCacheData("TPsParentSims", .Object, hash)
        }
    }
    
    return(.Object)
})

#' @templateVar class transformationProductsStructure
#' @template convertToMFDB
#' @export
setMethod("convertToMFDB", "transformationProductsStructure", function(TPs, out, includeParents = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertPathForOutput(out, overwrite = TRUE, add = ac) # NOTE: assert doesn't work on Windows...
    checkmate::assertFlag(includeParents, add = ac)
    checkmate::reportAssertions(ac)

    prodAll <- rbindlist(products(TPs), idcol = "parent")

    doConvertToMFDB(prodAll, parents(TPs), out, includeParents)
})

setMethod("linkParentsToFGroups", "transformationProductsStructure", function(TPs, fGroups)
{
    return(screenInfo(fGroups)[name %in% names(TPs), c("name", "group"), with = FALSE])
})

#' @describeIn transformationProductsStructure Performs rule-based filtering. Useful to simplify and clean-up the data.
#'
#' @param removeDuplicates If \code{TRUE} then the TPs of a parent with duplicate structures (\acronym{SMILES}) are
#'   removed. Such duplicates may occur when different transformation pathways yield the same TPs. The first TP
#'   candidate with duplicate structure will be kept.
#' @param removeParentIsomers If \code{TRUE} then TPs with an equal formula as their parent (isomers) are removed.
#' @param removeTPIsomers If \code{TRUE} then all TPs with equal formula as any sibling TPs (isomers) are removed.
#'   Unlike \code{removeDuplicates}, \emph{all} TP candidates are removed (including the first match). This filter
#'   automatically sets \code{removeDuplicates=TRUE} so that TPs are only removed if with different structure.
#' @param minSimilarity Minimum structure similarity (\samp{0-1}) that a TP should have relative to its parent. This
#'   data is only available if the \code{calcSims} argument to \code{\link{generateTPs}} was set to \code{TRUE}. May be
#'   useful under the assumption that parents and TPs who have a high structural similarity, also likely have a high
#'   MS/MS spectral similarity (which can be evaluated after componentization with \code{\link{generateComponentsTPs}}.
#'   Any values that are \code{NA} are removed (which only occur when a consensus was made from objects that not all
#'   have similarity information).
#'
#' @inheritParams filter,transformationProducts-method
#'
#' @return \code{filter} returns a filtered \code{transformationProductsStructure} object.
#'
#' @export
setMethod("filter", "transformationProductsStructure", function(obj, ..., removeParentIsomers = FALSE,
                                                                removeTPIsomers = FALSE, removeDuplicates = FALSE,
                                                                minSimilarity = NULL, verbose = TRUE, negate = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertFlag, . ~ removeParentIsomers + removeTPIsomers + removeDuplicates + verbose + negate,
           fixed = list(add = ac))
    checkmate::assertNumber(minSimilarity, lower = 0, finite = TRUE, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (length(obj) == 0)
        return(obj)
    
    if (removeTPIsomers)
        removeDuplicates <- TRUE
    
    oldn <- length(obj)
    
    hash <- makeHash(obj, removeParentIsomers, removeTPIsomers, removeDuplicates, minSimilarity, negate)
    cache <- loadCacheData("filterTPs", hash)
    if (!is.null(cache))
        obj <- cache
    else
    {
        if (removeParentIsomers || removeTPIsomers || removeDuplicates)
        {
            mark <- if (negate) function(x) !x else function(x) x
            obj <- delete(obj, j = function(tab, par)
            {
                tab <- copy(tab)
                tab[, keep := TRUE]
                
                if (removeDuplicates)
                    tab[keep == TRUE, keep := mark(!duplicated(SMILES))]
                if (removeParentIsomers)
                {
                    pform <- parents(obj)[match(par, name)]$formula
                    tab[keep == TRUE, keep := mark(formula != pform)]
                }
                if (removeTPIsomers)
                {
                    df <- getDuplicatedStrings(tab[keep == TRUE]$formula)
                    tab[keep == TRUE, keep := mark(!formula %chin% df)]
                }
                return(!tab$keep)
            })
        }
        
        if (!is.null(minSimilarity) && length(obj) > 0)
        {
            if (is.null(obj[[1]][["similarity"]]))
                stop("Cannot filter on structural similarities: no similarities were calculated. ",
                     "Please set calcSims=TRUE when calling generateTPs().", call. = FALSE)
            pred <- if (negate)
                function(x) is.na(x) | x < minSimilarity
            else
                function(x) !is.na(x) & numGTE(x, minSimilarity)
            obj <- delete(obj, j = function(tab, ...) !pred(tab$similarity))
        }
        
        obj@products <- pruneList(obj@products, checkZeroRows = TRUE)
        obj@parents <- obj@parents[name %in% names(obj@products)]
        
        saveCacheData("filterTPs", obj, hash)
    }
    
    if (...length() > 0)
        obj <- callNextMethod(obj, ..., verbose = FALSE, negate = negate)
    
    if (verbose)
    {
        newn <- length(obj)
        printf("Done! Filtered %d (%.2f%%) TPs. Remaining: %d\n", oldn - newn, if (oldn == 0) 0 else (1-(newn/oldn))*100, newn)
    }
    
    return(obj)
})

#' @templateVar class transformationProductsStructure
#' @template plotGraph-TPs
#' 
#' @param structuresMax An \code{integer} with the maximum number of structures to plot. Setting a maximum is mainly
#'   done to avoid long times needed to construct the graph.
#'   
#' @template plotGraph
#'
#' @export
setMethod("plotGraph", "transformationProductsStructure", function(obj, which, components = NULL, structuresMax = 25,
                                                                   prune = TRUE, onlyCompletePaths = FALSE,
                                                                   width = NULL, height = NULL)
{
    checkmate::assert(
        checkmate::checkSubset(which, names(obj), empty.ok = FALSE),
        checkmate::checkIntegerish(which, lower = 1, upper = nrow(parents(obj)), any.missing = FALSE, min.len = 1),
        .var.name = "which"
    )
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(components, "componentsTPs", null.ok = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ prune + onlyCompletePaths, fixed = list(add = ac))
    checkmate::assertCount(structuresMax, add = ac)
    checkmate::reportAssertions(ac)

    doPlotTPGraph(as.data.table(obj[which]), parents(obj),
                  cmpTab = if (!is.null(components)) as.data.table(components) else NULL,
                  structuresMax = structuresMax, prune = prune, onlyCompletePaths = onlyCompletePaths,
                  width = width, height = height)
})

#' @describeIn transformationProductsStructure plots a Venn diagram (using \pkg{\link{VennDiagram}}) outlining unique and shared
#'   candidates of up to five different \code{featureAnnotations} objects.
#'
#' @inheritParams plotVenn,featureAnnotations-method
#'
#' @template plotvenn-ret
#' 
#' @export
setMethod("plotVenn", "transformationProductsStructure", function(obj, ..., commonParents = FALSE,
                                                                  labels = NULL, vennArgs = NULL)
{
    # UNDONE: this method is mostly a copy of the featureAnnotations method, merge?
    
    allTPs <- c(list(obj), list(...))
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertList(allTPs, types = "transformationProductsStructure", min.len = 2, any.missing = FALSE,
                          unique = TRUE, .var.name = "...", add = ac)
    checkmate::assertFlag(commonParents, add = ac)
    checkmate::assertCharacter(labels, min.chars = 1, len = length(allTPs), null.ok = TRUE, add = ac)
    checkmate::assertList(vennArgs, names = "unique", null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (is.null(labels))
        labels <- make.unique(sapply(allTPs, algorithm))
    if (is.null(vennArgs))
        vennArgs <- list()

    if (commonParents)
    {
        commonPars <- Reduce(intersect, lapply(allTPs, names))
        allTPs <- lapply(allTPs, "[", commonPars)
    }
    
    allTPTabs <- lapply(allTPs, as.data.table)
    unLengths <- sapply(allTPs, function(TPs) if (length(TPs) == 0) 0 else sum(sapply(products(TPs),
                                                                                      function(p) uniqueN(p$InChIKey))))
    do.call(makeVennPlot, c(list(allTPTabs, labels, unLengths, function(obj1, obj2)
    {
        if (length(obj1) == 0 || length(obj2) == 0)
            return(data.table())
        fintersect(obj1[, c("parent", "InChIKey")], obj2[, c("parent", "InChIKey")])
    }, nrow), vennArgs))
})

#' @describeIn transformationProductsStructure Plots an UpSet diagram (using the \code{\link[UpSetR]{upset}} function)
#'   outlining unique and shared TPs between different \code{transformationProductsStructure} objects.
#' @templateVar withArgs TRUE
#' @template plotUpSet
#' @export
setMethod("plotUpSet", "transformationProductsStructure", function(obj, ..., commonParents = FALSE, labels = NULL,
                                                                   nsets = length(list(...)) + 1, nintersects = NA,
                                                                   upsetArgs = NULL)
{
    # UNDONE: this method is mostly a copy of the featureAnnotations method, merge?
    
    allTPs <- c(list(obj), list(...))
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertList(allTPs, types = "transformationProductsStructure", min.len = 2, any.missing = FALSE,
                          unique = TRUE, .var.name = "...", add = ac)
    checkmate::assertFlag(commonParents, add = ac)
    checkmate::assertCharacter(labels, min.chars = 1, len = length(allTPs), null.ok = TRUE, add = ac)
    checkmate::assertList(upsetArgs, names = "unique", null.ok = TRUE, add = ac)
    checkmate::assertCount(nsets, positive = TRUE)
    checkmate::assertCount(nintersects, positive = TRUE, na.ok = TRUE)
    checkmate::reportAssertions(ac)
    
    if (is.null(labels))
        labels <- make.unique(sapply(allTPs, algorithm))
    
    if (commonParents)
    {
        commonPars <- Reduce(intersect, lapply(allTPs, names))
        allTPs <- lapply(allTPs, "[", commonPars)
    }
    
    allTPTabs <- mapply(allTPs, labels, SIMPLIFY = FALSE, FUN = function(f, l)
    {
        ret <- as.data.table(f)
        if (nrow(ret) == 0)
            ret <- data.table(parent = character(), InChIKey = character())
        ret <- unique(ret[, c("parent", "InChIKey")])[, (l) := 1]
    })
    
    TPTab <- Reduce(function(f1, f2)
    {
        merge(f1, f2, by = c("parent", "InChIKey"), all = TRUE)
    }, allTPTabs)
    
    TPTab <- TPTab[, labels, with = FALSE]
    for (j in seq_along(TPTab))
        set(TPTab, which(is.na(TPTab[[j]])), j, 0)
    
    if (sum(sapply(TPTab, function(x) any(x>0))) < 2)
        stop("Need at least two non-empty objects to plot")
    
    do.call(UpSetR::upset, c(list(TPTab, nsets = nsets, nintersects = nintersects), upsetArgs))
})

#' @rdname transformationProductsStructure-class
transformationProductsStructureConsensus <- setClass("transformationProductsStructureConsensus",
                                                     contains = "transformationProductsStructure")

#' @describeIn transformationProductsStructure Generates a consensus from different
#'   \code{transformationProductsStructure} objects. Currently this removes any hierarchical data, and all TPs are
#'   considered to originate from the same (original) parent.
#'
#' @templateVar what TPs
#' @template consensus-common-args
#'
#' @param labels A \code{character} with names to use for labelling. If \code{NULL} labels are automatically generated.
#'
#' @return \code{consensus} returns a \code{transformationProductsStructure} object that is produced by merging results
#'   from multiple \code{transformationProductsStructure} objects.
#'
#' @note \code{consensus}: If the \code{retDir} values differs between matched TPs it will be set to \samp{0}. If
#'   structure similarity data is available (\emph{i.e.} \code{calcSims=TRUE} to \code{generateTPs}) then the mean
#'   similarity is calculated.
#'
#' @export
setMethod("consensus", "transformationProductsStructure", function(obj, ..., absMinAbundance = NULL,
                                                                   relMinAbundance = NULL, uniqueFrom = NULL,
                                                                   uniqueOuter = FALSE, labels = NULL)
{
    allTPs <- c(list(obj), list(...))
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertList(allTPs, types = "transformationProductsStructure", min.len = 2, any.missing = FALSE,
                          unique = TRUE, .var.name = "...", add = ac)
    checkmate::assertCharacter(labels, min.chars = 1, len = length(allTPs), null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (is.null(labels))
        labels <- make.unique(sapply(allTPs, algorithm))
    
    assertConsCommonArgs(absMinAbundance, relMinAbundance, uniqueFrom, uniqueOuter, labels)
    
    relMinAbundance <- max(NULLToZero(absMinAbundance) / length(allTPs), NULLToZero(relMinAbundance))
    if (!is.null(uniqueFrom) && !is.character(uniqueFrom))
        uniqueFrom <- labels[uniqueFrom]
    names(allTPs) <- labels
    
    # merge all TPs in a large table, collapse hierarchies, only keep algorithm independent columns
    # NOTE: name and IDs will be re-assigned later
    TPCols <- c("parent", "SMILES", "InChI", "InChIKey", "formula", "neutralMass", "retDir", "similarity")
    allTPsTab <- rbindlist(lapply(allTPs, function(TPs)
    {
        tab <- as.data.table(TPs)
        return(unique(tab[, intersect(TPCols, names(tab)), with = FALSE]))
    }), idcol = "mergedBy", fill = TRUE)
    
    if (nrow(allTPsTab) == 0)
        mergedTPs <- list()
    else
    {
        byTPCols <- c("parent", "InChIKey")
        
        # set retDir to 0 if there are conflicts
        allTPsTab[, retDir := if (!allSame(retDir)) 0 else retDir, by = byTPCols]
        
        # just average similarities
        if (!is.null(allTPsTab[["similarity"]]))
            allTPsTab[, similarity := if (all(is.na(similarity))) NA_real_ else mean(similarity, na.rm = TRUE),
                      by = byTPCols]
        
        if (!is.null(uniqueFrom))
        {
            allTPsTab[, keep := all(mergedBy %chin% uniqueFrom) && (!uniqueOuter || .N == 1), by = byTPCols]
            allTPsTab <- allTPsTab[keep == TRUE][, keep := NULL]
        }
        
        allTPsTab[, coverage := uniqueN(mergedBy) / length(allTPs), by = byTPCols]
        
        allTPsTab[, mergedBy := paste0(unique(mergedBy), collapse = ","), by = byTPCols]
        allTPsTab <- unique(allTPsTab, by = byTPCols)
        
        if (relMinAbundance > 0)
            allTPsTab <- allTPsTab[coverage >= relMinAbundance]
        
        allTPsTab[, c("parent_ID", "generation") := .(NA_integer_, 1)] # just assume all TPs are from the same parent
        allTPsTab[, ID := seq_len(.N), by = "parent"]
        allTPsTab[, c("chem_ID", "name") := .(ID, paste0(parent, "-TP", ID))]
        
        setcolorder(allTPsTab, c("name", "ID", "parent_ID", "chem_ID"))
        setcolorder(allTPsTab, setdiff(c(names(allTPsTab), "generation"), c("coverage", "mergedBy"))) # move coverage/mergedBy to end
        
        mergedTPs <- split(allTPsTab, by = "parent", keep.by = FALSE)
    }
    
    # same for parents: combine tables, only keep columns not specific to screening, remove duplicates and re-order
    parCols <- c("name", "SMILES", "InChI", "InChIKey", "formula", "neutralMass")
    allParentsTab <- rbindlist(lapply(allTPs, function(TPs) parents(TPs)[, intersect(parCols, names(parents(TPs))),
                                                                         with = FALSE]), idcol = "mergedBy", fill = TRUE)
    allParentsTab[, mergedBy := paste0(unique(mergedBy), collapse = ","), by = "name"]
    allParentsTab <- unique(allParentsTab, by = "name")
    allParentsTab <- allParentsTab[match(name, names(mergedTPs))] # sync order
    setcolorder(allParentsTab, parCols) # move mergedBy to end
    
    # HACK: dummy values for fpType/fpSimMethod, they are not really used as calcSims is always FALSE
    return(transformationProductsStructureConsensus(calcSims = FALSE, fpType = "extended", fpSimMethod = "tanimoto",
                                                    parents = allParentsTab, products = mergedTPs,
                                                    algorithm = paste0(unique(sapply(allTPs, algorithm)), collapse = ",")))
})
