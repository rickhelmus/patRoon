#' @include main.R
#' @include TP.R
NULL

#' @export
transformationProductsStructure <- setClass("transformationProductsStructure", contains = "transformationProducts")

setMethod("initialize", "transformationProductsStructure", function(.Object, calcSims, fpType, fpSimMethod, ...)
{
    .Object <- callNextMethod(.Object, ...)
    
    if (length(.Object) > 0 && calcSims)
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

#' @describeIn transformationProductsStructure Performs rule-based filtering of the \command{BioTransformer}
#'   predictions. Useful to simplify and clean-up the data.
#'
#' @param \dots Further argument passed to the base \code{\link[=filter,transformationProducts-method]{filter method}}.
#' @param removeDuplicates If \code{TRUE} then the TPs of a parent with duplicate structures (\acronym{SMILES}) are
#'   removed. Such duplicates may occur when different transformation pathways yield the same TPs. The first TP
#'   candidate with duplicate structure will be kept.
#' @param removeParentIsomers If \code{TRUE} then TPs with an equal formula as their parent (isomers) are removed.
#' @param removeTPIsomers If \code{TRUE} then all TPs with equal formula as any sibling TPs (isomers) are removed.
#'   Unlike \code{removeDuplicates}, \emph{all} TP candidates are removed (including the first match). This filter
#'   automatically sets \code{removeDuplicates=TRUE} to avoid complete removal of TPs with equal structure.
#' @param minSimilarity Minimum structure similarity (\samp{0-1}) that a TP should have relative to its parent. For
#'   details on how these similarities are calculated, see the \code{\link{generateTPsBioTransformer}} function. May be
#'   useful under the assumption that parents and TPs who have a high structural similarity, also likely have a high
#'   MS/MS spectral similarity (which can be evaluated after componentization with \code{\link{generateComponentsTPs}}.
#' @param negate If \code{TRUE} then filters are performed in opposite manner.
#'
#' @return \code{filter} returns a filtered \code{transformationProductsStructure} object.
#'
#' @export
setMethod("filter", "transformationProductsStructure", function(obj, ..., removeParentIsomers = FALSE,
                                                                removeTPIsomers = FALSE, removeDuplicates = FALSE,
                                                                minSimilarity = NULL, negate = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertFlag, . ~ removeParentIsomers + removeTPIsomers + removeDuplicates + negate,
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
            # NOTE: obj@products should be first arg to Map to keep names...
            obj@products <- Map(obj@products, parents(obj)$formula, f = function(prod, pform)
            {
                if (removeDuplicates)
                    prod <- if (negate) prod[duplicated(SMILES)] else prod[!duplicated(SMILES)]
                if (removeParentIsomers)
                    prod <- if (negate) prod[formula == pform] else prod[formula != pform]
                if (removeTPIsomers)
                {
                    df <- getDuplicatedStrings(prod$formula)
                    prod <- if (negate) prod[formula %chin% df] else prod[!formula %chin% df]
                }
                return(prod)
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
            obj@products <- lapply(obj@products, function(p) p[pred(similarity)])
        }
        
        obj@products <- pruneList(obj@products, checkZeroRows = TRUE)
        obj@parents <- obj@parents[name %in% names(obj@products)]
        
        saveCacheData("filterTPs", obj, hash)
    }
    
    if (...length() > 0)
        obj <- callNextMethod(obj, ..., negate = negate)
    
    newn <- length(obj)
    printf("Done! Filtered %d (%.2f%%) TPs. Remaining: %d\n", oldn - newn, if (oldn == 0) 0 else (1-(newn/oldn))*100, newn)
    
    return(obj)
})

#' @export
setMethod("plotGraph", "transformationProductsStructure", function(obj, which, components = NULL, structuresMax = 25,
                                                                   prune = TRUE, onlyCompletePaths = FALSE)
{
    # UNDONE: don't make name unique, but use IDs?
    
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
    
    obj <- obj[which]
    if (length(obj) == 0)
        stop("No TPs to plot", call. = FALSE)
    
    TPTab <- copy(as.data.table(obj))
    TPTab[, c("name_orig", "name") := .(name, make.unique(name))]
    TPTab[, parent_name := fifelse(is.na(parent_ID), parent, name[match(parent_ID, ID)]), by = "parent"]
    
    if (!is.null(components))
    {
        cmpTab <- as.data.table(components)
        TPTab <- TPTab[parent %chin% cmpTab$parent_name] # omit missing root parents
        TPTab[, present := name_orig %chin% cmpTab$TP_name]
        
        TPTab[, childPresent := FALSE]
        markChildPresent <- function(TPNames)
        {
            if (length(TPNames) == 0)
                return()
            TPTab[name %chin% TPNames, childPresent := TRUE]
            pars <- TPTab[name %chin% TPNames]$parent_name
            markChildPresent(pars[TPTab[name %chin% pars]$childPresent == FALSE])
        }
        markChildPresent(TPTab[present == TRUE]$parent_name)
        
        if (prune)
            TPTab <- TPTab[present == TRUE | childPresent == TRUE]
        if (onlyCompletePaths)
        {
            TPTab <- TPTab[present == TRUE]
            # keep removing TPs without parent until no change
            oldn <- nrow(TPTab)
            repeat
            {
                TPTab <- TPTab[parent_name == parent | parent_name %chin% name]
                newn <- nrow(TPTab)
                if (oldn == newn)
                    break
                oldn <- newn
            }
        }
    }
    
    pars <- parents(obj)
    TPTab[, parent_formula := fifelse(is.na(parent_ID),
                                      pars$formula[match(parent_name, pars$name)],
                                      formula[match(parent_name, name)])]
    TPTab[, formulaDiff := mapply(formula, parent_formula, FUN = function(f, pf)
    {
        sfl <- splitFormulaToList(subtractFormula(f, pf))
        ret <- ""
        subfl <- sfl[sfl < 0]
        if (length(subfl) > 0)
            ret <- paste0("-", formulaListToString(abs(subfl)))
        addfl <- sfl[sfl > 0]
        if (length(addfl) > 0)
            ret <- if (nzchar(ret)) paste0(ret, " +", formulaListToString(addfl)) else paste0("+", formulaListToString(addfl))
        return(ret)
    })]
    
    nodes <- data.table(id = union(TPTab$parent, TPTab$name))
    nodes[, isTP := id %chin% TPTab$name]
    nodes[isTP == TRUE, label := paste0("TP", TPTab$chem_ID[match(id, TPTab$name)])]
    nodes[isTP == FALSE, label := id]
    nodes[, group := if (.N > 1) label else "unique", by = "label"]
    nodes[, present := isTP == FALSE | TPTab$present[match(id, TPTab$name)]]
    nodes[present == TRUE, shapeProperties := list(list(list(useBorderWithImage = TRUE)))]
    nodes[present == FALSE, shapeProperties := list(list(list(useBorderWithImage = FALSE)))]
    nodes[, present := NULL]
    nodes[isTP == FALSE, level := 0]
    nodes[isTP == TRUE, level := TPTab$generation[match(id, TPTab$name)]]
    
    if (nrow(nodes) <= structuresMax && nrow(nodes) > 0)
    {
        # UNDONE: make util?
        imgf <- tempfile(fileext = ".png") # temp file is re-used
        getURIFromSMILES <- function(SMILES)
        {
            mol <- getMoleculesFromSMILES(SMILES, emptyIfFails = TRUE)[[1]]
            withr::with_png(imgf, withr::with_par(list(mar = rep(0, 4)), plot(getRCDKStructurePlot(mol, 150, 150))))
            return(knitr::image_uri(imgf))
        }
        nodes[, shape := "image"]
        nodes[, SMILES := fifelse(isTP, TPTab$SMILES[match(id, TPTab$name)], pars$SMILES[match(id, pars$name)])]
        nodes[, image := getURIFromSMILES(SMILES[1]), by = "SMILES"]
        nodes[, SMILES := NULL]
    }
    else
        nodes[, shape := "ellipse"]
    
    TPCols <- intersect(c("name", "name_lib", "SMILES", "formula", "generation", "accumulation", "production",
                          "globalAccumulation", "likelihood", "Lipinski_Violations", "Insecticide_Likeness_Violations",
                          "Post_Em_Herbicide_Likeness_Violations", "transformation", "transformation_ID", "enzyme",
                          "biosystem", "evidencedoi", "evidencedref", "sourcecomment", "datasetref", "similarity",
                          "mergedBy", "coverage"), names(TPTab))
    nodes[isTP == TRUE, title := sapply(id, function(TP)
    {
        TPTabSub <- TPTab[name == TP, TPCols, with = FALSE]
        return(paste0(names(TPTabSub), ": ", TPTabSub, collapse = "<br>"))
    })]
    
    edges <- data.table(from = TPTab$parent_name, to = TPTab$name, label = TPTab$formulaDiff)
    
    visNetwork::visNetwork(nodes = nodes, edges = edges) %>%
        visNetwork::visNodes(shapeProperties = list(useBorderWithImage = FALSE)) %>%
        visNetwork::visEdges(arrows = "to", font = list(align = "top", size = 12)) %>%
        visNetwork::visOptions(highlightNearest = list(enabled = TRUE, hover = TRUE, algorithm = "hierarchical"),
                               selectedBy = list(variable = "group", main = "Select duplicate TPs",
                                                 values = unique(nodes$group[nodes$group != "unique"]))) %>%
        visNetwork::visHierarchicalLayout(enabled = TRUE, sortMethod = "directed")
})

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

transformationProductsStructureConsensus <- setClass("transformationProductsStructureConsensus",
                                                     contains = "transformationProductsStructure")

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
                                                                         with = FALSE]), idcol = "mergedBy")
    allParentsTab[, mergedBy := paste0(unique(mergedBy), collapse = ","), by = "name"]
    allParentsTab <- unique(allParentsTab, by = "name")
    allParentsTab <- allParentsTab[match(name, names(mergedTPs))] # sync order
    setcolorder(allParentsTab, parCols) # move mergedBy to end
    
    # HACK: dummy values for fpType/fpSimMethod, they are not really used as calcSims is always FALSE
    return(transformationProductsStructureConsensus(calcSims = FALSE, fpType = "extended", fpSimMethod = "tanimoto",
                                                    parents = allParentsTab, products = mergedTPs,
                                                    algorithm = paste0(unique(sapply(allTPs, algorithm)), collapse = ",")))
})
