#' @include main.R
#' @include workflow-step.R
NULL

featureAnnotations <- setClass("featureAnnotations",
                               slots = c(groupAnnotations = "list", scoreTypes = "character", scoreRanges = "list"),
                               contains = c("workflowStep", "VIRTUAL"))

setMethod("initialize", "featureAnnotations", function(.Object, ...)
{
    .Object <- callNextMethod(.Object, ...)
    .Object@groupAnnotations <- makeEmptyListNamed(.Object@groupAnnotations)
    return(.Object)
})

setMethod("annotations", "featureAnnotations", function(obj) obj@groupAnnotations)

setMethod("groupNames", "featureAnnotations", function(obj) names(obj@groupAnnotations))

setMethod("length", "featureAnnotations", function(x) if (length(x@groupAnnotations) > 0) sum(sapply(x@groupAnnotations, nrow)) else 0)

setMethod("[", c("featureAnnotations", "ANY", "missing", "missing"), function(x, i, ...)
{
    if (!missing(i))
    {
        i <- assertSubsetArgAndToChr(i, groupNames(x))
        x <- delete(x, setdiff(groupNames(x), i))
    }
    
    return(x)
})

setMethod("[[", c("featureAnnotations", "ANY", "missing"), function(x, i, j)
{
    assertExtractArg(i)
    return(x@groupAnnotations[[i]])
})

setMethod("$", "featureAnnotations", function(x, name)
{
    eval(substitute(x@groupAnnotations$NAME_ARG, list(NAME_ARG = name)))
})

setMethod("as.data.table", "featureAnnotations", function(x, fGroups = NULL, fragments = FALSE, countElements = NULL,
                                                          countFragElements = NULL, OM = FALSE,
                                                          normalizeScores = "none",
                                                          excludeNormScores = defaultExclNormScores(x))
{
    # NOTE: keep args in sync with formulas/compounds methods
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", null.ok = TRUE, add = ac)
    checkmate::assertCharacter(countElements, min.chars = 1, any.missing = FALSE, null.ok = TRUE, add = ac)
    checkmate::assertCharacter(countFragElements, min.chars = 1, any.missing = FALSE, null.ok = TRUE, add = ac)
    checkmate::assertFlag(OM, add = ac)
    checkmate::assertFlag(fragments, add = ac)
    assertNormalizationMethod(normalizeScores, add = ac)
    checkmate::assertCharacter(excludeNormScores, min.chars = 1, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (!is.null(countFragElements))
        fragments <- TRUE
    
    mcn <- mergedConsensusNames(x)
    annTable <- annotations(x)
    if (normalizeScores != "none")
    {
        annTable <- Map(annTable, x@scoreRanges, f = normalizeAnnScores,
                        MoreArgs = list(scoreCols = annScoreNames(x, TRUE), mConsNames = mcn, normalizeScores == "minmax",
                                      exclude = excludeNormScores))
    }
    
    if (fragments)
    {
        ret <- rbindlist(lapply(annTable, function(ct)
        {
            ct <- copy(ct)
            ct[, row := seq_len(nrow(ct))]
            
            fragTab <- rbindlist(ct$fragInfo, idcol = "row", fill = TRUE)
            fragTab[, PLIndex := NULL]
            cnames <- setdiff(names(fragTab), "row")
            setnames(fragTab, cnames, paste0("frag_", cnames))
            
            return(merge(ct, fragTab, by = "row", all.x = TRUE)[, -"row"])
        }), idcol = "group", fill = TRUE)
    }
    else
        ret <- rbindlist(annTable, idcol = "group", fill = TRUE)
    
    if (!is.null(fGroups))
    {
        ret[, c("ret", "group_mz") := groupInfo(fGroups)[group, c("rts", "mzs")]]
        setcolorder(ret, c("group", "ret", "group_mz"))
    }
    
    ret <- addElementInfoToAnnTable(ret, countElements, countFragElements, OM, TRUE)
    
    if (!is.null(ret[["fragInfo"]]))
        return(ret[, -"fragInfo"]) # not there if empty results
    
    return(ret)
})

#' @export
setMethod("delete", "featureAnnotations", function(obj, i = NULL, j = NULL, ...)
{
    # NOTE: this is ~ a c/p from the features method
    
    ac <- checkmate::makeAssertCollection()
    i <- assertDeleteArgAndToChr(i, groupNames(obj), add = ac)
    checkmate::assert(
        checkmate::checkIntegerish(j, any.missing = FALSE, null.ok = TRUE),
        checkmate::checkFunction(j, null.ok = TRUE),
        .var.name = "j"
    )
    checkmate::reportAssertions(ac)
    
    if (length(i) == 0 || (!is.null(j) && length(j) == 0))
        return(obj) # nothing to remove...
    
    # UNDONE: NULL for i and j will remove all?
    
    # i = NULL; j = vector: remove from all groups
    # i = vector; j = NULL: remove specified groups
    # j = function: remove specific results from given groups (or all if i=NULL)
    
    if (!is.function(j))
    {
        if (is.null(j))
            obj@groupAnnotations <- obj@groupAnnotations[setdiff(groupNames(obj), i)]
        else
        {
            obj@groupAnnotations[i] <- lapply(obj@groupAnnotations[i], function(at)
            {
                inds <- j[j <= nrow(at)]
                return(if (length(inds) > 0) at[-inds] else at)
            })
        }
    }
    else
    {
        obj@groupAnnotations[i] <- Map(obj@groupAnnotations[i], i, f = function(at, grp)
        {
            rm <- j(at, grp, ...)
            if (is.logical(rm))
                return(at[!rm])
            return(at[setdiff(seq_len(nrow(at)), rm)])
        })
    }
    
    obj@groupAnnotations <- pruneList(obj@groupAnnotations, checkZeroRows = TRUE)
    
    obj@scoreRanges <- obj@scoreRanges[names(obj@groupAnnotations)]

    return(obj)
})

setMethod("filter", "featureAnnotations", function(obj, minExplainedPeaks = NULL, scoreLimits = NULL, elements = NULL,
                                                   fragElements = NULL, lossElements = NULL, topMost = NULL, OM = FALSE,
                                                   negate = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertCount, . ~ minExplainedPeaks + topMost, positive = c(FALSE, TRUE),
           null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertList(scoreLimits, null.ok = TRUE, types = "numeric", add = ac)
    if (!is.null(scoreLimits))
    {
        scCols <- annScoreNames(obj, FALSE)
        checkmate::assertNames(names(scoreLimits), type = "unique", subset.of = scCols, add = ac)
        checkmate::qassertr(scoreLimits, "N2")
    }
    aapply(checkmate::assertCharacter, . ~ elements + fragElements + lossElements,
           min.chars = 1, min.len = 1, null.ok = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertFlag, . ~ OM + negate, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    cat("Filtering annotations... ")
    
    mConsNames <- mergedConsensusNames(obj)
    
    if (!is.null(minExplainedPeaks))
        scoreLimits <- modifyList(if (is.null(scoreLimits)) list() else scoreLimits,
                                  list(explainedPeaks = c(minExplainedPeaks, Inf)))
    mark <- if (negate) function(x) !x else function(x) x
    
    oldn <- length(obj)
    obj <- delete(obj, j = function(annTable, ...)
    {
        annTable <- copy(annTable)
        annTable[, keep := TRUE]
        
        if (!is.null(scoreLimits))
        {
            for (sc in names(scoreLimits))
            {
                cols <- getAllMergedConsCols(sc, names(annTable), mConsNames)
                if (length(cols) == 0)
                    next
                
                annTable[keep == TRUE, keep :=
                             mark(do.call(pmax, c(.SD, list(na.rm = TRUE))) >= scoreLimits[[sc]][1] &
                                      do.call(pmin, c(.SD, list(na.rm = TRUE))) <= scoreLimits[[sc]][2]),
                         .SDcols = cols]
            }
        }
        
        if (!is.null(elements))
            annTable[keep == TRUE, keep := sapply(neutral_formula, checkFormula, elements, negate)]
        
        if (!is.null(fragElements) || !is.null(lossElements))
        {
            annTable[keep == TRUE, keep := sapply(fragInfo, function(fi)
            {
                if (nrow(fi) == 0)
                    return(FALSE)
                if (!is.null(fragElements) && !any(sapply(fi$ion_formula, checkFormula, fragElements, negate)))
                    return(FALSE)
                if (!is.null(lossElements) && !any(sapply(fi$neutral_loss, checkFormula, lossElements, negate)))
                    return(FALSE)
                return(TRUE)
            })]
        }
        
        if (OM)
        {
            annTable <- addElementInfoToAnnTable(annTable, NULL, NULL, OM = TRUE, classify = FALSE)
            annTable[keep == TRUE, keep := mark( 
                         # rules from Kujawinski & Behn, 2006 (10.1021/ac0600306)
                         H >= 1/3 * C &
                         H <= ((2 * C) + N + 2) &
                         (H + N) %% 2 == 0 &
                         N <= C &
                         O <= C &
                         P <= 2 &
                         S <= 2 &
                         
                         # rules from Koch & dittmar 2006 (10.1002/rcm.2386)
                         sapply(DBE_AI, checkmate::checkInt) &
                         HC <= 2.2 &
                         OC <= 1.2 &
                         NC <= 0.5)]
        }
        
        return(!annTable$keep)
    })

    if (!is.null(topMost))
    {
        if (negate)
            obj <- delete(obj, j = seq_len(topMost))
        else
            obj <- delete(obj, j = function(at, ...) seq_len(nrow(at)) > topMost)
    }
    
    newn <- length(obj)
    printf("Done! Filtered %d (%.2f%%) annotations. Remaining: %d\n", oldn - newn, if (oldn == 0) 0 else (1-(newn/oldn))*100, newn)
    return(obj)
})

setMethod("plotVenn", "featureAnnotations", function(obj, ..., labels = NULL, vennArgs = NULL)
{
    allFeatAnnotations <- c(list(obj), list(...))
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertList(allFeatAnnotations, types = "featureAnnotations", min.len = 2, any.missing = FALSE,
                          unique = TRUE, .var.name = "...", add = ac)
    checkmate::assertCharacter(labels, min.chars = 1, len = length(allFeatAnnotations), null.ok = TRUE, add = ac)
    checkmate::assertList(vennArgs, names = "unique", null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (is.null(labels))
        labels <- make.unique(sapply(allFeatAnnotations, algorithm))
    if (is.null(vennArgs))
        vennArgs <- list()
    
    allAnnTabs <- lapply(allFeatAnnotations, as.data.table)
    do.call(makeVennPlot, c(list(allAnnTabs, labels, lengths(allFeatAnnotations), function(obj1, obj2)
    {
        if (length(obj1) == 0 || length(obj2) == 0)
            return(data.table())
        fintersect(obj1[, c("group", "UID")], obj2[, c("group", "UID")])
    }, nrow), vennArgs))
})

setMethod("plotUpSet", "featureAnnotations", function(obj, ..., labels = NULL, nsets = length(list(...)) + 1,
                                                      nintersects = NA, upsetArgs = NULL)
{
    allFeatAnnotations <- c(list(obj), list(...))
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertList(allFeatAnnotations, types = "featureAnnotations", min.len = 2, any.missing = FALSE,
                          unique = TRUE, .var.name = "...", add = ac)
    checkmate::assertCharacter(labels, min.chars = 1, len = length(allFeatAnnotations), null.ok = TRUE, add = ac)
    checkmate::assertList(upsetArgs, names = "unique", null.ok = TRUE, add = ac)
    checkmate::assertCount(nsets, positive = TRUE)
    checkmate::assertCount(nintersects, positive = TRUE, na.ok = TRUE)
    checkmate::reportAssertions(ac)
    
    if (is.null(labels))
        labels <- make.unique(sapply(allFeatAnnotations, algorithm))
    
    allAnnTabs <- mapply(allFeatAnnotations, labels, SIMPLIFY = FALSE, FUN = function(f, l)
    {
        ret <- as.data.table(f)
        if (length(ret) == 0)
            ret <- data.table(group = character(), InChIKey1 = character())
        ret <- unique(ret[, c("group", "UID")])[, (l) := 1]
    })
    
    annTab <- Reduce(function(f1, f2)
    {
        merge(f1, f2, by = c("group", "UID"), all = TRUE)
    }, allAnnTabs)
    
    annTab <- annTab[, labels, with = FALSE]
    for (j in seq_along(annTab))
        set(annTab, which(is.na(annTab[[j]])), j, 0)
    
    if (sum(sapply(annTab, function(x) any(x>0))) < 2)
        stop("Need at least two non-empty objects to plot")
    
    do.call(UpSetR::upset, c(list(annTab, nsets = nsets, nintersects = nintersects), upsetArgs))
})

setMethod("prepareConsensusLabels", "featureAnnotations", function(obj, ..., labels)
{
    if (is.null(labels))
        labels <- sapply(list(obj, ...), algorithm)
    
    # in case names are (still) duplicated
    labels <- make.unique(labels)
    
    return(labels)
