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
        x@groupAnnotations <- x@groupAnnotations[i]
        x@scoreRanges <- x@scoreRanges[i]
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
    
    mcn <- mergedConsensusNames(x)
    annTable <- annotations(x)
    if (normalizeScores != "none")
    {
        annTable <- Map(annTable, x@scoreRanges, f = normalizeAnnScores,
                        MoreArgs = list(scoreCols = annScoreNames(obj), mConsNames = mcn, normalizeScores == "minmax",
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

setMethod("filter", "featureAnnotations", function(obj, minExplainedPeaks = NULL, scoreLimits = NULL, elements = NULL,
                                                   fragElements = NULL, lossElements = NULL, topMost = NULL,
                                                   negate = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertCount, . ~ minExplainedPeaks + topMost, positive = c(FALSE, TRUE),
           null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertList(scoreLimits, null.ok = TRUE, types = "numeric", add = ac)
    if (!is.null(scoreLimits))
    {
        scCols <- annScoreNames(obj)
        checkmate::assertNames(names(scoreLimits), type = "unique", subset.of = scCols, add = ac)
        checkmate::qassertr(scoreLimits, "N2")
    }
    aapply(checkmate::assertCharacter, . ~ elements + fragElements + lossElements,
           min.chars = 1, min.len = 1, null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertFlag(negate, add = ac)
    checkmate::reportAssertions(ac)
    
    cat("Filtering annotations... ")
    
    mConsNames <- mergedConsensusNames(obj)
    
    # UNDONE: minExplainedPeaks (include in scoreLimits?)
    
    oldn <- length(obj)
    obj@groupAnnotations <- sapply(obj@groupAnnotations, function(annTable)
    {
        if (!is.null(scoreLimits))
        {
            for (sc in names(scoreLimits))
            {
                cols <- getAllMergedConsCols(sc, names(annTable), mConsNames)
                if (length(cols) == 0)
                    next
                
                keep <- annTable[, do.call(pmin, c(.SD, list(na.rm = TRUE))) >= scoreLimits[[sc]][1] &
                                     do.call(pmax, c(.SD, list(na.rm = TRUE))) <= scoreLimits[[sc]][2],
                                 .SDcols = cols]
                if (negate)
                    keep <- !keep
                annTable <- annTable[keep]
            }
        }
        
        if (nrow(annTable) == 0)
            return(annTable)
        
        if (!is.null(elements))
        {
            keep <- sapply(annTable$formula, checkFormula, elements, negate = negate)
            annTable <- annTable[keep]
        }
        if (!is.null(fragElements) || !is.null(lossElements))
        {
            keep <- sapply(annTable$fragInfo, function(fi)
            {
                if (nrow(fi) == 0)
                    return(FALSE)
                if (!is.null(fragElements) && !any(sapply(fi$formula, checkFormula, fragElements, negate = negate)))
                    return(FALSE)
                if (!is.null(lossElements) && !any(sapply(fi$neutral_loss, checkFormula, lossElements, negate = negate)))
                    return(FALSE)
                return(TRUE)
            })
            annTable <- annTable[keep]
        }
        
        if (!is.null(topMost))
        {
            if (negate)
                annTable <- tail(annTable, topMost)
            else
                annTable <- head(annTable, topMost)
        }
        
        return(annTable)
    }, simplify = FALSE)
    
    if (length(obj) > 0)
        obj@groupAnnotations <- obj@groupAnnotations[sapply(obj@groupAnnotations, function(cm) !is.null(cm) && nrow(cm) > 0)]
    
    obj@scoreRanges <- obj@scoreRanges[names(obj@groupAnnotations)]
    
    newn <- length(obj)
    printf("Done! Filtered %d (%.2f%%) annotations. Remaining: %d\n", oldn - newn, if (oldn == 0) 0 else (1-(newn/oldn))*100, newn)
    return(obj)
})

setMethod("plotScores", "featureAnnotations", function(obj, index, groupName, normalizeScores = "max",
                                                       excludeNormScores = defaultExclNormScores(x),
                                                       onlyUsed = TRUE, useGGPlot2 = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertCount(index, positive = TRUE, add = ac)
    checkmate::assertString(groupName, min.chars = 1, add = ac)
    checkmate::assertChoice(normalizeScores, c("none", "max", "minmax"))
    checkmate::assertCharacter(excludeNormScores, min.chars = 1, null.ok = TRUE, add = ac)
    checkmate::assertFlag(onlyUsed, add = ac)
    checkmate::assertFlag(useGGPlot2, add = ac)
    checkmate::reportAssertions(ac)
    
    annTable <- annotations(obj)[[groupName]]
    
    if (is.null(annTable) || nrow(annTable) == 0 || index > nrow(annTable))
        return(NULL)
    
    mcn <- mergedConsensusNames(obj)
    
    if (normalizeScores != "none")
        annTable <- normalizeAnnScores(annTable, annScoreNames(obj), obj@scoreRanges[[groupName]], mcn,
                                       normalizeScores == "minmax", excludeNormScores)
    
    scoreCols <- getAllMergedConsCols(c(getCompScoreColNames(), getCompSuspectListColNames()), names(annTable), mcn)
    if (onlyUsed)
        scoreCols <- intersect(scoreCols, obj@scoreTypes)
    
    makeScoresPlot(annTable[index, scoreCols, with = FALSE], mcn, useGGPlot2)
})

setMethod("plotScoresHash", "featureAnnotations", function(obj, index, groupName, normalizeScores = "max",
                                                           excludeNormScores = defaultExclNormScores(x),
                                                           onlyUsed = TRUE, useGGPlot2 = FALSE)
{
    annTable <- annotations(obj)[[groupName]]
    if (is.null(annTable) || nrow(annTable) == 0 || index > nrow(annTable))
        annTable <- NULL
    else if (normalizeScores == "none")
        annTable <- annTable[index]
    
    return(makeHash(index, annTable, normalizeScores, excludeNormScores, onlyUsed, useGGPlot2))
})

setMethod("annotatedPeakList", "featureAnnotations", function(obj, index, groupName, MSPeakLists, onlyAnnotated = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertCount(index, positive = TRUE, add = ac)
    assertChoiceSilent(groupName, groupNames(obj), add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertFlag(onlyAnnotated, add = ac)
    checkmate::reportAssertions(ac)
    
    spec <- MSPeakLists[[groupName]][["MSMS"]]
    if (is.null(spec))
        return(NULL)
    
    spec <- copy(spec)
    spec[, PLIndex := seq_len(nrow(spec))] # for merging
    
    annTable <- annotations(obj)[[groupName]]
    fragInfo <- NULL
    if (!is.null(annTable) && nrow(annTable) >= index)
    {
        compr <- annTable[index]
        fragInfo <- compr$fragInfo[[1]]
        
        if (!is.null(fragInfo))
            spec <- merge(spec, fragInfo[, -c("intensity", "mz")], all.x = TRUE, by = "PLIndex")
        
        spec <- spec[, PLIndex := NULL]
    }
    
    if (onlyAnnotated)
    {
        if (is.null(spec[["formula"]]))
            spec <- spec[0]
        else
            spec <- spec[!is.na(formula)]
    }
    
    return(spec[])
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
    
    allCompoundTabs <- lapply(allFeatAnnotations, as.data.table)
    do.call(makeVennPlot, c(list(allCompoundTabs, labels, lengths(allFeatAnnotations), function(obj1, obj2)
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
    
    allCompsTabs <- mapply(allFeatAnnotations, labels, SIMPLIFY = FALSE, FUN = function(f, l)
    {
        ret <- as.data.table(f)
        if (length(ret) == 0)
            ret <- data.table(group = character(), InChIKey1 = character())
        ret <- unique(ret[, c("group", "UID")])[, (l) := 1]
    })
    
    compTab <- Reduce(function(f1, f2)
    {
        merge(f1, f2, by = c("group", "UID"), all = TRUE)
    }, allCompsTabs)
    
    compTab <- compTab[, labels, with = FALSE]
    for (j in seq_along(compTab))
        set(compTab, which(is.na(compTab[[j]])), j, 0)
    
    if (sum(sapply(compTab, function(x) any(x>0))) < 2)
        stop("Need at least two non-empty objects to plot")
    
    do.call(UpSetR::upset, c(list(compTab, nsets = nsets, nintersects = nintersects), upsetArgs))
})

