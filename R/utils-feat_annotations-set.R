#' @include main.R
#' @include formulas.R
#' @include compounds.R
NULL

makeFeatAnnSetConsensus <- function(setObjects, origFGNames, setThreshold, setThresholdAnn, mConsNames)
{
    # generate consensus by...
    # - checking setThreshold/setThresholdAnn
    # - merging by UID
    # - average scores
    # - merge fragInfos and update PLIndex
    
    
    sumMergedScoreRows <- function(sc, m) .rowSums(unlist(sc), na.rm = TRUE, m = m, n = 2)
    
    # used for merging below
    otherCols <- function(col) if (length(col) > 0) paste0("i.", col) else character()
    
    # get all annotated fGroup names with original order
    allAnnGNames <- intersect(origFGNames, unique(unlist(lapply(setObjects, groupNames))))
    
    sCount <- length(setObjects)
    
    scoreCols <- if (sCount > 0) annScoreNames(setObjects[[1]], TRUE)
    
    cons <- sapply(allAnnGNames, function(gName)
    {
        allResults <- pruneList(sapply(setObjects, "[[", gName, simplify = FALSE))
        
        # init some columns before merging
        allResults <- Map(allResults, names(allResults), f = function(ct, s)
        {
            ct <- copy(ct)
            rname <- paste0("rank-", s)
            ct[, c("set", "setCoverage", "setCoverageAnn", rname) := .(s, 1, 1, seq_len(.N))]
            
            # rename cols that are specific to a set or algo consensus or should otherwise not be combined
            cols <- getAllMergedConsCols(c("rank", "mergedBy", "coverage", "explainedPeaks"), names(ct), mConsNames)
            if (length(cols) > 0)
                setnames(ct, cols, paste0(cols, "-", s))
            
            return(ct)
        })
        
        if (length(allResults) == 1)
            return(allResults[[1]])
        
        Reduce(x = allResults, f = function(left, right)
        {
            scoreColsLeft <- getAllMergedConsCols(scoreCols, names(left), mConsNames)
            scoreColsRight <- getAllMergedConsCols(scoreCols, names(right), mConsNames)
            scoreColsBoth <- intersect(scoreColsLeft, scoreColsRight)
            # UNDONE: below cols are specific to compounds --> make method?
            combineCols <- c("compoundName", "compoundName2", "identifier", "relatedCIDs")
            combineColsBoth <- intersect(getAllMergedConsCols(combineCols, names(left), mConsNames),
                                         getAllMergedConsCols(combineCols, names(right), mConsNames))
            otherColsBoth <- setdiff(intersect(names(left), names(right)),
                                     c(scoreColsBoth, combineColsBoth, "set", "fragInfo"))
            colsOnlyRight <- setdiff(names(right), names(left))
            
            left[right, (c(scoreColsBoth, combineColsBoth, otherColsBoth, colsOnlyRight, "set", "fragInfo")) :=
                     # sum scores (they are averaged later below)
                     c(lapply(scoreColsBoth, function(col) sumMergedScoreRows(mget(c(col, otherCols(col)),
                                                                                   inherits = TRUE), .N)),
                       # combine some other overlapping columns
                       lapply(combineColsBoth, function(col)
                       {
                           vals <- mget(c(col, otherCols(col)), inherits = TRUE)
                           vals <- lapply(vals, strsplit, ";")
                           return(mapply(vals[[1]], vals[[2]], FUN = function(l, r)
                           {
                               l <- l[!is.na(l)]; r <- r[!is.na(r)]
                               if (length(l) == 0 && length(r) == 0)
                                   return(NA)
                               return(paste0(union(l, r), collapse = ";"))
                           }))
                       }),
                       # handle all other overlapping columns: select non NA
                       fifelse(is.na(mget(otherColsBoth)), mget(otherColsBoth), mget(otherCols(otherColsBoth))),
                       # add missing columns from right (if any)
                       mget(otherCols(colsOnlyRight)),
                       # mark set presence
                       .(paste0(set, ",", i.set)),
                       # combine fragInfos
                       .(Map(fragInfo, i.fragInfo, f = rbind, MoreArgs = list(fill = TRUE)))),
                 on = "UID"]
            
            # add missing candidates from right
            left <- rbind(left, right[!UID %chin% left$UID], fill = TRUE)
            
            left[, setCoverage := setCoverage + 1]
            left[UID %chin% right$UID, setCoverageAnn := setCoverageAnn + 1]
            
            return(left)
        })
    }, simplify = FALSE)
    
    # update fragInfos, average scores and convert absolute merge counts to coverage
    cons <- lapply(cons, function(ct)
    {
        ct[, fragInfo := lapply(fragInfo, function(fi)
        {
            fi <- copy(fi) # avoid DT warning/bug
            fi[, c("PLIndex", "PLIndexSet", "PLIndexOrig") := .(PLIndexSet, NULL, PLIndex)]
        })]
        
        scCols <- getAllMergedConsCols(scoreCols, names(ct), mConsNames)
        ct[, (scCols) := lapply(.SD, function(x) x / setCoverageAnn), .SDcols = scCols]
        
        # re-sort by avg rank scores
        rnames <- getAllMergedConsCols("rank", names(ct), names(setObjects))
        ncand <- nrow(ct)
        ct[, rankScore := {
            invRanks <- (ncand - (unlist(.SD) - 1)) / ncand
            invRanks[is.na(invRanks)] <- 0
            mean(invRanks)
        }, .SDcols = rnames, by = seq_len(ncand)]
        
        setorderv(ct, "rankScore", -1)
        ct[, rankScore := NULL]
        
        ct[, c("setCoverageAnn", "setCoverage") := .(setCoverageAnn / setCoverage, setCoverage / sCount)]
        return(ct)
    })
    
    if (setThreshold > 0 || setThresholdAnn > 0)
        cons <- pruneList(lapply(cons,
                                 function(ct) ct[setCoverage >= setThreshold & setCoverageAnn >= setThresholdAnn]),
                          checkZeroRows = TRUE)
    
    return(cons)
}

makeAnnSetScorings <- function(setObjects, origFGNames)
{
    scTypes <- character()
    if (length(setObjects) > 0)
        scTypes <- unique(unlist(lapply(setObjects, slot, "scoreTypes")))
    
    scRanges <- list()
    if (length(setObjects) > 0)
    {
        scRanges <- Reduce(x = lapply(setObjects, slot, "scoreRanges"), f = function(left, right)
        {
            # change ranges for overlap
            groupsLR <- intersect(names(left), names(right))
            ret <- mapply(left[groupsLR], right[groupsLR], SIMPLIFY = FALSE, FUN = function(rangesL, rangesR)
            {
                scLR <- names(rangesL) # should be same for left/right
                mapply(rangesL[scLR], rangesR[scLR], FUN = range, SIMPLIFY = FALSE)
            })
            
            # add unique from left
            groupsOnlyL <- setdiff(names(left), names(right))
            ret[groupsOnlyL] <- left[groupsOnlyL]
            
            # add unique from right
            groupsOnlyR <- setdiff(names(right), names(left))
            ret[groupsOnlyR] <- right[groupsOnlyR]
            
            return(ret[intersect(origFGNames, names(ret))]) # order
        })
    }
    
    return(list(scTypes = scTypes, scRanges = scRanges))
}

syncAnnSetObjects <- function(obj, makeCons)
{
    if (length(setObjects(obj)) >= 1)
    {
        if (makeCons)
            obj@groupAnnotations <- makeFeatAnnSetConsensus(obj@setObjects, obj@origFGNames,
                                                                     obj@setThreshold, obj@setThresholdAnn,
                                                                     mergedConsensusNames(obj))
        else
        {
            # sync available feature groups
            allFGroups <- unique(unlist(lapply(setObjects(obj), groupNames)))
            obj@groupAnnotations <- obj@groupAnnotations[intersect(groupNames(obj), allFGroups)]
            
            # only keep results from sets still present
            spat <- paste0(sets(obj), collapse = "|")
            obj@groupAnnotations <- lapply(obj@groupAnnotations, function(ct) ct[grepl(spat, set)])
        }
    }
    else
        obj@groupAnnotations <- list()
    
    obj@scoreRanges <- obj@scoreRanges[groupNames(obj)]
    
    # update scoreTypes/scoreRanges
    sc <- makeAnnSetScorings(setObjects(obj), obj@origFGNames)
    obj@scoreTypes <- sc$scTypes
    obj@scoreRanges <- sc$scRanges
    
    return(obj)
}

doFeatAnnConsensusSets <- function(allAnnObjs, origFGNames, labels, setThreshold, setThresholdAnn, ...)
{
    # make consensus of shared setObjects
    # add unique setObjects
    # make 'regular' set consensus from new results
    
    if (!allSame(lapply(allAnnObjs, sets)))
        stop("All objects must have the same sets.")
    if (!allSame(lapply(allAnnObjs, slot, "origFGNames")))
        stop("All objects must have been generated from the same feature groups.")
    
    # NOTE: don't want to keep -set suffix
    labels <- if (!is.null(labels)) labels else sub("\\-set$", "", sapply(allAnnObjs, algorithm))
    
    setObjects <- sapply(sets(allAnnObjs[[1]]), function(set)
    {
        return(do.call(consensus, c(lapply(lapply(allAnnObjs, setObjects), "[[", set), list(..., labels = labels))))
    }, simplify = FALSE)
    
    cons <- makeFeatAnnSetConsensus(setObjects, origFGNames, setThreshold, setThresholdAnn, labels)
    sc <- makeAnnSetScorings(setObjects, origFGNames)
    
    return(list(setObjects = setObjects, groupAnnotations = cons, scoreTypes = sc$scTypes, scoreRanges = sc$scRanges,
                algorithm = paste0(unique(sapply(allAnnObjs, algorithm)), collapse = ","),
                mergedConsensusNames = labels))
}

# HACK: formulas/compounds share a lot of methods, but there is no clean and proper way to do this with multiple
# inheritance. Instead, simply automatically define the same method for both

setMethodMult("[", list(c("formulasSet", "ANY", "missing", "missing"), c("compoundsSet", "ANY", "missing", "missing")),
              function(x, i, j, ..., sets = NULL, drop = TRUE)
{
    assertSets(x, sets, TRUE)
    
    if (!is.null(sets))
        x@setObjects <- x@setObjects[sets]
    
    if (!missing(i))
    {
        # NOTE: assume that subsetting with non-existing i will not result in errors
        i <- assertSubsetArgAndToChr(i, groupNames(x))
        x@setObjects <- lapply(x@setObjects, "[", i = i)
    }
    
    if (!is.null(sets) || !missing(i))
        x <- syncAnnSetObjects(x, FALSE)
    
    return(x)
})

setMethodMult("annotatedPeakList", c("formulasSet", "compoundsSet"), function(obj, ...)
{
    ret <- callNextMethod()
    if (!is.null(ret[["set.x"]]))
    {
        # remove duplicate set column resulting from merging peak lists and fragInfo
        setnames(ret, "set.x", "set")
        ret[, set.y := NULL]
    }
    if (!is.null(ret[["PLIndexOrig"]]))
        ret[, PLIndexOrig := NULL]
    
    return(ret[])
})

setMethodMult("filter", c("formulasSet", "compoundsSet"), function(obj, ..., negate = FALSE, sets = NULL)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFlag(negate, add = ac)
    assertSets(obj, sets, TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (!is.null(sets) && length(sets) > 0)
    {
        if (negate)
            sets <- setdiff(get("sets", pos = 2)(obj), sets)
        obj <- obj[, sets = sets]
    }
    
    if (...length() > 0)
    {
        # filter set objects and re-generate annotation consensus
        
        obj@setObjects <- lapply(obj@setObjects, filter, ..., negate = negate)
        obj@setObjects <- pruneList(obj@setObjects, checkEmptyElements = TRUE)
        
        # synchronize other objects
        cat("Synchronizing set objects...\n")
        obj <- syncAnnSetObjects(obj, TRUE)
        cat("Done!\n")
    }
    
    return(obj)
})
