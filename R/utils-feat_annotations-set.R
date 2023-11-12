#' @include main.R
#' @include formulas.R
#' @include compounds.R
NULL

featAnnSetSpecificScoreCols <- function()
{
    return(c("score", "fragScore", "isoScore", "metFusionScore", "individualMoNAScore", "peakFingerprintScore",
             "lossFingerprintScore", "libMatch", "formulaScore", "RF_SIRFP", "LC50_SIRFP", "combMatch", "isoScore",
             "mSigma", "MSMSScore"))
}

makeFeatAnnSetConsensus <- function(setObjects, origFGNames, setThreshold, setThresholdAnn, setAvgSpecificScores,
                                    mConsNames)
{
    # generate consensus by...
    # - checking setThreshold/setThresholdAnn
    # - merging by UID
    # - average scores
    
    
    sumMergedScoreRows <- function(sc, m) .rowSums(unlist(sc), na.rm = TRUE, m = m, n = 2)
    
    # used for merging below
    otherCols <- function(col) if (length(col) > 0) paste0("i.", col) else character()
    
    # get all annotated fGroup names with original order
    allAnnGNames <- intersect(origFGNames, unique(unlist(lapply(setObjects, groupNames))))
    
    sCount <- length(setObjects)
    
    allScoreCols <- unique(unlist(lapply(setObjects, annScoreNames, TRUE)))
    
    cons <- sapply(allAnnGNames, function(gName)
    {
        allResults <- pruneList(sapply(setObjects, "[[", gName, simplify = FALSE))
        
        # init some columns before merging
        allResults <- Map(allResults, names(allResults), f = function(ct, s)
        {
            ct <- copy(ct)
            rname <- paste0("rank-", s)
            ct[, c("set", "setsCount", "setsMergedCount", rname) := .(s, 1, 1, seq_len(.N))]
            
            # rename cols that are specific to a set or algo consensus or should otherwise not be combined
            cols <- getAllMergedConsCols(c(
                "rank", "mergedBy", "coverage", "explainedPeaks", "ion_formula", "ion_formula_mz", "precursorType",
                "libPeaksCompared", "libPeaksTotal", "annSim", "estIDLevel"
            ), names(ct), mConsNames)
            if (!setAvgSpecificScores)
                cols <- c(cols, getAllMergedConsCols(featAnnSetSpecificScoreCols(), names(ct), mConsNames))
            if (length(cols) > 0)
                setnames(ct, cols, paste0(cols, "-", s))
            
            return(ct)
        })
        
        if (length(allResults) == 1)
            return(allResults[[1]])
        
        Reduce(x = allResults, f = function(left, right)
        {
            scoreColsLeft <- intersect(allScoreCols, names(left))
            scoreColsRight <- intersect(allScoreCols, names(right))
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
                       fifelse(!is.na(mget(otherColsBoth)), mget(otherColsBoth), mget(otherCols(otherColsBoth))),
                       # add missing columns from right (if any)
                       mget(otherCols(colsOnlyRight)),
                       # mark set presence
                       .(paste0(set, ",", i.set)),
                       # combine fragInfos
                       .(Map(fragInfo, i.fragInfo, f = rbind, MoreArgs = list(fill = TRUE)))),
                 on = "UID"]
            
            left[UID %chin% right$UID, setsMergedCount := setsMergedCount + 1]
            
            # add missing candidates from right
            left <- rbind(left, right[!UID %chin% left$UID], fill = TRUE)
            
            left[, setsCount := setsCount + 1]
            
            return(left)
        })
    }, simplify = FALSE)
    
    # update average scores and convert absolute merge counts to coverage
    cons <- lapply(cons, function(ct)
    {
        scCols <- intersect(allScoreCols, names(ct))
        if (length(scCols) > 0)
            ct[, (scCols) := lapply(.SD, function(x) x / setsMergedCount), .SDcols = scCols]
        
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
        
        ct[, c("setCoverageAnn", "setCoverage") := .(setsMergedCount / setsCount, setsMergedCount / sCount)]
        ct[, c("setsCount", "setsMergedCount") := NULL]
        return(ct)
    })
    
    if (setThreshold > 0 || setThresholdAnn > 0)
        cons <- pruneList(lapply(cons,
                                 function(ct) ct[setCoverage >= setThreshold & setCoverageAnn >= setThresholdAnn]),
                          checkZeroRows = TRUE)
    
    return(cons)
}

makeAnnSetScorings <- function(setObjects, setAvgSpecificScores, origFGNames)
{
    scTypes <- character()
    scRanges <- list()
    
    if (length(setObjects) > 0)
    {
        allScTypes <- sapply(setObjects, slot, "scoreTypes", simplify = FALSE)
        allScRanges <- sapply(setObjects, slot, "scoreRanges", simplify = FALSE)
        
        if (!setAvgSpecificScores)
        {
            # rename set specific scorings
            renameSc <- function(v, sn)
            {
                wh <- grepl(paste0(featAnnSetSpecificScoreCols(), collapse = "|"), v)
                v[wh] <- paste0(v[wh], "-", sn)
                return(v)
            }
            allScTypes <- Map(allScTypes, names(allScTypes), f = renameSc)
            allScRanges <- Map(allScRanges, names(allScRanges), f = function(sr, sn)
            {
                lapply(sr, function(r) setNames(r, renameSc(names(r), sn)))
            })
        }
        
        scTypes <- unique(unlist(allScTypes))
            
        scRanges <- Reduce(x = allScRanges, f = function(left, right)
        {
            # change ranges for overlap
            groupsLR <- intersect(names(left), names(right))
            ret <- Map(left[groupsLR], right[groupsLR], f = function(rangesL, rangesR)
            {
                # overlap
                scLR <- intersect(names(rangesL), names(rangesR))
                rret <- Map(range, rangesL[scLR], rangesR[scLR])
                
                # unique left
                scL <- setdiff(names(rangesL), names(rangesR))
                rret[scL] <- rangesL[scL]
                
                # unique right
                scR <- setdiff(names(rangesR), names(rangesL))
                rret[scR] <- rangesR[scR]
                
                return(rret)
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

doUpdateSetConsensus <- function(obj)
{
    if (length(setObjects(obj)) >= 1)
    {
        obj@groupAnnotations <- makeFeatAnnSetConsensus(obj@setObjects, obj@origFGNames,
                                                        obj@setThreshold, obj@setThresholdAnn,
                                                        obj@setAvgSpecificScores,
                                                        mergedConsensusNames(obj, FALSE))
    }
    else
        obj@groupAnnotations <- list()
    
    obj@scoreRanges <- obj@scoreRanges[groupNames(obj)]
    
    return(obj)
}

initSetFragInfos <- function(setObjects, MSPeakListsSet)
{
    # update fragInfos
    for (s in names(setObjects))
    {
        for (fg in groupNames(setObjects[[s]]))
        {
            setObjects[[s]]@groupAnnotations[[fg]][, fragInfo := lapply(fragInfo, function(fi)
            {
                fi <- copy(fi) # BUG: avoid warning that somehow was incorrectly copied (invalid .internal.selfref)
                fi[, set := s]
                return(fi)
            })]
        }
    }
    
    return(setObjects)
}

doFeatAnnConsensusSets <- function(allAnnObjs, MSPeakLists, specSimParams, labels, setThreshold, setThresholdAnn,
                                   setAvgSpecificScores, rankWeights)
{
    # make consensus of shared setObjects
    # add unique setObjects
    # make 'regular' set consensus from new results
    
    if (!allSame(lapply(allAnnObjs, sets)))
        stop("All objects must have the same sets.")
    
    allOrigFGNames <- pruneList(lapply(allAnnObjs, slot, "origFGNames"), checkEmptyElements = TRUE)
    if (!allSame(allOrigFGNames))
        stop("All objects must have been generated from the same feature groups.")
    origFGNames <- if (length(allOrigFGNames) == 0) character() else allAnnObjs[[1]]@origFGNames

    # NOTE: filtering (thresholds, unique) is not performed here: this is done afterwards in the consensus methods, as
    # it makes more sense to filter the end result instead of those from set objects
    consArgs <- list(specSimParams = specSimParams, rankWeights = rankWeights, labels = labels, absMinAbundance = NULL,
                     relMinAbundance = NULL, uniqueFrom = NULL)
    unsetMSPeakLists <- checkAndUnSetOther(sets(allAnnObjs[[1]]), MSPeakLists, "MSPeakLists")
    setObjects <- sapply(sets(allAnnObjs[[1]]), function(set)
    {
        return(do.call(consensus, c(lapply(lapply(allAnnObjs, setObjects), "[[", set),
                                    list(MSPeakLists = unsetMSPeakLists[[set]]), consArgs)))
    }, simplify = FALSE)

    cons <- makeFeatAnnSetConsensus(setObjects, origFGNames, setThreshold, setThresholdAnn, setAvgSpecificScores,
                                    labels)
    cons <- lapply(cons, function(at)
    {
        at <- copy(at)
        cols <- getAllMergedConsCols("mergedBy", names(at), names(setObjects))
        at[, mergedBy := {
            # collapse all cols with comma, split by comma, take unique, re-collapse with comma
            allMB <- unlist(mget(cols))
            allMB <- allMB[!is.na(allMB)]
            allMB <- unique(unlist(strsplit(allMB, ",")))
            paste0(allMB, collapse = ",")
        }, by = seq_len(nrow(at))]
        at[, coverage := sapply(mergedBy, function(mb) (countCharInStr(mb, ",") + 1) / length(allAnnObjs))]
        return(at)
    })
    
    return(list(setObjects = setObjects, groupAnnotations = cons,
                algorithm = paste0(unique(sapply(allAnnObjs, algorithm)), collapse = ","),
                origFGNames = origFGNames, mergedConsensusNames = labels))
}

doFeatAnnUnset <- function(obj, set)
{
    hash <- makeHash(obj, set)
    cd <- loadCacheData("featAnnUnset", hash)
    if (!is.null(cd))
        return(cd)
    
    obj <- obj[, sets = set]
    
    ann <- lapply(annotations(obj), copy)
    
    # get rid of sets specific columns
    ann <- lapply(ann, data.table::set, j = c("set", "setCoverage", "setCoverageAnn"), value = NULL)
    
    # ... and in fragInfos
    ann <- lapply(ann, function(tbl)
    {
        tbl[, fragInfo := lapply(fragInfo, function(fi)
        {
            fi <- copy(fi)
            set(fi, j = "set", value = NULL)
        })]
    })
    
    pat <- paste0("\\-", set, "$")
    
    # restore sets specific columns
    ann <- lapply(ann, function(a)
    {
        cols <- grep(pat, names(a), value = TRUE)
        setnames(a, cols, sub(pat, "", cols))
        a[, rank := NULL] # no need for this anymore
        return(a)
    })

    # restore scorings
    scTypes <- sub(pat, "", obj@scoreTypes)
    scRanges <- lapply(obj@scoreRanges, function(sc) setNames(sc, gsub(pat, "", names(sc))))
        
    ret <- list(annotations = ann, scoreTypes = scTypes, scoreRanges = scRanges)
    saveCacheData("featAnnUnset", ret, hash)
    
    return(ret)
}

doFeatAnnDeleteSets <- function(obj, i, j, ...)
{
    old <- obj
    obj <- callNextMethod()
    
    # sync setObjects
    annTab <- annotations(obj); annTabOld <- annotations(old)
    obj@setObjects <- lapply(obj@setObjects, function(so)
    {
        delete(so, j = function(atso, grp)
        {
            if (is.null(annTab[[grp]]))
                return(TRUE) # fully removed group
            
            # remove removed...
            return(atso$UID %chin% setdiff(annTabOld[[grp]]$UID, annTab[[grp]]$UID))
        })
    })
    
    # update ranks
    obj@groupAnnotations <- Map(obj@groupAnnotations, groupNames(obj), f = function(at, grp)
    {
        rankCols <- getAllMergedConsCols("rank", names(at), mergedConsensusNames(obj))
        at[, (rankCols) := lapply(rankCols, function(rc)
        {
            s <- sub("^rank\\-", "", rc)
            atso <- setObjects(obj)[[s]][[grp]]
            if (is.null(atso))
                return(NA_integer_)
            return(match(UID, atso$UID, nomatch = NA_integer_))
        })][]
    })
    
    return(obj)
}

doFeatAnnSubsetSets <- function(x, i, j, ..., sets = NULL, updateConsensus = FALSE, drop = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    assertSets(x, sets, TRUE, add = ac)
    checkmate::assertFlag(updateConsensus, add = ac)
    checkmate::reportAssertions(ac)
    
    if (!is.null(sets))
    {
        oldSets <- sets(x)
        x@setObjects <- x@setObjects[sets]
        if (!updateConsensus) # update sets result; otherwise done by updateSetConsensus() when new consensus is made
        {
            rmSets <- setdiff(oldSets, sets(x))
            if (length(rmSets) > 0)
            {
                pat <- paste0("\\-(", paste0(rmSets, collapse = "|"), ")$")
                x@groupAnnotations <- lapply(x@groupAnnotations, function(ct)
                {
                    cols <- grep(pat, names(ct), value = TRUE)
                    if (length(cols) == 0)
                        return(ct) # set already not present
                    
                    ct <- copy(ct)
                    
                    # remove set specific columns
                    ct[, (cols) := NULL]
                    
                    # update set column
                    rankCols <- intersect(paste0("rank-", sets(x)), names(x))
                    ct[, set := {
                        s <- unlist(strsplit(set, ","))
                        paste0(setdiff(s, rmSets), collapse = ",")
                    }, by = seq_len(nrow(ct))]
                    
                    # update fragInfos
                    ct[, fragInfo := lapply(fragInfo, function(fi) fi[!set %chin% rmSets])]
                    
                    return(ct[])       
                })
                
                # remove results from removed sets --> those are now without set assignment
                x <- delete(x, j = function(at, ...) !nzchar(at$set))
                
                # update scorings
                x@scoreTypes <- x@scoreTypes[!grepl(pat, x@scoreTypes)]
                x@scoreRanges <- lapply(x@scoreRanges, function(sr)
                {
                    sr <- sr[!grepl(pat, names(sr))]
                    return(sr)
                })
            }
        }
    }
    
    if (!missing(i))
    {
        if (updateConsensus)
        {
            # NOTE: assume that subsetting with non-existing i will not result in errors
            i <- assertSubsetArgAndToChr(i, groupNames(x))
            x@setObjects <- lapply(x@setObjects, "[", i = i)
        }
        else
            x <- callNextMethod()
        
    }
    
    if ((!is.null(sets) || !missing(i)) && updateConsensus)
        x <- updateSetConsensus(x)
    
    return(x)
}

doFeatAnnFilterSets <- function(obj, ..., sets = NULL, updateConsensus = FALSE, negate = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    assertSets(obj, sets, TRUE, add = ac)
    checkmate::assertFlag(negate, add = ac)
    checkmate::assertFlag(updateConsensus, add = ac)
    checkmate::reportAssertions(ac)
    
    if (!is.null(sets) && length(sets) > 0)
    {
        if (negate)
            sets <- setdiff(get("sets", pos = 2)(obj), sets)
        obj <- obj[, sets = sets, updateConsensus = updateConsensus]
    }
    
    if (...length() > 0)
    {
        if (updateConsensus)
        {
            # filter set objects and re-generate annotation consensus
            
            obj@setObjects <- lapply(obj@setObjects, filter, ..., negate = negate)
            
            # synchronize other objects
            cat("Synchronizing set objects...\n")
            obj <- updateSetConsensus(obj)
            cat("Done!\n")
        }
        else
            obj <- callNextMethod(obj, ..., negate = negate)
    }
    
    return(obj)
}

doFeatAnnMCNSets <- function(obj, sets) return(if (sets) patRoon:::sets(obj) else character())

doFeatAnnMCNSetsCons <- function(obj, sets)
{
    if (sets)
        return(c(sets(obj), obj@mergedConsensusNames, sapply(obj@mergedConsensusNames, paste0, "-", sets(obj))))
    return(obj@mergedConsensusNames)
}

doFeatAnnPredictRFSets <- function(obj, fGroups, calibrants, ...)
{
    checkmate::assertClass(fGroups, "featureGroupsSet")
    checkmate::assertList(calibrants, any.missing = FALSE, types = "data.frame", len = length(sets(obj)))
    
    unsetFGs <- checkAndUnSetOther(sets(obj), fGroups, "fGroups")
    obj@setObjects <- Map(setObjects(obj), unsetFGs, calibrants, f = predictRespFactors, MoreArgs = list(...))
    obj <- updateSetConsensus(obj)
    
    return(obj)
}

doFeatAnnPredictToxSets <- function(obj, ...)
{
    obj@setObjects <- lapply(setObjects(obj), predictTox, ...)
    return(updateSetConsensus(obj))
}
