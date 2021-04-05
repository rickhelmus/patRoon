normalizeAnnScores <- function(res, scoreCols, scoreRanges, mConsNames, minMaxNormalization, exclude = NULL)
{
    res <- copy(res)
    columns <- names(res)
    scoreCols <- getAllMergedConsCols(scoreCols, columns, mConsNames)
    
    if (!is.null(exclude))
        scoreCols <- scoreCols[!scoreCols %in% getAllMergedConsCols(exclude, columns, mConsNames)]
    
    if (length(scoreCols) > 0)
    {
        scoreRanges <- scoreRanges[scoreCols]
        res[, (scoreCols) := Map(.SD, scoreRanges, f = function(sc, scr) normalize(sc, minMaxNormalization, scr)),
            .SDcols = scoreCols]
    }
    
    return(res)
}

addElementInfoToAnnTable <- function(annTable, elements, fragElements, OM, classify)
{
    annTable <- copy(annTable)
    
    # ensure CHNOPS counts are present
    if (OM)
        elements <- union(elements, c("C", "H", "N", "O", "P", "S"))
    
    if (!is.null(elements) && length(elements) > 0)
    {
        # Retrieve element lists from formulas
        el <- getElements(annTable$neutral_formula, elements)
        annTable[, names(el) := el]
    }
    if (!is.null(fragElements) && !is.null(annTable[["frag_formula"]]) && length(fragElements) > 0)
    {
        el <- getElements(annTable$frag_formula, fragElements)
        annTable[, (paste0("frag_", names(el))) := el]
    }
    
    if (OM)
    {
        # add element ratios commonly used for plotting
        elrat <- function(el1, el2) ifelse(el2 == 0, 0, el1 / el2)
        annTable[, c("OC", "HC", "NC", "PC", "SC") :=
                     .(elrat(O, C), elrat(H, C), elrat(N, C), elrat(P, C), elrat(S, C))]
        
        # aromaticity index and related DBE (see Koch 2016, 10.1002/rcm.7433)
        annTable[, DBE_AI := 1 + C - O - S - 0.5 * (N + P + H)]
        getAI <- function(dbe, cai) ifelse(cai == 0, 0, dbe / cai)
        annTable[, AI := getAI(DBE_AI, (C - O - N - S - P))]
        
        if (classify)
            annTable[, classification := Vectorize(classifyFormula)(OC, HC, NC, AI)]
    }
    
    return(annTable)
}

doAnnotatePeakList <- function(spec, annTable, index, onlyAnnotated)
{
    if (is.null(spec))
        return(NULL)
    
    spec <- copy(spec)
    spec[, PLIndex := seq_len(nrow(spec))] # for merging
    
    if (!is.null(annTable) && nrow(annTable) >= index)
    {
        compr <- annTable[index]
        fragInfo <- compr$fragInfo[[1]]
        
        if (!is.null(fragInfo))
            spec <- merge(spec, fragInfo[, -"mz"], all.x = TRUE, by = "PLIndex")
        
        spec <- spec[, PLIndex := NULL]
    }
    
    if (onlyAnnotated)
    {
        if (is.null(spec[["ion_formula"]]))
            spec <- spec[0]
        else
            spec <- spec[!is.na(formula)]
    }
    
    return(spec[])
}

doFeatAnnConsensus <- function(obj, ..., absMinAbundance, relMinAbundance, uniqueFrom, uniqueOuter,
                               minMaxNormalization, rankWeights, annNames, uniqueCols)
{
    # NOTE: keep args in sync with sets methods
    
    allFeatAnnotations <- c(list(obj), list(...))
    
    rankWeights <- rep(rankWeights, length.out = length(allFeatAnnotations))
    relMinAbundance <- max(NULLToZero(absMinAbundance) / length(allFeatAnnotations), NULLToZero(relMinAbundance))
    
    # initialize all objects for merge: copy them, rename columns to
    # avoid duplicates and set merged by field of fragInfo.
    allAnnTables <- lapply(seq_along(allFeatAnnotations), function(anni)
    {
        mergedBy <- annNames[[anni]]
        
        return(lapply(annotations(allFeatAnnotations[[anni]]), function(ct)
        {
            ret <- copy(ct)
            
            for (r in seq_len(nrow(ret)))
            {
                fi <- ret[["fragInfo"]][[r]]
                if (!is.null(fi) && nrow(fi) > 0 && is.null(fi[["mergedBy"]]))
                {
                    fi <- copy(fi)
                    set(fi, j = "mergedBy", value = mergedBy)
                    set(ret, r, "fragInfo", list(list(fi)))
                }
            }
            
            ret[, mergedBy := annNames[anni]]
            ret[, rank := seq_len(.N)]
            
            setnames(ret, paste0(names(ret), "-", annNames[[anni]]))
            
            return(ret)
        }))
    })
    
    leftName <- annNames[[1]]
    mAnnList <- allAnnTables[[1]]
    for (compIndex in seq(2, length(allAnnTables)))
    {
        rightName <- annNames[[compIndex]]
        
        printf("Merging %s with %s... ", paste0(annNames[seq_len(compIndex-1)], collapse = ","), rightName)
        
        rightTable <- allAnnTables[[compIndex]]
        
        for (grp in union(names(mAnnList), names(rightTable)))
        {
            if (is.null(rightTable[[grp]]))
                next # nothing to merge
            else if (is.null(mAnnList[[grp]])) # not yet present
            {
                mAnnotations <- rightTable[[grp]]
                
                # rename columns that should be unique from right to left
                unCols <- c(uniqueCols, "fragInfo", "UID", "mergedBy")
                unCols <- unCols[sapply(unCols, function(uc) !is.null(mAnnotations[[paste0(uc, "-", rightName)]]))]
                setnames(mAnnotations, paste0(unCols, "-", rightName), paste0(unCols, "-", leftName))
            }
            else
            {
                mAnnotations <- merge(mAnnList[[grp]], rightTable[[grp]], by.x = paste0("UID-", leftName),
                                      by.y = paste0("UID-", rightName), all = TRUE)
                
                # merge fragment info
                fiColLeft <- paste0("fragInfo-", leftName)
                fiColRight <- paste0("fragInfo-", rightName)
                
                if (!is.null(mAnnotations[[fiColLeft]]))
                {
                    for (r in seq_len(nrow(mAnnotations)))
                    {
                        # use copy as a workaround for buggy nested data.tables
                        fiRLeft <- copy(mAnnotations[[fiColLeft]][[r]])
                        fiRRight <- copy(mAnnotations[[fiColRight]][[r]])
                        hasLeft <- length(fiRLeft) > 0 && nrow(fiRLeft) > 0
                        hasRight <- length(fiRRight) > 0 && nrow(fiRRight) > 0
                        
                        if (hasLeft && hasRight)
                        {
                            # both have fraginfo
                            fiMerged <- mergeFragInfo(fiRLeft, fiRRight, leftName, rightName)
                            set(mAnnotations, r, fiColLeft, list(list(fiMerged)))
                        }
                        else if (hasRight) # only right
                            set(mAnnotations, r, fiColLeft, list(list(fiRRight)))
                    }
                    
                    mAnnotations[, (fiColRight) := NULL]
                }
                
                # remove duplicate columns that shouldn't
                for (col in uniqueCols)
                {
                    colLeft <- paste0(col, "-", leftName)
                    colRight <- paste0(col, "-", rightName)
                    if (!is.null(mAnnotations[[colRight]]))
                    {
                        if (is.null(mAnnotations[[colLeft]]))
                            setnames(mAnnotations, colRight, colLeft)
                        else
                        {
                            mAnnotations[, (colLeft) := ifelse(!is.na(get(colLeft)), get(colLeft), get(colRight))]
                            mAnnotations[, (colRight) := NULL]
                        }
                    }
                }
                
                # collapse mergedBy
                ml <- paste0("mergedBy-", leftName); mr <- paste0("mergedBy-", rightName)
                mAnnotations[!is.na(get(ml)), (ml) := ifelse(!is.na(get(mr)), paste(get(ml), get(mr), sep = ","), get(ml))]
                mAnnotations[is.na(get(ml)), (ml) := get(mr)]
                mAnnotations[, (mr) := NULL]
            }
            
            mAnnList[[grp]] <- mAnnotations
        }
        
        cat("Done!\n")
    }
    
    printf("Determining coverage and final scores... ")
    
    # Determine coverage of between objects and the merged score. The score column can be
    # used for the former as there is guaranteed to be one for each merged object.
    for (grpi in seq_along(mAnnList))
    {
        # fix up de-duplicated column names
        deDupCols <- c(uniqueCols, c("fragInfo", "UID", "mergedBy"))
        leftCols <- paste0(deDupCols, "-", leftName)
        deDupCols <- deDupCols[leftCols %in% names(mAnnList[[grpi]])]
        leftCols <- leftCols[leftCols %in% names(mAnnList[[grpi]])]
        if (length(leftCols) > 0)
            setnames(mAnnList[[grpi]], leftCols, deDupCols)
        
        mAnnList[[grpi]][, coverage := sapply(mergedBy, function(mb) (countCharInStr(mb, ",") + 1) / length(allFeatAnnotations))]
        
        if (relMinAbundance > 0)
            mAnnList[[grpi]] <- mAnnList[[grpi]][coverage >= relMinAbundance]
        else if (!is.null(uniqueFrom))
        {
            if (!is.character(uniqueFrom))
                uniqueFrom <- annNames[uniqueFrom]
            
            keep <- function(mergedBy)
            {
                mbs <- unlist(strsplit(mergedBy, ","))
                return(all(mbs %in% uniqueFrom) && (!uniqueOuter || length(mbs) == 1))
            }
            
            mAnnList[[grpi]] <- mAnnList[[grpi]][mAnnList[[grpi]][, sapply(mergedBy, keep)]]
        }
        
        rnames <- getAllMergedConsCols("rank", names(mAnnList[[grpi]]), annNames)
        # get relevant weights with correct order
        rwInds <- unlist(lapply(annNames, grep, rnames)) # unlist: in case of no matches, sapply would yield list
        rWeights <- rankWeights[rwInds]
        ncand <- nrow(mAnnList[[grpi]])
        mAnnList[[grpi]][, rankscore := {
            invRanks <- (ncand - (unlist(.SD) - 1)) / ncand
            invRanks[is.na(invRanks)] <- 0
            weighted.mean(invRanks, rWeights)
        }, .SDcols = rnames, by = seq_len(ncand)]
        
        setorderv(mAnnList[[grpi]], "rankscore", order = -1)
        mAnnList[[grpi]][, c(rnames, "rankscore") := NULL]
    }
    
    cat("Done!\n")
    
    # prune empty/NULL results
    if (length(mAnnList) > 0)
        mAnnList <- mAnnList[sapply(mAnnList, function(r) !is.null(r) && nrow(r) > 0, USE.NAMES = FALSE)]
    
    return(mAnnList)
}