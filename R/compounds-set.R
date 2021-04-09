#' @include main.R
#' @include compounds.R
#' @include workflow-step-set.R
NULL

makeCompoundsSetConsensus <- function(setObjects, origFGNames, setThreshold, setThresholdAnn, mConsNames)
{
    # generate consensus by...
    # - checking setThreshold/setThresholdAnn
    # - merging by IK1
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
                   on = "InChIKey1"]

            # add missing candidates from right
            left <- rbind(left, right[!InChIKey1 %chin% left$InChIKey1], fill = TRUE)
            
            left[, setCoverage := setCoverage + 1]
            left[InChIKey1 %chin% right$InChIKey1, setCoverageAnn := setCoverageAnn + 1]
            
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

makeCompoundsSetScorings <- function(setObjects, origFGNames)
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

syncCompoundsSetObjects <- function(compoundsSet, makeCons)
{
    if (length(setObjects(compoundsSet)) >= 1)
    {
        if (makeCons)
            compoundsSet@groupAnnotations <- makeCompoundsSetConsensus(compoundsSet@setObjects, compoundsSet@origFGNames,
                                                                       compoundsSet@setThreshold, compoundsSet@setThresholdAnn,
                                                                       mergedConsensusNames(compoundsSet))
        else
        {
            # sync available feature groups
            allFGroups <- unique(unlist(lapply(setObjects(compoundsSet), groupNames)))
            compoundsSet@groupAnnotations <- compoundsSet@groupAnnotations[intersect(groupNames(compoundsSet), allFGroups)]
            
            # only keep results from sets still present
            spat <- paste0(sets(compoundsSet), collapse = "|")
            compoundsSet@groupAnnotations <- lapply(compoundsSet@groupAnnotations, function(ct) ct[grepl(spat, set)])
        }
    }
    else
        compoundsSet@groupAnnotations <- list()
    
    compoundsSet@scoreRanges <- compoundsSet@scoreRanges[groupNames(compoundsSet)]
    
    # update scoreTypes/scoreRanges
    sc <- makeCompoundsSetScorings(setObjects(compoundsSet), compoundsSet@origFGNames)
    compoundsSet@scoreTypes <- sc$scTypes
    compoundsSet@scoreRanges <- sc$scRanges
    
    return(compoundsSet)
}

compoundsSet <- setClass("compoundsSet", slots = c(setThreshold = "numeric", setThresholdAnn = "numeric",
                                                   origFGNames = "character"),
                        contains = c("compounds", "workflowStepSet"))

compoundsConsensusSet <- setClass("compoundsConsensusSet", slots = c(mergedConsensusNames = "character"),
                                  contains = "compoundsSet")
setMethod("mergedConsensusNames", "compoundsConsensusSet", function(obj) obj@mergedConsensusNames)


#' @describeIn compoundsSet Shows summary information for this object.
#' @export
setMethod("show", "compoundsSet", function(object)
{
    callAllNextMethods(object, show, firstClass = "compounds", startFrom = "compoundsSet")
})

setMethod("[", c("compoundsSet", "ANY", "missing", "missing"), function(x, i, j, ..., sets = NULL, drop = TRUE)
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
        x <- syncCompoundsSetObjects(x, FALSE)
    
    return(x)
})

#' @export
setMethod("filter", "compoundsSet", function(obj, ..., negate = FALSE, sets = NULL)
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
        obj <- syncCompoundsSetObjects(obj, TRUE)
        cat("Done!\n")
    }
    
    return(obj)
})

#' @export
setMethod("plotSpectrum", "compoundsSet", function(obj, index, groupName, MSPeakLists, formulas = NULL,
                                                   plotStruct = TRUE, title = NULL,
                                                   specSimParams = getDefSpecSimParams(), useGGPlot2 = FALSE,
                                                   mincex = 0.9, xlim = NULL, ylim = NULL, maxMolSize = c(0.2, 0.4),
                                                   molRes = c(100, 100), perSet = TRUE, mirror = TRUE, ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertIntegerish(index, lower = 1, min.len = 1, max.len = 2, any.missing = FALSE, add = ac)
    checkmate::assertCharacter(groupName, min.len = 1, max.len = 2, min.chars = 1, add = ac)
    if (length(index) != length(groupName))
        stop("Lengths of index and groupName should be equal.")
    assertSpecSimParams(specSimParams, add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakListsSet", add = ac)
    checkmate::assertClass(formulas, "formulasSet", null.ok = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ plotStruct + useGGPlot2 + perSet + mirror, fixed = list(add = ac))
    checkmate::assertNumber(mincex, lower = 0, finite = TRUE, add = ac)
    assertXYLim(xlim, ylim, add = ac)
    aapply(checkmate::assertNumeric, . ~ maxMolSize + molRes, finite = TRUE, len = 2, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    if (!perSet || length(sets(obj)) == 1)
        return(callNextMethod(obj, index, groupName, MSPeakLists, formulas, plotStruct, title, specSimParams,
                              useGGPlot2, mincex, xlim, ylim, maxMolSize, molRes, ...))

    if (length(groupName) == 1)
    {
        if (index > nrow(obj[[groupName]]))
            stop(sprintf("Specified candidate index out of range %d/%d", index, nrow(obj[[groupName]])), call. = FALSE)
        
        compr <- obj[[groupName]][index, ]
        mol <- NULL
        if (plotStruct)
        {
            mol <- getMoleculesFromSMILES(compr$SMILES)
            if (!isValidMol(mol))
                mol <- NULL
        }
        
        if (is.null(title))
            title <- getCompoundsSpecPlotTitle(compr$compoundName, compr$neutral_formula)
        
        spec <- annotatedPeakList(obj, index, groupName, MSPeakLists, formulas)
        if (is.null(spec))
            return(NULL)
        
        specs <- split(spec, by = "set")
        # UNDONE: this will overwrite consensus algo if present, OK?
        specs <- lapply(specs, function(x) x[!is.na(formula), mergedBy := set])
        
        plotData <- getMSPlotDataOverlay(specs, mirror, TRUE, 1, NULL)
        return(makeMSPlotOverlay(plotData, title, mincex, xlim, ylim, useGGPlot2, ...,  mol = mol,
                                 maxMolSize = maxMolSize, molRes = molRes))
    }
    else
    {
        if (plotStruct)
            stop("Cannot plot structure when comparing spectra") # UNDONE?
        
        for (i in seq_len(2))
        {
            if (index[i] > nrow(obj[[groupName[i]]]))
                stop(sprintf("Specified candidate index out of range %d/%d", index[i], nrow(obj[[groupName[i]]])),
                     call. = FALSE)
        }
        
        if (is.null(title))
        {
            compr1 <- obj[[groupName[1]]][index[1], ]; compr2 <- obj[[groupName[2]]][index[2], ]
            title <- getCompoundsSpecPlotTitle(compr1$compoundName, compr1$neutral_formula, compr2$compoundName,
                                               compr2$neutral_formula)
        }
        
        theSets <- sets(obj)
        
        usObj <- sapply(theSets, unset, obj = obj, simplify = FALSE)
        
        # check which sets actually contain requested data
        theSets <- theSets[sapply(theSets, function(s) all(groupName %in% groupNames(usObj[[s]])))]
        if (length(theSets) == 0)
            return(NULL)
        usObj <- usObj[theSets]
        
        usMSPL <- checkAndUnSetOther(theSets, MSPeakLists, "MSPeakLists")
        binnedPLs <- Map(usMSPL, theSets, f = getBinnedPLPair,
                         MoreArgs = list(groupNames = groupName, analyses = NULL, MSLevel = 2,
                                         specSimParams = specSimParams, mustExist = FALSE))
        
        if (!is.null(formulas))
            usForm <- checkAndUnSetOther(theSets, formulas, "formulas")
        else
            usForm <- rep(list(NULL), length(theSets))
        
        mergeBinnedAnn <- function(nr)
        {
            # convert candidate index to non set version by using the ranks
            usInds <- lapply(theSets, function(s) obj[[groupName[nr]]][[paste0("rank-", s)]][index[nr]])
            binPLs <- sapply(binnedPLs, "[[", nr, simplify = FALSE)
            annPLs <- Map(usObj, usInds, usMSPL, usForm, f = annotatedPeakList,
                          MoreArgs = list(groupName = groupName[nr]))
            annPLs <- Map(mergeBinnedAndAnnPL, binPLs, annPLs, MoreArgs = list(which = nr))
            annPLs <- rbindlist(annPLs, idcol = "set", fill = TRUE)
            return(annPLs)
        }
        
        topSpec <- mergeBinnedAnn(1); bottomSpec <- mergeBinnedAnn(2)
        allSpectra <- rbind(topSpec, bottomSpec, fill = TRUE)
        
        specs <- split(allSpectra, by = "which")
        plotData <- getMSPlotDataOverlay(specs, mirror, FALSE, 2, "overlap")
        
        makeMSPlotOverlay(plotData, title, mincex, xlim, ylim, useGGPlot2, ...)
    }
})

setMethod("plotSpectrumHash", "compoundsSet", function(obj, index, groupName, MSPeakLists, formulas = NULL,
                                                       plotStruct = TRUE, title = NULL,
                                                       specSimParams = getDefSpecSimParams(), useGGPlot2 = FALSE,
                                                       mincex = 0.9, xlim = NULL, ylim = NULL,
                                                       maxMolSize = c(0.2, 0.4), molRes = c(100, 100),
                                                       perSet = TRUE, mirror = TRUE, ...)
{
    return(makeHash(callNextMethod(obj, index, groupName, MSPeakLists, formulas, plotStruct, title, specSimParams,
                                   useGGPlot2, mincex, xlim, ylim, maxMolSize, molRes, ...),
                    perSet, mirror))
})

setMethod("annotatedPeakList", "compoundsSet", function(obj, ...)
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

setMethod("addFormulaScoring", "compoundsSet", function(compounds, formulas, updateScore,
                                                        formulaScoreWeight)
{
    checkmate::assertClass(formulas, "formulasSet")
    
    unsetFormulas <- checkAndUnSetOther(sets(compounds), formulas, "formulas")
    compounds@setObjects <- Map(setObjects(compounds), unsetFormulas, f = addFormulaScoring,
                                MoreArgs = list(updateScore = updateScore, formulaScoreWeight = formulaScoreWeight))
    compounds <- syncCompoundsSetObjects(compounds, TRUE)
    
    return(compounds)
})

#' @export
setMethod("consensus", "compoundsSet", function(obj, ..., absMinAbundance = NULL, relMinAbundance = NULL,
                                                uniqueFrom = NULL, uniqueOuter = FALSE, minMaxNormalization = FALSE,
                                                rankWeights = 1, labels = NULL, setThreshold = 0,
                                                setThresholdAnn = 0.75)
{
    # make consensus of shared setObjects
    # add unique setObjects
    # make 'regular' set consensus from new results
    
    allCompounds <- c(list(obj), list(...))
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertList(allCompounds, types = "compoundsSet", min.len = 2, any.missing = FALSE,
                          unique = TRUE, .var.name = "...", add = ac)
    checkmate::assertNumeric(rankWeights, lower = 0, finite = TRUE, add = ac)
    checkmate::assertCharacter(labels, min.chars = 1, len = length(allCompounds), null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (!allSame(lapply(allCompounds, sets)))
        stop("All objects must have the same sets.")
    if (!allSame(lapply(allCompounds, slot, "origFGNames")))
        stop("All objects must have been generated from the same feature groups.")
    
    # NOTE: don't want to keep -set suffix
    compNames <- if (!is.null(labels)) labels else sub("\\-set$", "", sapply(allCompounds, algorithm))
    
    setObjects <- sapply(sets(obj), function(set)
    {
        return(do.call(consensus, c(lapply(lapply(allCompounds, setObjects), "[[", set),
                                    list(absMinAbundance = absMinAbundance, relMinAbundance = relMinAbundance,
                                         uniqueFrom = uniqueFrom, uniqueOuter = uniqueOuter,
                                         minMaxNormalization = minMaxNormalization, rankWeights = rankWeights,
                                         labels = compNames))))
    }, simplify = FALSE)
    
    cons <- makeCompoundsSetConsensus(setObjects, obj@origFGNames, setThreshold, setThresholdAnn, compNames)
    sc <- makeCompoundsSetScorings(setObjects, obj@origFGNames)
    
    return(compoundsConsensusSet(setObjects = setObjects, setThreshold = setThreshold, setThresholdAnn = setThresholdAnn,
                                 origFGNames = obj@origFGNames, groupAnnotations = cons, scoreTypes = sc$scTypes,
                                 scoreRanges = sc$scRanges,
                                 algorithm = paste0(unique(sapply(allCompounds, algorithm)), collapse = ","),
                                 mergedConsensusNames = compNames))
})


generateCompoundsSet <- function(fGroupsSet, MSPeakListsSet, generator, ..., setThreshold, setThresholdAnn)
{
    aapply(checkmate::assertNumber, . ~ setThreshold + setThresholdAnn, lower = 0, upper = 1.0, finite = TRUE)
    msplArgs <- assertAndGetMSPLSetsArgs(fGroupsSet, MSPeakListsSet)
    
    unsetFGroupsList <- sapply(sets(fGroupsSet), unset, obj = fGroupsSet, simplify = FALSE)
    setObjects <- Map(unsetFGroupsList, msplArgs,
                      f = function(fg, mspl) generator(fGroups = fg, MSPeakLists = mspl[[1]], ...))
    
    # update fragInfos
    for (s in names(setObjects))
    {
        for (fg in groupNames(setObjects[[s]]))
        {
            pl <- copy(MSPeakListsSet[[fg]][["MSMS"]]); pl[, PLIndex := seq_len(.N)]; pl <- pl[set == s]
            ct <- setObjects[[s]]@groupAnnotations[[fg]]
            ct[, fragInfo := lapply(fragInfo, function(fi)
            {
                fi <- copy(fi) # BUG: avoid warning that somehow was incorrectly copied (invalid .internal.selfref)
                if (nrow(fi) == 0)
                {
                    # otherwise it will be assigned as empty list, which messes up merging elsewhere
                    fi[, PLIndexSet := numeric()]
                }
                else
                    fi[, PLIndexSet := sapply(mz, function(fimz) pl[which.min(abs(fimz - mz))][["PLIndex"]])]
                fi[, set := s]
                return(fi)
            })]
            setObjects[[s]]@groupAnnotations[[fg]] <- ct            
        }
    }
    
    cons <- makeCompoundsSetConsensus(setObjects, names(fGroupsSet), setThreshold, setThresholdAnn, NULL)
    sc <- makeCompoundsSetScorings(setObjects, names(fGroupsSet))
    
    return(compoundsSet(setObjects = setObjects, setThreshold = setThreshold, setThresholdAnn = setThresholdAnn,
                        origFGNames = names(fGroupsSet), groupAnnotations = cons, scoreTypes = sc$scTypes,
                        scoreRanges = sc$scRanges, algorithm = makeSetAlgorithm(setObjects)))
}


compoundsUnset <- setClass("compoundsUnset", contains = "compounds")
setMethod("unset", "compoundsSet", function(obj, set)
{
    assertSets(obj, set, FALSE)
    obj <- obj[, sets = set]
    
    cList <- lapply(annotations(obj), copy)
    
    # get rid of sets specific columns
    cList <- lapply(cList, data.table::set, j = c("set", "setCoverage", "setCoverageAnn"), value = NULL)
    
    # ... and in fragInfos
    cList <- lapply(cList, function(ct)
    {
        ct[, fragInfo := lapply(fragInfo, function(fi)
        {
            fi <- copy(fi)
            fi <- fi[set == ..set]
            set(fi, j = c("set", "PLIndex"), value = NULL)
            setnames(fi, "PLIndexOrig", "PLIndex") # restore original index
        })]
    })
    
    return(compoundsUnset(groupAnnotations = cList, scoreTypes = obj@scoreTypes, scoreRanges = obj@scoreRanges,
                          algorithm = paste0(algorithm(obj), "_unset")))
})
