#' @include main.R
#' @include compounds.R
#' @include workflow-step-set.R
NULL

# NOTE: some methods are set in utils-feat_annotations-set.R


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

setMethod("addFormulaScoring", "compoundsSet", function(compounds, formulas, updateScore,
                                                        formulaScoreWeight)
{
    checkmate::assertClass(formulas, "formulasSet")
    
    unsetFormulas <- checkAndUnSetOther(sets(compounds), formulas, "formulas")
    compounds@setObjects <- Map(setObjects(compounds), unsetFormulas, f = addFormulaScoring,
                                MoreArgs = list(updateScore = updateScore, formulaScoreWeight = formulaScoreWeight))
    compounds <- syncAnnSetObjects(compounds, TRUE)
    
    return(compounds)
})

#' @export
setMethod("consensus", "compoundsSet", function(obj, ..., absMinAbundance = NULL, relMinAbundance = NULL,
                                                uniqueFrom = NULL, uniqueOuter = FALSE, minMaxNormalization = FALSE,
                                                rankWeights = 1, labels = NULL, setThreshold = 0,
                                                setThresholdAnn = 0.75)
{
    allAnnObjs <- c(list(obj), list(...))
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertList(allAnnObjs, types = "compoundsSet", min.len = 2, any.missing = FALSE,
                          unique = TRUE, .var.name = "...", add = ac)
    checkmate::assertCharacter(labels, min.chars = 1, len = length(allAnnObjs), null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    cons <- doFeatAnnConsensusSets(allAnnObjs, obj@origFGNames, labels, setThreshold, setThresholdAnn,
                                   absMinAbundance = absMinAbundance, relMinAbundance = relMinAbundance,
                                   uniqueFrom = uniqueFrom, uniqueOuter = uniqueOuter,
                                   minMaxNormalization = minMaxNormalization, rankWeights = rankWeights)
        

    return(compoundsConsensusSet(setObjects = cons$setObjects, setThreshold = setThreshold,
                                 setThresholdAnn = setThresholdAnn, origFGNames = obj@origFGNames,
                                 groupAnnotations = cons$groupAnnotations, scoreTypes = cons$scTypes,
                                 scoreRanges = cons$scRanges, algorithm = cons$algorithm,
                                 mergedConsensusNames = cons$mergedConsensusNames))
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
    
    cons <- makeFeatAnnSetConsensus(setObjects, names(fGroupsSet), setThreshold, setThresholdAnn, NULL)
    sc <- makeAnnSetScorings(setObjects, names(fGroupsSet))
    
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
