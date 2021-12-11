#' @include main.R
#' @include compounds.R
#' @include workflow-step-set.R
#' @include utils-feat_annotations-set.R
NULL

# NOTE: some methods are set in utils-feat_annotations-set.R

# NOTE: can't break the long line below
#' @templateVar class compoundsSet
#' @templateVar parent compounds
#' @templateVar generator generateCompounds
#' @templateVar classUnset compoundsUnset
#' @templateVar exObj compounds
#' @templateVar consAlgo2 MetFrag
#' @templateVar extraMethods \item \code{addFormulaScoring} Adds the formula scorings to the original data and re-creates the annotation set consensus (see below for implications).
#' @template featAnnSets-class
#' 
#' @rdname compounds-class
#' @export
compoundsSet <- setClass("compoundsSet", slots = c(setThreshold = "numeric", setThresholdAnn = "numeric",
                                                   origFGNames = "character"),
                        contains = c("compounds", "workflowStepSet"))

#' @rdname compounds-class
#' @export
compoundsConsensusSet <- setClass("compoundsConsensusSet", slots = c(mergedConsensusNames = "character"),
                                  contains = "compoundsSet")


setMethod("updateSetConsensus", "compoundsSet", function(obj)
{
    obj <- doUpdateSetConsensus(obj)
    
    # update scoreTypes/scoreRanges
    sc <- makeAnnSetScorings(setObjects(obj), obj@origFGNames)
    obj@scoreTypes <- sc$scTypes
    obj@scoreRanges <- sc$scRanges
    
    return(obj)
})

setMethod("mergedConsensusNames", "compoundsSet", doFeatAnnMCNSets)
setMethod("mergedConsensusNames", "compoundsConsensusSet", doFeatAnnMCNSetsCons)

#' @rdname compounds-class
#' @export
setMethod("show", "compoundsSet", function(object)
{
    callAllNextMethods(object, show, firstClass = "compounds", startFrom = "compoundsSet")
})

#' @rdname compounds-class
#' @export
setMethod("delete", "compoundsSet", doFeatAnnDeleteSets)

#' @param i,j,drop Passed to the \code{\link[=[,featureAnnotations,ANY,missing,missing-method]{featureAnnotations}}
#'   method.
#' @rdname compounds-class
#' @export
setMethod("[", c("compoundsSet", "ANY", "missing", "missing"), doFeatAnnSubsetSets)

#' @rdname compounds-class
#' @param negate Passed to the \code{\link[=filter,featureAnnotations-method]{featureAnnotations}} method.
#' @export
setMethod("filter", "compoundsSet", doFeatAnnFilterSets)

#' @rdname compounds-class
#' @export
setMethod("plotSpectrum", "compoundsSet", function(obj, index, groupName, MSPeakLists, formulas = NULL,
                                                   plotStruct = FALSE, title = NULL,
                                                   specSimParams = getDefSpecSimParams(),
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
    aapply(checkmate::assertFlag, . ~ plotStruct +  perSet + mirror, fixed = list(add = ac))
    checkmate::assertNumber(mincex, lower = 0, finite = TRUE, add = ac)
    assertXYLim(xlim, ylim, add = ac)
    aapply(checkmate::assertNumeric, . ~ maxMolSize + molRes, finite = TRUE, len = 2, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    if (!perSet || length(sets(obj)) == 1)
        return(callNextMethod(obj, index, groupName, MSPeakLists, formulas, plotStruct, title, specSimParams,
                              mincex, xlim, ylim, maxMolSize, molRes, ...))

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
        specs <- lapply(specs, function(x)
        {
            x[annotated == TRUE, mergedBy := set]
            return(x)
        })
        
        plotData <- getMSPlotDataOverlay(specs, mirror, TRUE, 1, NULL)
        return(makeMSPlotOverlay(plotData, title, mincex, xlim, ylim, ...,  mol = mol, maxMolSize = maxMolSize,
                                 molRes = molRes))
    }
    else
    {
        if (plotStruct)
            stop("Cannot plot structure when comparing spectra") # UNDONE?
        
        for (i in seq_len(2))
        {
            if (is.null(obj[[groupName[i]]]))
                stop(paste("No data for specified feature group:", groupName[i]), call. = FALSE)
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
        
        makeMSPlotOverlay(plotData, title, mincex, xlim, ylim, ...)
    }
})

setMethod("plotSpectrumHash", "compoundsSet", function(obj, index, groupName, MSPeakLists, formulas = NULL,
                                                       plotStruct = FALSE, title = NULL,
                                                       specSimParams = getDefSpecSimParams(),
                                                       mincex = 0.9, xlim = NULL, ylim = NULL,
                                                       maxMolSize = c(0.2, 0.4), molRes = c(100, 100),
                                                       perSet = TRUE, mirror = TRUE, ...)
{
    return(makeHash(callNextMethod(obj, index, groupName, MSPeakLists, formulas, plotStruct, title, specSimParams,
                                   mincex, xlim, ylim, maxMolSize, molRes, ...),
                    perSet, mirror))
})

#' @rdname compounds-class
#' @export
setMethod("addFormulaScoring", "compoundsSet", function(compounds, formulas, updateScore,
                                                        formulaScoreWeight)
{
    checkmate::assertClass(formulas, "formulasSet")
    
    unsetFormulas <- checkAndUnSetOther(sets(compounds), formulas, "formulas")
    compounds@setObjects <- Map(setObjects(compounds), unsetFormulas, f = addFormulaScoring,
                                MoreArgs = list(updateScore = updateScore, formulaScoreWeight = formulaScoreWeight))
    compounds <- updateSetConsensus(compounds)
    
    return(compounds)
})

#' @rdname compounds-class
#' @export
setMethod("annotatedPeakList", "compoundsSet", function(obj, index, groupName, MSPeakLists, formulas = NULL, ...)
{
    checkmate::assertClass(formulas, "formulasSet", null.ok = TRUE)
    allFGroups <- groupNames(obj)
    if (!is.null(formulas))
        allFGroups <- union(allFGroups, groupNames(formulas))
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertCount(index, positive = TRUE, add = ac)
    assertChoiceSilent(groupName, allFGroups, add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakListsSet", add = ac)
    checkmate::reportAssertions(ac)
    
    return(doAnnotatePeakListSet(obj, index, groupName, MSPeakLists, formulas, ...))
})

#' @rdname compounds-class
#' @export
setMethod("consensus", "compoundsSet", function(obj, ..., absMinAbundance = NULL, relMinAbundance = NULL,
                                                uniqueFrom = NULL, uniqueOuter = FALSE, rankWeights = 1, labels = NULL,
                                                filterSets = FALSE, setThreshold = 0, setThresholdAnn = 0)
{
    allAnnObjs <- c(list(obj), list(...))
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertList(allAnnObjs, types = "compoundsSet", min.len = 2, any.missing = FALSE,
                          unique = TRUE, .var.name = "...", add = ac)
    checkmate::assertCharacter(labels, min.chars = 1, len = length(allAnnObjs), null.ok = TRUE, add = ac)
    checkmate::assertFlag(filterSets, add = ac)
    aapply(checkmate::assertNumber, . ~ setThreshold + setThresholdAnn, lower = 0, upper = 1, finite = TRUE)
    checkmate::reportAssertions(ac)

    labels <- prepareConsensusLabels(obj, ..., labels = labels)
    
    assertConsCommonArgs(absMinAbundance, relMinAbundance, uniqueFrom, uniqueOuter, labels)
    
    cons <- doFeatAnnConsensusSets(allAnnObjs, labels, setThreshold, setThresholdAnn, rankWeights)
    sc <- makeAnnSetScorings(cons$setObjects, cons$origFGNames)
    

    ret <- compoundsConsensusSet(setObjects = cons$setObjects, setThreshold = setThreshold,
                                 setThresholdAnn = setThresholdAnn, origFGNames = cons$origFGNames,
                                 groupAnnotations = cons$groupAnnotations, scoreTypes = sc$scTypes,
                                 scoreRanges = sc$scRanges, algorithm = cons$algorithm,
                                 mergedConsensusNames = cons$mergedConsensusNames)
    
    ret <- filterFeatAnnConsensus(ret, absMinAbundance, relMinAbundance, uniqueFrom, uniqueOuter, filterSets)
    
    return(ret)
})


generateCompoundsSet <- function(fGroupsSet, MSPeakListsSet, adduct, generator, ..., setThreshold, setThresholdAnn)
{
    aapply(checkmate::assertNumber, . ~ setThreshold + setThresholdAnn, lower = 0, upper = 1.0, finite = TRUE)
    msplArgs <- assertAndGetMSPLSetsArgs(fGroupsSet, MSPeakListsSet)
    verifyNoAdductIonizationArg(adduct)
    
    unsetFGroupsList <- sapply(sets(fGroupsSet), unset, obj = fGroupsSet, simplify = FALSE)
    setObjects <- Map(unsetFGroupsList, msplArgs,
                      f = function(fg, mspl) generator(fGroups = fg, MSPeakLists = mspl[[1]], adduct = NULL, ...))
    setObjects <- initSetFragInfos(setObjects, MSPeakListsSet)

    cons <- makeFeatAnnSetConsensus(setObjects, names(fGroupsSet), setThreshold, setThresholdAnn, NULL)
    sc <- makeAnnSetScorings(setObjects, names(fGroupsSet))
    
    return(compoundsSet(setObjects = setObjects, setThreshold = setThreshold, setThresholdAnn = setThresholdAnn,
                        origFGNames = names(fGroupsSet), groupAnnotations = cons, scoreTypes = sc$scTypes,
                        scoreRanges = sc$scRanges, algorithm = makeSetAlgorithm(setObjects)))
}


#' @rdname compounds-class
#' @export
compoundsUnset <- setClass("compoundsUnset", contains = "compounds")

#' @rdname compounds-class
#' @export
setMethod("unset", "compoundsSet", function(obj, set)
{
    assertSets(obj, set, FALSE)
    uann <- doFeatAnnUnset(obj, set)
    return(compoundsUnset(groupAnnotations = uann, scoreTypes = obj@scoreTypes, scoreRanges = obj@scoreRanges,
                          algorithm = paste0(algorithm(obj), "_unset")))
})

#' @rdname compounds-class
#' @export
setMethod("unset", "compoundsConsensusSet", function(obj, set)
{
    # get rid of overall consensus cols, as they interfere when set specific are renamed in the parent unset method
    obj@groupAnnotations <- lapply(obj@groupAnnotations, function(annTable)
    {
        annTable <- copy(annTable)
        annTable[, c("coverage", "mergedBy") := NULL]
        return(annTable)
    })
    
    return(callNextMethod(obj, set))
})
