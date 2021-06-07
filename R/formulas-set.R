#' @include main.R
#' @include formulas.R
#' @include workflow-step-set.R
#' @include utils-feat_annotations-set.R
NULL

# NOTE: some methods are set in utils-feat_annotations-set.R

#' @param set \setsWF The name of the set.
#' @param sets \setsWF For \code{[} and \code{filter}: a \code{character} with name(s) of the sets to keep (or remove if
#'   \code{negate=TRUE}). Note: if \code{updateConsensus=FALSE} then the \code{setCoverage} column of the annotation
#'   results is not updated.
#' @param updateConsensus \setsWF If \code{TRUE} then the annonation consensus among set results is updated. See the
#'   \verb{Sets workflows} section for more details.
#' @param perSet,mirror \setsWF If \code{perSet=TRUE} then the set specific mass peaks are annotated separately.
#'   Furthermore, if \code{mirror=TRUE} (and there are two sets in the object) then a mirror plot is generated.
#' @param setThreshold,setThresholdAnn \setsWF Thresholds used to create the annotation set consensus. See
#'   \code{\link{generateFormulas}}.
#'
#' @slot setThreshold,setThresholdAnn \setsWF A copy of the equally named arguments that were passed when this object
#'   was created by \code{\link{generateFormulas}}.
#' @slot origFGNames \setsWF The original (order of) names of the \code{\link{featureGroups}} object that was used to
#'   create this object.
#'
#' @section Sets workflows: \setsWFClass{formulasSet}{formulas}
#'
#'   \setsWFNewMethodsSO{formulasUnset}{Only the annotation results that are present in the specified set are kept
#'   (based on the set consensus, see below for implications).}
#'
#'   \setsWFChangedMethods{
#'
#'   \item \code{filter} and the subset operator (\code{[}) Can be used to select data that is only present for selected
#'   sets. Depending on the \code{updateConsenus}, both either operate on set consensus or original data (see below for
#'   implications).
#'
#'   \item \code{annotatedPeakList} Returns a combined annotation table with all sets.
#'
#'   \item \code{plotSpectrum} Is able to highlight set specific mass peaks (\code{perSet} and \code{mirror} arguments).
#'
#'   \item \code{consensus} Creates the algorithm consensus based on the original annotation data (see below for
#'   implications). Furthermore, like the sets workflow method for \code{\link{generateFormulas}}, this method supports
#'   the \code{setThreshold} and \code{setThresholdAnn} arguments.
#'
#'   }
#'
#'   Two types of annotation data are stored in a \code{formulasSet} object: \enumerate{
#'
#'   \item Annotations that are produced from a consensus between set results (see \code{generateFormulas}).
#'
#'   \item The 'original' annotation data per set, prior to when the set consensus was made. This includes candidates
#'   that were filtered out because of the thresholds set by \code{setThreshold} and \code{setThresholdAnn}. However,
#'   when \code{filter} or subsetting (\code{[}) operations are performed, the original data is also updated.
#'
#'   }
#'
#'   In most cases the first data is used. However, in a few cases the original annotation data is used (as indicated
#'   above), for instance, to re-create the set consensus. It is important to realize that the original annotation data
#'   may have \emph{additional} candidates, and a newly created set consensus may therefore have 'new' candidates. For
#'   instance, when the object consists of the sets \code{"positive"} and \code{"negative"} and \code{setThreshold=1}
#'   was used to create it, then \code{formulas[, sets = "positive", updateConsensus = TRUE]} may now have additional
#'   candidates, \emph{i.e.} those that were not present in the \code{"negative"} set and were previously removed due to
#'   the consensus threshold filter.
#'
#' @rdname formulas-class
#' @export
formulasSet <- setClass("formulasSet", slots = c(setThreshold = "numeric",
                                                 setThresholdAnn = "numeric",
                                                 origFGNames = "character"),
                        contains = c("formulas", "workflowStepSet"))

#' @rdname formulas-class
#' @export
formulasConsensusSet <- setClass("formulasConsensusSet", slots = c(mergedConsensusNames = "character"),
                                 contains = "formulasSet")


setMethod("updateSetConsensus", "formulasSet", function(obj)
{
    obj <- doUpdateSetConsensus(obj)
    
    # update feature formulas
    obj@featureFormulas <- pruneList(lapply(obj@featureFormulas, function(ff) pruneList(ff[groupNames(obj)])), TRUE)
    
    return(obj)
})

setMethod("mergedConsensusNames", "formulasSet", doFeatAnnMCNSets)
setMethod("mergedConsensusNames", "formulasConsensusSet", doFeatAnnMCNSetsCons)

#' @rdname formulas-class
#' @export
setMethod("show", "formulasSet", function(object)
{
    callAllNextMethods(object, show, firstClass = "formulas", startFrom = "formulasSet")
})

#' @rdname formulas-class
#' @export
setMethod("delete", "formulasSet", doFeatAnnDeleteSets)

#' @rdname formulas-class
#' @export
setMethod("[", "formulasSet", doFeatAnnSubsetSets)

#' @rdname formulas-class
#' @export
setMethod("filter", "formulasSet", doFeatAnnFilterSets)

#' @rdname formulas-class
#' @export
setMethod("plotSpectrum", "formulasSet", function(obj, index, groupName, analysis = NULL, MSPeakLists,
                                                  title = NULL, specSimParams = getDefSpecSimParams(),
                                                  useGGPlot2 = FALSE, mincex = 0.9, xlim = NULL,
                                                  ylim = NULL, perSet = TRUE, mirror = TRUE, ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertIntegerish(index, lower = 1, min.len = 1, max.len = 2, any.missing = FALSE, add = ac)
    checkmate::assertCharacter(groupName, min.len = 1, max.len = 2, min.chars = 1, add = ac)
    checkmate::assertCharacter(analysis, min.len = 1, max.len = 2, min.chars = 1, null.ok = TRUE, add = ac)
    if (length(index) != length(groupName))
        stop("Lengths of index and groupName should be equal.")
    if (!is.null(analysis) && length(analysis) != length(groupName))
        stop("Lengths of analysis and groupName should be equal.")
    assertSpecSimParams(specSimParams, add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakListsSet", add = ac)
    checkmate::assertString(title, null.ok = TRUE, add = ac)
    checkmate::assertNumber(mincex, lower = 0, finite = TRUE, add = ac)
    assertXYLim(xlim, ylim, add = ac)
    aapply(checkmate::assertFlag, . ~ useGGPlot2 + perSet + mirror, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    if (!perSet || length(sets(obj)) == 1 || !is.null(analysis))
        return(callNextMethod(obj, index, groupName, analysis, MSPeakLists, title, specSimParams = specSimParams,
                              useGGPlot2, mincex, xlim, ylim, ...))
    
    if (length(groupName) == 1)
    {
        spec <- annotatedPeakList(obj, index, groupName, analysis, MSPeakLists)
        if (is.null(spec))
            return(NULL)

        if (index > nrow(obj[[groupName]]))
            stop(sprintf("Specified candidate index out of range %d/%d", index, nrow(obj[[groupName]])), call. = FALSE)

        if (is.null(title))
            title <- subscriptFormula(obj[[groupName]]$neutral_formula[index])
        
        specs <- split(spec, by = "set")
        # UNDONE: this will overwrite consensus algo if present, OK?
        specs <- lapply(specs, function(x)
        {
            if (!is.null(x[["ion_formula"]]))
                x[!is.na(ion_formula), mergedBy := set]
            return(x)
        })
        
        plotData <- getMSPlotDataOverlay(specs, mirror, TRUE, 1, NULL)
        return(makeMSPlotOverlay(plotData, title, mincex, xlim, ylim, useGGPlot2, ...))
    }
    else
    {
        for (i in seq_len(2))
        {
            if (index[i] > nrow(obj[[groupName[i]]]))
                stop(sprintf("Specified candidate index out of range %d/%d", index[i], nrow(obj[[groupName[i]]])),
                     call. = FALSE)
        }
        
        if (is.null(title))
            title <- subscriptFormula(obj[[groupName[1]]]$neutral_formula[index[1]],
                                      formulas2 = obj[[groupName[2]]]$neutral_formula[index])        
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

        mergeBinnedAnn <- function(nr)
        {
            binPLs <- sapply(binnedPLs, "[[", nr, simplify = FALSE)
            annPLs <- Map(usObj, usMSPL, f = annotatedPeakList,
                          MoreArgs = list(index = index[nr], groupName = groupName[nr], analysis = analysis[nr]))
            annPLs <- Map(mergeBinnedAndAnnPL, binPLs, annPLs, MoreArgs = list(which = nr))
            annPLs <- rbindlist(annPLs, idcol = "set")
            return(annPLs)
        }
        
        topSpec <- mergeBinnedAnn(1); bottomSpec <- mergeBinnedAnn(2)
        allSpectra <- rbind(topSpec, bottomSpec)
        
        specs <- split(allSpectra, by = "which")
        plotData <- getMSPlotDataOverlay(specs, mirror, FALSE, 2, "overlap")
        
        makeMSPlotOverlay(plotData, title, mincex, xlim, ylim, useGGPlot2, ...)
    }
})

setMethod("plotSpectrumHash", "formulasSet", function(obj, index, groupName, analysis = NULL, MSPeakLists,
                                                      title = NULL, specSimParams = getDefSpecSimParams(),
                                                      useGGPlot2 = FALSE, mincex = 0.9, xlim = NULL, ylim = NULL,
                                                      perSet = TRUE, mirror = TRUE, ...)
{
    return(makeHash(callNextMethod(obj, index, groupName, analysis, MSPeakLists,
                                   title, specSimParams, useGGPlot2, mincex, xlim, ylim, ...),
                    perSet, mirror))
})

#' @rdname formulas-class
#' @export
setMethod("annotatedPeakList", "formulasSet", function(obj, index, groupName, analysis = NULL, MSPeakLists, ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertCount(index, positive = TRUE, add = ac)
    assertChoiceSilent(groupName, groupNames(obj), add = ac)
    checkmate::assertString(analysis, min.chars = 1, null.ok = TRUE, add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakListsSet", add = ac)
    checkmate::reportAssertions(ac)

    return(doAnnotatePeakListSet(obj, index, groupName, MSPeakLists, formulas = NULL, analysis = analysis, ...))
})

#' @rdname formulas-class
#' @export
setMethod("consensus", "formulasSet", function(obj, ..., absMinAbundance = NULL, relMinAbundance = NULL,
                                               uniqueFrom = NULL, uniqueOuter = FALSE, rankWeights = 1, labels = NULL,
                                               filterSets = FALSE, setThreshold = 0, setThresholdAnn = 0.75)
{
    allAnnObjs <- c(list(obj), list(...))
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertList(allAnnObjs, types = "formulasSet", min.len = 2, any.missing = FALSE,
                          unique = TRUE, .var.name = "...", add = ac)
    checkmate::assertCharacter(labels, min.chars = 1, len = length(allAnnObjs), null.ok = TRUE, add = ac)
    checkmate::assertFlag(filterSets, add = ac)
    aapply(checkmate::assertNumber, . ~ setThreshold + setThresholdAnn, lower = 0, upper = 1, finite = TRUE)
    checkmate::reportAssertions(ac)
    
    labels <- prepareConsensusLabels(obj, ..., labels = labels)
    
    assertConsCommonArgs(absMinAbundance, relMinAbundance, uniqueFrom, uniqueOuter, labels)
    
    cons <- doFeatAnnConsensusSets(allAnnObjs, labels, setThreshold, setThresholdAnn, rankWeights)
    combFormulas <- Reduce(modifyList, lapply(cons$setObjects, annotations, features = TRUE))
    
    ret <- formulasConsensusSet(setObjects = cons$setObjects, setThreshold = setThreshold,
                                setThresholdAnn = setThresholdAnn, origFGNames = cons$origFGNames,
                                groupAnnotations = cons$groupAnnotations, featureFormulas = combFormulas,
                                algorithm = cons$algorithm, mergedConsensusNames = cons$mergedConsensusNames)
    
    ret <- filterFeatAnnConsensus(ret, absMinAbundance, relMinAbundance, uniqueFrom, uniqueOuter, filterSets)
    
    return(ret)
})


generateFormulasSet <- function(fGroupsSet, MSPeakListsSet, generator, ..., setThreshold, setThresholdAnn)
{
    aapply(checkmate::assertNumber, . ~ setThreshold + setThresholdAnn, lower = 0, upper = 1, finite = TRUE)
    msplArgs <- assertAndGetMSPLSetsArgs(fGroupsSet, MSPeakListsSet)
    
    # UNDONE: mention that adduct argument is automatically set

    unsetFGroupsList <- sapply(sets(fGroupsSet), unset, obj = fGroupsSet, simplify = FALSE)
    setObjects <- Map(unsetFGroupsList, msplArgs,
                      f = function(fg, mspl) generator(fGroups = fg, MSPeakLists = mspl[[1]], ...))
    setObjects <- initSetFragInfos(setObjects, MSPeakListsSet)
    
    combFormulas <- Reduce(modifyList, lapply(setObjects, annotations, features = TRUE))
    
    cons <- makeFeatAnnSetConsensus(setObjects, names(fGroupsSet), setThreshold, setThresholdAnn, NULL)

    return(formulasSet(setObjects = setObjects, origFGNames = names(fGroupsSet), setThreshold = setThreshold,
                       setThresholdAnn = setThresholdAnn, groupAnnotations = cons, featureFormulas = combFormulas,
                       algorithm = makeSetAlgorithm(setObjects)))
}

#' @rdname formulas-class
#' @export
formulasUnset <- setClass("formulasUnset", contains = "formulas")

#' @rdname formulas-class
#' @export
setMethod("unset", "formulasSet", function(obj, set)
{
    assertSets(obj, set, FALSE)
    uann <- doFeatAnnUnset(obj, set)
    return(formulasUnset(groupAnnotations = uann, featureFormulas = annotations(obj, features = TRUE),
                         algorithm = paste0(algorithm(obj), "_unset")))
})

#' @rdname formulas-class
#' @export
setMethod("unset", "formulasConsensusSet", function(obj, set)
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
