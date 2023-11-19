#' @include main.R
#' @include mspeaklists.R
NULL

syncMSPeakListsSetObjects <- function(MSPeakListsSet)
{
    # update/initialize from setObjects
    if (length(setObjects(MSPeakListsSet)) >= 2)
    {
        MSPeakListsSet@peakLists <- Reduce(modifyList, lapply(MSPeakListsSet@setObjects, peakLists))
        MSPeakListsSet@averagedPeakLists <- averageMSPeakLists(MSPeakListsSet)
    }
    else if (length(setObjects(MSPeakListsSet)) == 1)
    {
        MSPeakListsSet@peakLists <- peakLists(MSPeakListsSet@setObjects[[1]])
        MSPeakListsSet@averagedPeakLists <- averagedPeakLists(MSPeakListsSet@setObjects[[1]])
    }
    else
        MSPeakListsSet@peakLists <- MSPeakListsSet@averagedPeakLists <- list()
    
    MSPeakListsSet@analysisInfo <- MSPeakListsSet@analysisInfo[analysis %in% names(peakLists(MSPeakListsSet))]
    
    return(MSPeakListsSet)
}

#' @param set \setsWF The name of the set.
#' @param sets \setsWF A \code{character} with name(s) of the sets to keep (or remove if \code{negate=TRUE}).
#'
#' @slot analysisInfo \setsWF  \link[=analysis-information]{Analysis information}. Use the \code{analysisInfo} method
#'   for access.
#'
#' @section Sets workflows: \setsWFClass{MSPeakListsSet}{MSPeakLists}
#'
#'   \setsWFNewMethodsSOExtra{MSPeakListsUnset}{Only the MS peaks that are present in the specified set are kept.}{
#'
#'   \item \code{analysisInfo} Returns the  \link[=analysis-information]{analysis info} for this object.
#'
#'   }
#'
#'   \setsWFChangedMethods{
#'
#'   \item \code{filter} and the subset operator (\code{[}) Can be used to select data that is only present for selected
#'   sets. The \code{filter} method is applied for each set individually, and afterwards the results are combined again
#'   (see \code{\link{generateMSPeakLists}}). Note that this has important implications for \emph{e.g.} relative
#'   intensity filters (\code{relMSIntThr}/\code{relMSMSIntThr}), \code{topMSPeaks}/\code{topMSMSPeaks} and
#'   \code{minMSMSPeaks}. Similarly, when the \code{annotatedBy} filter is applied, each set specific MS peak list is
#'   filtered by the annotation results from only that set.
#'
#'   \item \code{plotSpectrum} Is able to highlight set specific mass peaks (\code{perSet} and \code{mirror} arguments).
#'
#'   \item \code{spectrumSimilarity} First calculates similarities for each spectral pair per set (\emph{e.g.} all
#'   positive mode spectra are compared and then all negative mode spectra are compared). This data is then combined
#'   into an overall similarity value. How this combination is performed depends on the \code{setCombineMethod} field of
#'   the \code{\link{specSimParams}} argument.
#'
#'   }
#'
#' @rdname MSPeakLists-class
#' @export
MSPeakListsSet <- setClass("MSPeakListsSet",
                           slots = c(analysisInfo = "data.table"),
                           contains = c("MSPeakLists", "workflowStepSet"))

setMethod("initialize", "MSPeakListsSet", function(.Object, ...) callNextMethod(.Object, ..., setIDs = FALSE))
    
setMethod("averageMSPeakLists", "MSPeakListsSet", function(obj)
{
    # create 'averaged' peak lists by simply merging the averaged lists from the setObjects
    
    hash <- makeHash(lapply(obj@setObjects, averagedPeakLists))
    avgPLists <- loadCacheData("MSPeakListsSetAvg", hash)
    
    gNames <- unique(unlist(sapply(obj@setObjects, groupNames, simplify = FALSE), use.names = FALSE))
    gNames <- intersect(obj@origFGNames, gNames) # sort to original order
    
    if (length(gNames) == 0)
        avgPLists <- makeEmptyListNamed(list())
    else if (is.null(avgPLists))
    {
        avgPLists <- sapply(gNames, function(gName)
        {
            PLMS <- rbindlist(pruneList(lapply(obj@setObjects, function(mspl) mspl[[gName]][["MS"]])), idcol = "set")
            PLMSMS <- rbindlist(pruneList(lapply(obj@setObjects, function(mspl) mspl[[gName]][["MSMS"]])), idcol = "set")
            
            return(pruneList(list(MS = if (nrow(PLMS) > 0) PLMS else NULL,
                                  MSMS = if (nrow(PLMSMS) > 0) PLMSMS else NULL)))
        }, simplify = FALSE)
        avgPLists <- pruneList(avgPLists, checkEmptyElements = TRUE)
        saveCacheData("MSPeakListsSetAvg", avgPLists, hash)
    }

    return(avgPLists)
})


#' @rdname MSPeakLists-class
#' @export
setMethod("analysisInfo", "MSPeakListsSet", function(obj, df)
{
    checkmate::assertFlag(df)
    return(if (df) as.data.frame(obj@analysisInfo) else obj@analysisInfo)
})

#' @rdname MSPeakLists-class
#' @export
setMethod("show", "MSPeakListsSet", function(object)
{
    callAllNextMethods(object, show, firstClass = "MSPeakLists", startFrom = "MSPeakListsSet")
})

#' @rdname MSPeakLists-class
#' @export
setMethod("[", c("MSPeakListsSet", "ANY", "ANY", "missing"), function(x, i, j, ..., reAverage = FALSE,
                                                                      sets = NULL, drop = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    assertSets(x, sets, TRUE, add = ac)
    checkmate::assertFlag(reAverage, add = ac)
    checkmate::reportAssertions(ac)

    if (!is.null(sets))
        i <- mergeAnaSubsetArgWithSets(i, sets, analysisInfo(x))
    
    if (!missing(i) || !missing(j))
    {
        args <- list(reAverage = reAverage)
        if (!missing(i))
            args <- c(args, list(i = assertSubsetArgAndToChr(i, analyses(x))))
        if (!missing(j))
            args <- c(args, list(j = assertSubsetArgAndToChr(j, groupNames(x))))
        
        if (!is.null(sets))
            x@setObjects <- x@setObjects[sets]
        
        # NOTE: assume that subsetting with non-existing i/j will not result in errors
        x@setObjects <- lapply(x@setObjects, function(o) do.call("[", args = c(list(x = o), args)))
        
        x <- syncMSPeakListsSetObjects(x)
    }
    
    return(x)
})

#' @rdname MSPeakLists-class
#' @export
setMethod("as.data.table", "MSPeakListsSet", function(x, fGroups = NULL, averaged = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroupsSet", null.ok = TRUE, add = ac)
    checkmate::assertFlag(averaged, add = ac)
    checkmate::reportAssertions(ac)
    
    ret <- callNextMethod(x, fGroups = fGroups, averaged = averaged)
    
    if (!averaged) # add set column
    {
        anaInfo <- analysisInfo(x)
        ret[, set := anaInfo$set[match(analysis, anaInfo$analysis)]]
        setcolorder(ret, "set")
    }
    
    return(ret[])
})

#' @rdname MSPeakLists-class
#' @export
setMethod("delete", "MSPeakListsSet", function(obj, ...)
{
    obj <- callNextMethod()
    
    # sync set objects
    obj@setObjects <- Map(sets(obj), obj@setObjects, f = function(sn, so)
    {
        # remove removed analyses
        so <- delete(so, k = !analyses(so) %chin% analyses(obj))
        
        # remove removed groups
        so <- delete(so, i = !groupNames(so) %chin% groupNames(obj))
        
        # remove removed peaks
        delete(so, j = function(PL, grp, ana, t)
        {
            setPL <- if (is.null(ana)) obj[[grp]][[t]] else obj[[ana, grp]][[t]]
            if (is.null(setPL))
                return(TRUE)
            if (is.null(ana))
                setPL <- setPL[set == sn]
            return(!PL$ID %in% setPL$ID)
        })
    })
    
    return(obj)
})

#' @rdname MSPeakLists-class
#' @export
setMethod("filter", "MSPeakListsSet", function(obj, ..., annotatedBy = NULL, retainPrecursorMSMS = TRUE,
                                               reAverage = FALSE, negate = FALSE, sets = NULL)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFlag(negate, add = ac)
    checkmate::assert(
        checkmate::checkNull(annotatedBy),
        checkmate::checkClass(annotatedBy, "formulasSet"),
        checkmate::checkClass(annotatedBy, "compoundsSet"),
        checkmate::checkList(annotatedBy, c("formulasSet", "compoundsSet"), any.missing = FALSE, min.len = 1,
                             unique = TRUE),
        .var.name = "annotatedBy"
    )
    assertSets(obj, sets, TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    annotatedByList <- NULL
    if (!is.null(annotatedBy))
    {
        # unset objects before passing them to parent method
        
        if (!is.list(annotatedBy))
            annotatedByList <- lapply(sets(obj), unset, obj = annotatedBy)
        else
            annotatedByList <- lapply(sets(obj), function(s) lapply(annotatedBy, unset, set = s))
    }
    
    if (!is.null(sets) && length(sets) > 0)
    {
        if (negate)
            sets <- setdiff(get("sets", pos = 2)(obj), sets)
        obj <- obj[, sets = sets]
    }
    
    if (...length() > 0 || !is.null(annotatedBy))
    {
        if (is.null(annotatedBy))
            obj@setObjects <- lapply(obj@setObjects, filter, ..., retainPrecursorMSMS = retainPrecursorMSMS,
                                     reAverage = reAverage, negate = negate)
        else
        {
            obj@setObjects <- Map(obj@setObjects, annotatedByList, f = function(so, ab)
            {
                filter(so, ..., annotatedBy = ab, retainPrecursorMSMS = retainPrecursorMSMS,
                       reAverage = reAverage, negate = negate)
            })
        }
        
        # synchronize other objects
        cat("Synchronizing set objects...\n")
        obj <- syncMSPeakListsSetObjects(obj)
        cat("Done!\n")
    }
    
    return(obj)
})

#' @param perSet,mirror \setsWF If \code{perSet=TRUE} then the set specific mass peaks are annotated separately.
#'   Furthermore, if \code{mirror=TRUE} (and there are two sets in the object) then a mirror plot is generated.
#' @rdname MSPeakLists-class
#' @export
setMethod("plotSpectrum", "MSPeakListsSet", function(obj, groupName, analysis = NULL, MSLevel = 1, title = NULL,
                                                     specSimParams = getDefSpecSimParams(),
                                                     xlim = NULL, ylim = NULL, perSet = TRUE, mirror = TRUE, ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertCharacter(groupName, min.len = 1, max.len = 2, min.chars = 1, add = ac)
    checkmate::assertCharacter(analysis, min.len = 1, max.len = 2, min.chars = 1, null.ok = TRUE, add = ac)
    if (!is.null(analysis) && length(analysis) != length(groupName))
        stop("Lengths of analysis and groupName should be equal.")
    assertSpecSimParams(specSimParams, add = ac)
    checkmate::assertChoice(MSLevel, 1:2, add = ac)
    assertXYLim(xlim, ylim, add = ac)
    aapply(checkmate::assertFlag, . ~ perSet, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    argsParent <- list(groupName = groupName, analysis = analysis, MSLevel = MSLevel, title = title,
                       specSimParams = specSimParams, xlim = xlim, ylim = ylim, ...)
    
    if (!perSet || length(sets(obj)) == 1 || !is.null(analysis))
        return(do.call(callNextMethod, c(list(obj), argsParent)))
    
    setTitle <- is.null(title)
    if (setTitle)
        title <- getMSPeakListPlotTitle(MSLevel, analysis, groupName)
    
    if (length(groupName) == 1)
    {
        spec <- getSpec(obj, groupName, MSLevel, analysis)
        if (is.null(spec))
            return(NULL)
        
        specs <- split(spec, by = "set")
        specs <- lapply(specs, setnames, "set", "mergedBy")
        
        plotData <- getMSPlotDataOverlay(specs, mirror, TRUE, 1, NULL)
        return(makeMSPlotOverlay(plotData, title, 1, xlim, ylim, ...))
    }
    else
    {
        if (setTitle)
        {
            sim <- spectrumSimilarity(obj, groupName[1], groupName[2], analysis[1], analysis[2],
                                      MSLevel, specSimParams, NAToZero = TRUE, drop = TRUE)
            title <- c(title, sprintf("Similarity: %.2f", sim))
        }
        
        theSets <- sets(obj)
        usObj <- sapply(theSets, unset, obj = obj, simplify = FALSE)
        
        binnedPLs <- Map(usObj, theSets, f = getBinnedPLPair,
                         MoreArgs = list(groupNames = groupName, analyses = analysis, MSLevel = MSLevel,
                                         specSimParams = specSimParams, mustExist = FALSE))
        if (all(unlist(lapply(binnedPLs, sapply, nrow)) == 0))
            return(NULL)

        topSpec <- rbindlist(sapply(binnedPLs, "[[", 1, simplify = FALSE), idcol = "set")
        bottomSpec <- rbindlist(sapply(binnedPLs, "[[", 2, simplify = FALSE), idcol = "set")
        
        plotData <- getMSPlotDataOverlay(list(topSpec, bottomSpec), mirror, FALSE, 2, "overlap")
        makeMSPlotOverlay(plotData, title, 1, xlim, ylim, ...)
    }
})

setMethod("plotSpectrumHash", "MSPeakListsSet", function(obj, groupName, analysis = NULL, MSLevel = 1, title = NULL,
                                                         specSimParams = getDefSpecSimParams(), xlim = NULL,
                                                         ylim = NULL, perSet = TRUE, mirror = TRUE, ...)
{
    return(makeHash(callNextMethod(obj, groupName, analysis, MSLevel, title, specSimParams, xlim, ylim, ...),
                    perSet, mirror))
})

#' @rdname MSPeakLists-class
#' @export
setMethod("spectrumSimilarity", "MSPeakListsSet", function(obj, groupName1, groupName2 = NULL, analysis1 = NULL,
                                                           analysis2 = NULL, MSLevel = 1,
                                                           specSimParams = getDefSpecSimParams(), NAToZero = FALSE,
                                                           drop = TRUE)
{
    if (length(obj) == 0)
        return(NULL)
    
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertSubset, . ~ groupName1 + groupName2, empty.ok = c(FALSE, TRUE),
           fixed = list(choices = groupNames(obj), add = ac))
    aapply(checkmate::assertFlag, . ~ NAToZero + drop, fixed = list(add = ac))
    checkmate::reportAssertions(ac)

    sims <- pruneList(lapply(setObjects(obj), function(so)
    {
        wh <- groupName1 %in% groupNames(so)
        if (!is.null(analysis1))
            wh <- wh & (analysis1 %in% analyses(so))
        gn1 <- groupName1[wh]
        ana1 <- if (is.null(analysis1)) NULL else analysis1[wh]
        
        if (is.null(groupName2))
            gn2 <- ana2 <- NULL
        else
        {
            wh <- groupName2 %in% groupNames(so)
            if (!is.null(analysis2))
                wh <- wh & (analysis2 %in% analyses(so))
            gn2 <- groupName2[wh]
            ana2 <- if (is.null(analysis2)) NULL else analysis2[wh]
        }
        
        if (length(gn1) == 0 || (!is.null(groupName2) && length(gn2) == 0))
            ret <- as.matrix(numeric()) # will be expanded below
        else
        {
            # NOTE: don't drop NAs/dimensions here yet
            ret <- spectrumSimilarity(so, gn1, gn2, ana1, ana2, MSLevel, specSimParams, NAToZero = FALSE, drop = FALSE)
        }
        return(expandFillSpecSimilarities(ret, groupName1, if (is.null(groupName2)) groupName1 else groupName2))
    }))

    if (length(sims) == 0)
        return(NA_real_)
        
    if (length(sims) > 1)
    {
        if (specSimParams$setCombineMethod == "mean")
        {
            # deal with NAs
            simsNONA <- lapply(sims, function(s) { s[is.na(s)] <- 0; return(s) })
            noNACounts <- lapply(sims, function(s) matrix(!is.na(s), nrow(s), ncol(s)))
            noNACounts <- Reduce("+", noNACounts)
            
            sims <- Reduce("+", simsNONA) / noNACounts
            sims[is.nan(sims)] <- NA # may be introduced for noNACounts==0
        }
        else
            sims <- do.call(if (specSimParams$setCombineMethod == "min") pmin else pmax, c(sims, list(na.rm = TRUE)))
    }
    else
        sims <- sims[[1]]
    
    if (NAToZero)
        sims[is.na(sims)] <- 0
    
    return(if (drop && length(sims) == 1) drop(sims) else sims)
})


generateMSPeakListsSet <- function(fGroupsSet, generator, ...)
{
    # unset all fGroups sets, calculate MS peak lists for each set and store in setObjects
    # store combined setObject results in peakLists
    # store merged averaged peak lists in averagedPeakLists
    
    unsetFGroupsList <- sapply(sets(fGroupsSet), unset, obj = fGroupsSet, simplify = FALSE)
    setObjects <- sapply(unsetFGroupsList, generator, ..., simplify = FALSE)
    
    # combine non averaged (per analysis) MSPeakLists
    combPL <- Reduce(modifyList, lapply(setObjects, peakLists))

    # UNDONE: set metadata?
    ret <- MSPeakListsSet(setObjects = setObjects,
                          analysisInfo = copy(analysisInfo(fGroupsSet)),
                          peakLists = combPL, metadata = list(),
                          origFGNames = names(fGroupsSet),
                          algorithm = makeSetAlgorithm(setObjects))
    
    return(ret)
}

#' @rdname MSPeakLists-class
#' @export
MSPeakListsUnset <- setClass("MSPeakListsUnset", contains = "MSPeakLists")

#' @rdname MSPeakLists-class
#' @export
setMethod("unset", "MSPeakListsSet", function(obj, set)
{
    assertSets(obj, set, FALSE)
    obj <- obj[, sets = set]

    hasSO <- length(obj@setObjects) > 0    
    avArgs <- if (hasSO) obj@setObjects[[1]]@avgPeakListArgs else list()
    
    # # only re-average if >1 sets, otherwise just copy the setObject
    # UNDONE: re-enable if unset will support >1 sets again...
    # if (length(obj@setObjects) > 1)
    #     return(MSPeakListsUnset(peakLists = obj@peakLists, metadata = list(), avgPeakListArgs = avArgs,
    #                             origFGNames = obj@origFGNames, algorithm = paste0(algorithm(obj), "_unset")))
    
    return(MSPeakListsUnset(peakLists = if (hasSO) obj@setObjects[[1]]@peakLists else list(),
                            averagedPeakLists = if (hasSO) averagedPeakLists(obj@setObjects[[1]]) else list(),
                            metadata = list(), avgPeakListArgs = avArgs, setIDs = FALSE, doAverage = FALSE,
                            origFGNames = obj@origFGNames, algorithm = paste0(algorithm(obj), "_unset")))
})

