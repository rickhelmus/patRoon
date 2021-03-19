#' @include main.R
#' @include formulas.R
#' @include workflow-step-set.R
NULL

makeFormulasSetConsensus <- function(setObjects, origFGNames, setThreshold, setThresholdAnn, mConsNames)
{
    groupFormsList <- sapply(setObjects, formulaTable, features = FALSE, simplify = FALSE)
    
    # prep set form tables
    groupFormsList <- Map(groupFormsList, names(setObjects), f = function(fTableG, set)
    {
        ft <- lapply(fTableG, function(ft)
        {
            ft <- copy(ft)
            rname <- paste0("rank-", set)
            uflen <- length(unique(ft$neutral_formula))
            ranks <- seq_len(uflen)
            ft[, (rname) := ranks[.GRP], by = "neutral_formula"]
            ft[, rankScore := (uflen - (get(rname) - 1)) / uflen] # for sorting
            
            # rename cols that are specific to a set or algo consensus or should otherwise not be combined
            # UNDONE: needed? can be inferred from set_from (although not so clear...)
            cols <- getAllMergedConsCols(c("rank", "mergedBy", "coverage"), names(ft), mConsNames)
            if (length(cols) > 0)
                setnames(ft, cols, paste0(cols, "-", set))
            
            return(ft)
        })
    })
    
    mc <- setNames(rep(length(groupFormsList), length(origFGNames)), origFGNames)
    ret <- generateGroupFormulasByConsensus(groupFormsList, mc, setThreshold, setThresholdAnn, origFGNames, "set_from",
                                            "sets", "setCoverage", "setCoverageAnn", NULL, NULL, mConsNames,
                                            "rankScore")
    
    # fix ranking
    
    ret <- lapply(ret, function(ft)
    {
        setorderv(ft, "rankScore", -1)
        # ft[, rankScore := NULL]
        return(ft)
    })
    
    return(ret[])
}

syncFormulasSetObjects <- function(formulasSet, makeCons)
{
    # update/initialize from setObjects
    if (length(setObjects(formulasSet)) >= 1)
    {
        if (length(setObjects(formulasSet)) > 1)
            formulasSet@featureFormulas <- Reduce(modifyList, lapply(formulasSet@setObjects, formulaTable, features = TRUE))
        else
            formulasSet@featureFormulas <- formulaTable(setObjects(formulasSet)[[1]], features = TRUE)
        
        if (makeCons)
            formulasSet@formulas <- makeFormulasSetConsensus(setObjects(formulasSet), formulasSet@setThreshold,
                                                             formulasSet@setThresholdAnn, formulasSet@origFGNames,
                                                             mergedConsensusNames(formulasSet))
        else
        {
            # sync available feature groups
            allFGroups <- unique(unlist(lapply(setObjects(formulasSet), groupNames)))
            formulasSet@formulas <- formulasSet@formulas[intersect(groupNames(formulasSet), allFGroups)]
            
            # only keep results from sets still present
            spat <- paste0(sets(formulasSet), collapse = "|")
            formulasSet@formulas <- lapply(formulasSet@formulas, function(ft) ft[grepl(spat, sets)])
        }
    }
    else
        formulasSet@featureFormulas <- formulasSet@formulas <- list()
    
    formulasSet@scoreRanges <- formulasSet@scoreRanges[groupNames(formulasSet)]
    
    return(formulasSet)
}

formulasSet <- setClass("formulasSet", slots = c(setThreshold = "numeric",
                                                 setThresholdAnn = "numeric",
                                                 origFGNames = "character"),
                        contains = c("formulas", "workflowStepSet"))

#' @describeIn formulasSet Shows summary information for this object.
#' @export
setMethod("show", "formulasSet", function(object)
{
    callAllNextMethods(object, show, firstClass = "formulas", startFrom = "formulasSet")
})

setMethod("[", c("formulasSet", "ANY", "missing", "missing"), function(x, i, j, ..., sets = NULL, drop = TRUE)
{
    assertSets(x, sets, TRUE)
    
    if (!is.null(sets))
        x@setObjects <- x@setObjects[sets]
    
    if (!missing(i))
    {
        i <- assertSubsetArgAndToChr(i, groupNames(x))
        # NOTE: assume that subsetting with non-existing i will not result in errors
        x@setObjects <- lapply(x@setObjects, "[", i = i)
    }

    if (!is.null(sets) || !missing(i))
        x <- syncFormulasSetObjects(x, FALSE)
    
    return(x)
})

setMethod("as.data.table", "formulasSet", function(x, fGroups = NULL, average = FALSE, ...)
{
    ret <- callNextMethod(x, fGroups = fGroups, average = average, ...)
    if (average)
        ret[, c("formula", "set_from") := NULL] # formula column doesn't make sense anymore, set_from column is left-over
    return(ret[])
})

#' @export
setMethod("filter", "formulasSet", function(obj, ..., negate = FALSE, sets = NULL)
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
        obj <- syncFormulasSetObjects(obj, TRUE)
        cat("Done!\n")
    }
    
    return(obj)
})

#' @export
setMethod("annotatedPeakList", "formulasSet", function(obj, precursor, groupName, analysis = NULL, MSPeakLists,
                                                       onlyAnnotated = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    assertChoiceSilent(groupName, groupNames(obj), add = ac)
    checkmate::assertString(analysis, min.chars = 1, null.ok = TRUE, add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakListsSet", add = ac)
    checkmate::reportAssertions(ac)
    
    ret <- rbindlist(Map(setObjects(obj), lapply(sets(obj), unset, obj = MSPeakLists), f = function(so, mspl)
    {
        if (!groupName %in% groupNames(so) || (!is.null(analysis) && !analysis %in% analyses(so)))
            return(NULL)
        annotatedPeakList(so, precursor, groupName, analysis, mspl, onlyAnnotated)
    }), fill = TRUE, idcol = "set")
    
    if (nrow(ret) == 0)
        return(NULL)
    
    return(ret)
})

#' @export
setMethod("plotSpectrum", "formulasSet", function(obj, precursor, groupName, analysis = NULL, MSPeakLists,
                                                  title = NULL, specSimParams = NULL, shift = "none",
                                                  useGGPlot2 = FALSE, mincex = 0.9, xlim = NULL,
                                                  ylim = NULL, perSet = TRUE, mirror = TRUE, ...)
{
    ac <- checkmate::makeAssertCollection()
    if (!is.null(specSimParams))
    {
        checkmate::assertCharacter(precursor, len = 2, min.chars = 1, add = ac)
        checkmate::assertCharacter(groupName, len = 2, min.chars = 1, add = ac)
        checkmate::assertCharacter(analysis, len = 2, min.chars = 1, null.ok = TRUE, add = ac)
        assertSpecSimParams(specSimParams, add = ac)
        checkmate::assertChoice(shift, c("none", "precursor", "both"), add = ac)
    }
    else
    {
        checkmate::assertString(precursor, min.chars = 1, add = ac)
        checkmate::assertString(groupName, min.chars = 1, add = ac)
        checkmate::assertString(analysis, min.chars = 1, null.ok = TRUE, add = ac)
    }
    checkmate::assertClass(MSPeakLists, "MSPeakListsSet", add = ac)
    checkmate::assertString(title, null.ok = TRUE, add = ac)
    checkmate::assertNumber(mincex, lower = 0, finite = TRUE, add = ac)
    assertXYLim(xlim, ylim, add = ac)
    aapply(checkmate::assertFlag, . ~ useGGPlot2 + perSet + mirror, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    if (!perSet || length(sets(obj)) == 1 || !is.null(analysis))
        return(callNextMethod(obj, precursor, groupName, analysis, MSPeakLists, title, specSimParams = NULL, shift,
                              useGGPlot2, mincex, xlim, ylim, ...))
    
    if (is.null(specSimParams))
    {
        if (is.null(title))
            title <- subscriptFormula(precursor)
        
        spec <- annotatedPeakList(obj, precursor, groupName, analysis, MSPeakLists)
        if (is.null(spec))
            return(NULL)
        
        specs <- split(spec, by = "set")
        specs <- lapply(specs, setnames, "set", "mergedBy")
        
        plotData <- getMSPlotDataOverlay(specs, mirror, TRUE, 1, NULL)
        return(makeMSPlotOverlay(plotData, title, mincex, xlim, ylim, useGGPlot2, ...))
    }
    else
    {
        if (is.null(title))
            title <- subscriptFormula(precursor[1], formulas2 = precursor[2])
        
        theSets <- sets(obj)
        
        usObj <- sapply(theSets, unset, obj = obj, simplify = FALSE)
        
        # check which sets actually contain requested data
        theSets <- theSets[sapply(theSets, function(s) all(groupName %in% groupNames(usObj[[s]])))]
        if (length(theSets) == 0)
            return(NULL)
        usObj <- usObj[theSets]
        
        usMSPL <- sapply(theSets, unset, obj = MSPeakLists, simplify = FALSE)
        binnedPLs <- Map(usMSPL, theSets, f = getBinnedPLPair,
                         MoreArgs = list(groupNames = groupName, analyses = NULL, MSLevel = 2,
                                         specSimParams = specSimParams, shift = shift, mustExist = FALSE))

        mergeBinnedAnn <- function(nr)
        {
            binPLs <- sapply(binnedPLs, "[[", nr, simplify = FALSE)
            annPLs <- Map(usObj, usMSPL, f = annotatedPeakList,
                          MoreArgs = list(precursor = precursor[nr], groupName = groupName[nr], analysis = analysis[nr]))
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

setMethod("plotSpectrumHash", "formulasSet", function(obj, precursor, groupName, analysis = NULL, MSPeakLists,
                                                      title = NULL, specSimParams = NULL, shift = "none",
                                                      useGGPlot2 = FALSE, mincex = 0.9, xlim = NULL, ylim = NULL,
                                                      perSet = TRUE, mirror = TRUE, ...)
{
    return(makeHash(callNextMethod(obj, precursor, groupName, analysis, MSPeakLists,
                                   title, specSimParams, shift, useGGPlot2, mincex, xlim, ylim, ...),
                    perSet, mirror))
})

#' @export
setMethod("consensus", "formulasSet", function(obj, ..., absMinAbundance = NULL,
                                               relMinAbundance = NULL,
                                               uniqueFrom = NULL, uniqueOuter = FALSE,
                                               rankWeights = 1, labels = NULL, setThreshold = 0,
                                               setThresholdAnn = 0.75)
{
    # make consensus of shared setObjects
    # add unique setObjects
    # make 'regular' set consensus from new results
    
    allFormulas <- c(list(obj), list(...))
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertList(allFormulas, types = "formulasSet", min.len = 2, any.missing = FALSE,
                          unique = TRUE, .var.name = "...", add = ac)
    checkmate::assertNumeric(rankWeights, lower = 0, finite = TRUE, add = ac)
    checkmate::assertCharacter(labels, min.chars = 1, len = length(allFormulas), null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (!allSame(lapply(allFormulas, sets)))
        stop("All objects must have the same sets.")
    
    # NOTE: don't want to keep -set suffix
    formNames <- if (!is.null(labels)) labels else sub("\\-set$", "", sapply(allFormulas, algorithm))
    
    setObjects <- sapply(sets(obj), function(set)
    {
        return(do.call(consensus, c(lapply(lapply(allFormulas, setObjects), "[[", set),
                                    list(absMinAbundance = absMinAbundance, relMinAbundance = relMinAbundance,
                                         uniqueFrom = uniqueFrom, uniqueOuter = uniqueOuter,
                                         rankWeights = rankWeights, labels = formNames))))
    }, simplify = FALSE)

    combFormulas <- Reduce(modifyList, lapply(setObjects, formulaTable, features = TRUE))
    
    gNames <- allFormulas[[1]]@origFGNames # UNDONE? at least verify all objects are equal
    groupForms <- makeFormulasSetConsensus(setObjects, gNames, setThreshold, setThresholdAnn, formNames)
    
    browser()
    return(formulasSet(setObjects = setObjects, origFGNames = names(fGroupsSet), setThreshold = setThreshold,
                       setThresholdAnn = setThresholdAnn, formulas = groupForms, featureFormulas = combFormulas,
                       algorithm = makeSetAlgorithm(setObjects)))
    
    cons <- makeCompoundsSetConsensus(setObjects, obj@origFGNames, setThreshold, setThresholdAnn, formNames)
    sc <- makeCompoundsSetScorings(setObjects, obj@origFGNames)
    
    return(compoundsConsensusSet(setObjects = setObjects, setThreshold = setThreshold, setThresholdAnn = setThresholdAnn,
                                 origFGNames = obj@origFGNames, compounds = cons, scoreTypes = sc$scTypes,
                                 scoreRanges = sc$scRanges,
                                 algorithm = paste0(unique(sapply(allFormulas, algorithm)), collapse = ","),
                                 mergedConsensusNames = formNames))
})


generateFormulasSet <- function(fGroupsSet, generator, ..., setArgs, setThreshold, setThresholdAnn)
{
    aapply(checkmate::assertNumber, . ~ setThreshold + setThresholdAnn, lower = 0, upper = 1, finite = TRUE)
    
    # UNDONE: mention that adduct argument is automatically set

    unsetFGroupsList <- sapply(sets(fGroupsSet), unset, obj = fGroupsSet, simplify = FALSE)
    setObjects <- Map(unsetFGroupsList, setArgs,
                      f = function(fg, sa) do.call(generator, c(list(fGroups = fg, ...), sa)))
    
    combFormulas <- Reduce(modifyList, lapply(setObjects, formulaTable, features = TRUE))
    
    groupForms <- makeFormulasSetConsensus(setObjects, names(fGroupsSet), setThreshold, setThresholdAnn, character())

    return(formulasSet(setObjects = setObjects, origFGNames = names(fGroupsSet), setThreshold = setThreshold,
                       setThresholdAnn = setThresholdAnn, formulas = groupForms, featureFormulas = combFormulas,
                       algorithm = makeSetAlgorithm(setObjects)))
}

formulasUnset <- setClass("formulasUnset", contains = "formulas")
setMethod("unset", "formulasSet", function(obj, set)
{
    assertSets(obj, set, FALSE)
    obj <- obj[, sets = set]
    
    groupForms <- lapply(formulaTable(obj), copy)
    groupForms <- lapply(groupForms, data.table::set, j = c("sets", "setCoverage"), value = NULL)
    
    return(formulasUnset(formulas = groupForms, featureFormulas = formulaTable(obj, features = TRUE),
                         algorithm = paste0(algorithm(obj), "_unset")))
})
