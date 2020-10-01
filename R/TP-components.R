#' @include main.R
#' @include TP.R
#' @include components.R
#' @include components-set.R
NULL

#' @export
TPPredictionsComponents <- setClass("TPPredictionsComponents", contains = "TPPredictions")

setMethod("initialize", "TPPredictionsComponents",
          function(.Object, ...) callNextMethod(.Object, algorithm = "components", ...))

setMethod("needsScreening", "TPPredictionsComponents", function(pred) FALSE)

setMethod("linkPrecursorsToFGroups", "TPPredictionsComponents", function(pred, fGroups)
{
    fg <- intersect(names(pred), names(fGroups))
    return(data.table(name = fg, group = fg))
})

setMethod("linkTPsToFGroups", "TPPredictionsComponents", function(pred, fGroups)
{
    TPNames <- intersect(as.data.table(pred)$name, names(fGroups))
    return(data.table(group = TPNames, TP_name = TPNames))
})

#' @export
setMethod("predictTPsComponents", "components", function(components, fGroupsPrecursors, fGroupsTPs)
{
    # UNDONE: mention this only works for components where a precursor
    # occurs not more than once in a component
    # UNDONE: cache
    
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertClass, . ~ fGroupsPrecursors + fGroupsTPs, "featureGroups", fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    gInfo <- groupInfo(fGroupsPrecursors)
    suspects <- data.table(name = rownames(gInfo), rt = gInfo$rts, mz = gInfo$mzs)
    cTab <- as.data.table(components)
    cTabTPs <- as.data.table(components[, groupNames(fGroupsTPs)])
    
    prog <- openProgBar(0, nrow(suspects))
    
    predictions <- pruneList(setNames(lapply(seq_len(nrow(suspects)), function(si)
    {
        precGName <- suspects$name[si]
        precCMP <- cTab[group == precGName]$name
        if (length(precCMP) == 0)
            return(NULL) # not found
        if (length(precCMP) > 1)
        {
            warning(sprintf("Found multiple components with precursor %s! Only first will be considered.", precGName))
            precCMP <- precCMP[1]
        }
        
        ret <- cTabTPs[name == precCMP]
        if (nrow(ret) == 0)
            return(NULL)
        
        ret <- ret[, c("group", "ret", "mz"), with = FALSE]
        setnames(ret, "group", "name")
        ret[, RTDir := 0]
        
        # UNDONE: more checks?
        
        setTxtProgressBar(prog, si)
        
        return(ret)
    }), suspects$name))
    
    setTxtProgressBar(prog, nrow(suspects))
    close(prog)
    
    suspects <- suspects[name %in% names(predictions)]
    
    return(TPPredictionsComponents(suspects = suspects, predictions = predictions))
})

#' @export
setMethod("predictTPsComponents", "componentsSet", function(components, fGroupsPrecursors, fGroupsTPs)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertList, . ~ fGroupsPrecursors + fGroupsTPs,
           types = "featureGroupsSet", min.len = 1, names = "unique", fixed = list(add = ac))
    checkmate::reportAssertions(ac)

    if (length(fGroupsPrecursors) != length(fGroupsTPs))
        stop("Different amount of precursor/TP featureGroups objects specified")
    if (any(names(fGroupsPrecursors) != names(fGroupsTPs)))
        stop("Different sets for precursor/TP featureGroups objects specified")
    if (!all(names(fGroupsPrecursors) %in% sets(components)))
        stop("Not all sets have components")
    
    theSets <- names(fGroupsPrecursors)
    
    ionizedFGroupsPrec <- Map(ionize, fGroupsPrecursors, names(fGroupsPrecursors))
    ionizedFGroupsTPs <- Map(ionize, fGroupsTPs, names(fGroupsTPs))
    
    # first calculate for each set
    predictSets <- Map(predictTPsComponents, setObjects(components)[theSets], ionizedFGroupsPrec, ionizedFGroupsTPs)
    TPListSets <- Map(predictSets, theSets, f = function(p, s)
    {
        ret <- lapply(predictions(p), function(tps)
        {
            tps <- copy(tps)
            tps[, (theSets) := as.list(theSets == s)]
            return(tps)
        })
        return(ret)
    })
    
    pred <- ReduceWithArgs(x = TPListSets, theSets, f = function(left, right, setLeft, setRight)
    {
        ret <- left
        over <- intersect(names(left), names(right))
        ret[over] <- Map(ret[over], right[over], f = function(TPsL, TPsR)
        {
            TPsL[name %in% TPsR$name, (setRight) := TRUE] # mark overlap
            TPsL <- rbind(TPsL, TPsR[!name %in% TPsL$name], fill = TRUE)
            return(TPsL)
        })
        
        return(c(ret, right[setdiff(names(right), names(left))]))
    })
    
    # convert to data.tables with group name, as data.frame rbind messes up rownames
    allGInfos <- rbindlist(lapply(fGroupsPrecursors, function(fg)
    {
        gi <- as.data.table(groupInfo(fg))
        gi[, group := names(fg)]
    }))
    suspects <- data.table(name = names(pred))
    suspects[, c("rt", "mz") := allGInfos[match(name, group), c("rts", "mzs")]]
    for (s in theSets)
        suspects[, (s) := name %in% names(fGroupsPrecursors[[s]])][]
    
    return(TPPredictionsComponents(suspects = suspects, predictions = pred))
})
