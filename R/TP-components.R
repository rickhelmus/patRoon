#' @include main.R
#' @include TP.R
NULL

#' @export
TPPredictionsComponents <- setClass("TPPredictionsComponents", contains = "TPPredictions")

setMethod("initialize", "TPPredictionsComponents",
          function(.Object, ...) callNextMethod(.Object, algorithm = "logic", ...))

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
predictTPsComponents <- function(fGroupsPrecursors, fGroupsTPs, components)
{
    # UNDONE: mention this only works for components where a precursor
    # occurs not more than once in a component
    # UNDONE: cache
    
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertClass, . ~ fGroupsPrecursors + fGroupsTPs + components,
           c("featureGroups", "featureGroups", "components"), fixed = list(add = ac))
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
}
