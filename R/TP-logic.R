#' @include main.R
#' @include TP.R
#' @include feature_groups-screening-set.R
NULL

# UNDONE: precursor --> parent?

getTPLogicTransformations <- function(transformations)
{
    if (is.null(transformations))
        ret <- patRoon:::TPsLogicReactions # stored inside R/sysdata.rda
    
    if (is.data.table(ret))
        ret <- copy(ret)
    else
        ret <- as.data.table(ret)
    
    ret[, deltaMZ := 0]
    ret[nzchar(add), deltaMZ := sapply(add, getFormulaMass)]
    ret[nzchar(sub), deltaMZ := deltaMZ - sapply(sub, getFormulaMass)]
    
    return(ret[])
}

doPredictTPsLogic <- function(fGroups, minMass, neutralMasses, transformations)
{
    gInfo <- groupInfo(fGroups)
    suspects <- data.table(name = names(fGroups), rt = gInfo$rts, neutralMass = neutralMasses)
    
    transformations <- getTPLogicTransformations(transformations)
    
    prog <- openProgBar(0, nrow(suspects))
    
    predictions <- lapply(seq_len(nrow(suspects)), function(si)
    {
        ret <- data.table(name = paste0(suspects$name[si], "-", transformations$reaction),
                          neutralMass = neutralMasses[si] + transformations$deltaMZ,
                          deltaMZ = transformations$deltaMZ,
                          reaction_add = transformations$add,
                          reaction_sub = transformations$sub,
                          RTDir = transformations$RTDir)
        ret <- ret[neutralMass >= minMass]
        
        setTxtProgressBar(prog, si)
        
        return(ret)
    })
    names(predictions) <- suspects$name
    
    setTxtProgressBar(prog, nrow(suspects))
    close(prog)
    
    return(list(suspects = suspects, predictions = predictions))
}

#' @export
TPPredictionsLogic <- setClass("TPPredictionsLogic", contains = "TPPredictions")

setMethod("initialize", "TPPredictionsLogic",
          function(.Object, ...) callNextMethod(.Object, algorithm = "logic", ...))

setMethod("linkPrecursorsToFGroups", "TPPredictionsLogic", function(pred, fGroups)
{
    fg <- intersect(names(pred), names(fGroups))
    return(data.table(name = fg, group = fg))
})

#' @export
setMethod("predictTPsLogic", "featureGroups", function(fGroups, minMass = 40, adduct = NULL, transformations = NULL)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertNumber(minMass, finite = TRUE, add = ac)
    assertLogicTransformations(transformations, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    adduct <- checkAndToAdduct(adduct, fGroups)
    fgAdd <- getFGroupAdducts(names(fGroups), annotations(fGroups), adduct, "generic")
    neutralMasses <- groupInfo(fGroups)[, "mzs"] - sapply(fgAdd$grpAdducts, adductMZDelta)
    
    res <- doPredictTPsLogic(fGroups, minMass, neutralMasses, transformations)
    
    return(TPPredictionsLogic(suspects = res$suspects, predictions = res$predictions))
})

#' @export
setMethod("predictTPsLogic", "featureGroupsSet", function(fGroups, minMass = 40, transformations = NULL)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertNumber(minMass, finite = TRUE, add = ac)
    assertLogicTransformations(transformations, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    res <- doPredictTPsLogic(fGroups, minMass, groupInfo(fGroups)[, "mzs"], transformations)
    
    return(TPPredictionsLogic(suspects = res$suspects, predictions = res$predictions))
})
