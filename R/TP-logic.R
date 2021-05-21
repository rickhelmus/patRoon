#' @include main.R
#' @include TP.R
#' @include feature_groups-set.R
NULL

getTPLogicTransformations <- function(transformations)
{
    if (is.null(transformations))
        transformations <- patRoon:::TPsLogicTransformations # stored inside R/sysdata.rda
    
    if (is.data.table(transformations))
        ret <- copy(transformations)
    else
        ret <- as.data.table(transformations)
    
    ret[, deltaMZ := 0]
    ret[nzchar(add), deltaMZ := sapply(add, getFormulaMass)]
    ret[nzchar(sub), deltaMZ := deltaMZ - sapply(sub, getFormulaMass)]
    
    return(ret[])
}

doGenerateTPsLogic <- function(fGroups, minMass, neutralMasses, transformations)
{
    if (length(fGroups) == 0)
        return(list(parents = data.table(name = character(), rt = numeric(), neutralMass = numeric()),
                    products = list()))
    
    gInfo <- groupInfo(fGroups)
    parents <- data.table(name = names(fGroups), rt = gInfo$rts, neutralMass = neutralMasses)
    
    transformations <- getTPLogicTransformations(transformations)
    
    prog <- openProgBar(0, nrow(parents))
    
    products <- lapply(seq_len(nrow(parents)), function(si)
    {
        ret <- data.table(name = paste0(parents$name[si], "-", transformations$transformation),
                          neutralMass = neutralMasses[si] + transformations$deltaMZ,
                          deltaMZ = transformations$deltaMZ,
                          trans_add = transformations$add,
                          trans_sub = transformations$sub,
                          retDir = transformations$retDir)
        ret <- ret[neutralMass >= minMass]
        
        setTxtProgressBar(prog, si)
        
        return(ret)
    })
    names(products) <- parents$name
    
    setTxtProgressBar(prog, nrow(parents))
    close(prog)
    
    return(list(parents = parents, products = products))
}

#' @export
transformationProductsLogic <- setClass("transformationProductsLogic", contains = "transformationProducts")

setMethod("initialize", "transformationProductsLogic",
          function(.Object, ...) callNextMethod(.Object, algorithm = "logic", ...))

setMethod("needsScreening", "transformationProductsComponents", function(TPs) FALSE)

setMethod("linkParentsToFGroups", "transformationProductsLogic", function(TPs, fGroups)
{
    fg <- intersect(names(TPs), names(fGroups))
    return(data.table(name = fg, group = fg))
})

#' @export
setMethod("generateTPsLogic", "featureGroups", function(fGroups, minMass = 40, adduct = NULL, transformations = NULL)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertNumber(minMass, lower = 0, finite = TRUE, add = ac)
    assertLogicTransformations(transformations, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    adduct <- checkAndToAdduct(adduct, fGroups)
    fgAdd <- getFGroupAdducts(names(fGroups), annotations(fGroups), adduct, "generic")
    neutralMasses <- if (length(fGroups) == 0) numeric () else groupInfo(fGroups)[, "mzs"] - sapply(fgAdd$grpAdducts,
                                                                                                    adductMZDelta)
    
    res <- doGenerateTPsLogic(fGroups, minMass, neutralMasses, transformations)
    
    return(transformationProductsLogic(parents = res$parents, products = res$products))
})

#' @export
setMethod("generateTPsLogic", "featureGroupsSet", function(fGroups, minMass = 40, transformations = NULL)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertNumber(minMass, lower = 0, finite = TRUE, add = ac)
    assertLogicTransformations(transformations, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    res <- doGenerateTPsLogic(fGroups, minMass, groupInfo(fGroups)[, "mzs"], transformations)
    
    return(transformationProductsLogic(parents = res$parents, products = res$products))
})
