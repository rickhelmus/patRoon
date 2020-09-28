#' @include main.R
#' @include TP.R
#' @include feature_groups-screening-set.R
NULL

# UNDONE: precursor --> parent?

getTPLogicTransformations <- function()
{
    # UNDONE: cite 10.1021/acs.analchem.5b02905
    
    ret <- rbindlist(list(
        list("hydroxylation", "O", "", -1),
        list("demethylation", "", "CH2", -1),
        list("deethylation", "", "C2H4", -1),
        list("dehydrogenation", "H2", "", -1), # UNDONE: isn't this hydrogenation (and vice versa)?
        list("hydrogenation", "", "H2", -1),
        list("dehydration", "", "H2O", -1),
        list("chlorine_reduction", "H", "Cl", -1),
        list("acetylation", "C2H2O", "", -1),
        list("deacetylation", "", "C2H2O", -1),
        list("glucuronidation", "C6H8O6", "", -1),
        list("deglucuronidation", "", "C8H8O6", 0), # UNDONE: or RTDir=1
        list("sulfonation", "SO3", "", -1),
        list("desulfonation", "", "SO3", -1)
    ))
    
    # RTDir: -1: <= precursor; 0: no check; 1: >= precursor
    setnames(ret, c("reaction", "add", "sub", "RTDir"))
    
    getMZ <- function(f) rcdk::get.formula(f)@mass
    
    ret[, deltaMZ := 0]
    ret[nzchar(add), deltaMZ := sapply(add, getMZ)]
    ret[nzchar(sub), deltaMZ := deltaMZ - sapply(sub, getMZ)]
    
    return(ret[])
}

doPredictTPsLogic <- function(fGroups, minMass, adduct)
{
    if (!is.null(adduct))
        adductMZ <- adductMZDelta(adduct)
    
    gInfo <- groupInfo(fGroups)
    neutralMasses <- if (!is.null(adduct)) gInfo$mzs - adductMZ else gInfo$mzs
    suspects <- data.table(name = rownames(gInfo), rt = gInfo$rts, neutralMass = neutralMasses)
    transformations <- getTPLogicTransformations()
    
    prog <- openProgBar(0, nrow(suspects))
    
    predictions <- lapply(seq_len(nrow(suspects)), function(si)
    {
        ret <- data.table(name = paste0(suspects$name[si], "-",
                                        transformations$reaction),
                          neutralMass = neutralMasses[si] + transformations$deltaMZ,
                          RTDir = transformations$RTDir)
        ret <- ret[neutralMass >= minMass]
        
        # UNDONE: more checks (e.g. formulas)
        
        setTxtProgressBar(prog, si)
        
        return(ret)
    })
    
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
setMethod("predictTPsLogic", "featureGroups", function(fGroups, minMass = 40, adduct)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertNumber(minMass, finite = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    adduct <- checkAndToAdduct(adduct)
    
    res <- doPredictTPsLogic(fGroups, minMass, adduct)
    
    return(TPPredictionsLogic(suspects = res$suspects, predictions = res$predictions))
})

#' @export
setMethod("predictTPsLogic", "featureGroupsSet", function(fGroups, minMass = 40)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertNumber(minMass, finite = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    res <- doPredictTPsLogic(fGroups, minMass, NULL)
    
    return(TPPredictionsLogic(suspects = res$suspects, predictions = res$predictions))
})
