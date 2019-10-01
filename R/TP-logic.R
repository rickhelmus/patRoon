#' @include main.R
#' @include TP.R
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

#' @export
TPPredictionsLogic <- setClass("TPPredictionsLogic", contains = "TPPredictions")

setMethod("initialize", "TPPredictionsLogic",
          function(.Object, ...) callNextMethod(.Object, algorithm = "logic", ...))

setMethod("linkPrecursorsToFGroups", "TPPredictionsLogic", function(pred, fGroups, adduct, mzWindow)
{
    # here we can just link directly (ie suspects are the groups)
    return(data.table(name = suspects(pred)$name,
                      group = suspects(pred)$name))
})

#' @export
predictTPsLogic <- function(fGroups, adduct, minMass = 40)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    checkmate::assertNumber(minMass, finite = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    adduct <- checkAndToAdduct(adduct)
    adductMZ <- adductMZDelta(adduct)
    
    gInfo <- groupInfo(fGroups)
    suspects <- data.table(name = rownames(gInfo), rt = gInfo$rts, mz = gInfo$mzs)
    transformations <- getTPLogicTransformations()

    prog <- openProgBar(0, nrow(suspects))
    
    predictions <- lapply(seq_len(nrow(suspects)), function(si)
    {
        mass <- suspects$mz[si] - adductMZ
        
        ret <- data.table(name = paste0(suspects$name[si], "-",
                                              transformations$reaction),
                          mass = mass + transformations$deltaMZ,
                          mz = suspects$mz[si] + transformations$deltaMZ,
                          RTDir = transformations$RTDir)
        ret <- ret[mass >= minMass]
        
        # UNDONE: more checks (e.g. formulas)

        setTxtProgressBar(prog, si)

        return(ret)
    })
    
    setTxtProgressBar(prog, nrow(suspects))
    close(prog)
    
    
    return(TPPredictionsLogic(suspects = suspects, predictions = predictions))
}
