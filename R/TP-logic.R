#' @include main.R
#' @include TP.R
NULL

# UNDONE: precursor --> parent?

getTPLogicTransformations <- function()
{
    # UNDONE: cite 10.1021/acs.analchem.5b02905
    
    ret <- rbindlist(list(
        list("hydroxylation", "O", "", TRUE),
        list("demethylation", "", "CH2", TRUE),
        list("deethylation", "", "C2H4", TRUE),
        list("dehydrogenation", "H2", "", TRUE), # UNDONE: isn't this hydrogenation (and vice versa)?
        list("hydrogenation", "", "H2", TRUE),
        list("dehydration", "", "H2O", TRUE),
        list("chlorine_reduction", "H", "Cl", TRUE),
        list("acetylation", "C2H2O", "", TRUE),
        list("deacetylation", "", "C2H2O", TRUE),
        list("glucuronidation", "C6H8O6", "", TRUE),
        list("deglucuronidation", "", "C8H8O6", FALSE),
        list("sulfonation", "SO3", "", TRUE),
        list("desulfonation", "", "SO3", TRUE)
    ))
    setnames(ret, c("reaction", "add", "sub", "lowerRT"))
    
    getMZ <- function(f) rcdk::get.formula(f)@mass
    
    ret[, delta_mz := 0]
    ret[nzchar(add), delta_mz := sapply(add, getMZ)]
    ret[nzchar(sub), delta_mz := delta_mz - sapply(sub, getMZ)]
    
    return(ret[])
}

#' @export
TPPredictionsLogic <- setClass("TPPredictionsLogic", contains = "TPPredictions")

setMethod("initialize", "TPPredictionsLogic",
          function(.Object, ...) callNextMethod(.Object, algorithm = "logic", ...))

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
                          mass = mass + transformations$delta_mz,
                          mz = suspects$mz[si] + transformations$delta_mz,
                          lowerRT = transformations$lowerRT)
        ret <- ret[mass >= minMass]
        
        # UNDONE: more checks (e.g. formulas)

        setTxtProgressBar(prog, si)

        return(ret)
    })
    
    setTxtProgressBar(prog, nrow(suspects))
    close(prog)
    
    
    return(TPPredictionsLogic(suspects = suspects, predictions = predictions))
}
