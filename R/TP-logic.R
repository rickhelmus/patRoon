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

#' @rdname transformationProducts-class
#' @export
transformationProductsLogic <- setClass("transformationProductsLogic", contains = "transformationProducts")

setMethod("initialize", "transformationProductsLogic",
          function(.Object, ...) callNextMethod(.Object, algorithm = "logic", ...))

setMethod("linkParentsToFGroups", "transformationProductsLogic", function(TPs, fGroups)
{
    fg <- intersect(names(TPs), names(fGroups))
    return(data.table(name = fg, group = fg))
})

#' Obtain transformation products (TPs) with metabolic logic
#'
#' Automatically calculate potential transformation products with \emph{metabolic logic}.
#'
#' @templateVar algo metabolic logic
#' @templateVar do obtain transformation products
#' @templateVar generic generateTPs
#' @templateVar algoParam logic
#' @template algo_generator
#'
#' @details With this algorithm TPs are predicted from common (environmental) chemical reactions, such as hydroxylation,
#'   demethylation etc. The generated TPs result from calculating the mass differences between a parent feature after it
#'   underwent the reaction. While this only results in little information on chemical properties of the TP, an
#'   advantage of this method is that it does not rely on structural information of the parent, which may be unknown in
#'   a full non-target analysis.
#'
#' @param fGroups A \code{\link{featureGroups}} object for which TPs should be calculated.
#' @param minMass A \code{numeric} that specifies the minimum mass of calculated TPs. If below this mass it will be
#'   removed.
#' @param transformations A \code{data.frame} with transformation reactions to be used for calculating the TPs (see
#'   details below). If \code{NULL}, a default table from Schollee \emph{et al.} is used (see references).
#'
#' @template adduct-arg
#'
#' @inherit generateTPs return
#'
#' @section Custom transformations: The \code{transformations} argument to \code{generateTPsLogic} is used to specify
#'   custom rules to calculate transformation products. This should be a \code{data.frame} with the following columns:
#'   \itemize{
#'
#'   \item \code{transformation} The name of the chemical transformation
#'
#'   \item \code{add} The elements that are added by this reaction (\emph{e.g.} \code{"O"}).
#'
#'   \item \code{sub} The elements that are removed by this reaction (\emph{e.g.} \code{"H2O"}).
#'
#'   \item \code{retDir} The expected retention time direction relative to the parent (assuming a reversed phase like LC
#'   separation). Valid values are: \samp{-1} (elutes before the parent), \samp{1} (elutes after the parent) or \samp{0}
#'   (no significant change or unknown).
#'
#'   }
#'
#' @section Source: The algorithm of \code{generateTPsLogic} is directly based on the work done by Schollee \emph{et
#'   al.} (see references).
#'
#' @references \insertRef{Scholle2015}{patRoon}
#'
#' @name generateTPsLogic
#' @export
setMethod("generateTPsLogic", "featureGroups", function(fGroups, minMass = 40, adduct = NULL, transformations = NULL)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertNumber(minMass, lower = 0, finite = TRUE, add = ac)
    assertLogicTransformations(transformations, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    adduct <- checkAndToAdduct(adduct, fGroups)
    fgAdd <- getFGroupAdducts(names(fGroups), annotations(fGroups), adduct, "generic")
    neutralMasses <- if (length(fGroups) == 0) numeric() else calculateMasses(groupInfo(fGroups)[, "mzs"],
                                                                              fgAdd$grpAdducts, type = "neutral")
    res <- doGenerateTPsLogic(fGroups, minMass, neutralMasses, transformations)
    
    return(transformationProductsLogic(parents = res$parents, products = res$products))
})

#' @rdname generateTPsLogic
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
 
