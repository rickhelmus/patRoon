# UNDONE: design common set base class for this, compounds etc?

#' @include main.R
NULL

neutralizeFeatures <- function(feat, adduct)
{
    adductMZ <- adductMZDelta(adduct)
    feat@features <- lapply(feat@features, function(fTab)
    {
        fTab <- copy(fTab)
        fTab[, mz := mz - adductMZ]
        fTab[, mzmin := mz - adductMZ]
        fTab[, mzmax := mz - adductMZ]
        return(fTab)
    })
    return(feat)
}

featuresSet <- setClass("featuresSet",
                        slots = c(adducts = "list", setObjects = "list",
                                  ionizedFeatures = "list"),
                        contains = "features")

setMethod("initialize", "featuresSet",
          function(.Object, ...) callNextMethod(.Object, algorithm = "set", ...))


# UNDONE: move (partially)
#' @describeIn featuresSet Shows summary information for this object.
#' @export
setMethod("show", "featuresSet", function(object)
{
    callNextMethod(object)
    printf("sets: %s\n", paste0(sets(object), collapse = ", "))
    printf("adducts: %s\n", paste0(sapply(object@adducts, as.character), collapse = ", "))
})

# UNDONE: export/move
setMethod("sets", "featuresSet", function(obj) names(obj@adducts))

#' @describeIn featuresSet Get table with feature information
#'
#' @return \code{featureTable}: A \code{list} containing a
#'   \code{\link{data.table}} for each analysis with feature data
#'
#' @export
setMethod("featureTable", "featuresSet", function(obj, neutralized = TRUE, set = NULL)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFlag(neutralized, add = ac)
    checkmate::assertChoice(set, sets(obj), null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    if (!is.null(set))
        obj <- obj[, set = set]
    
    if (neutralized)
        return(callNextMethod(obj))
    
    return(obj@ionizedFeatures)
})


#' @describeIn featuresSet Returns all feature data in a table.
#' @export
setMethod("as.data.table", "featuresSet", function(x, neutralized = TRUE, set = NULL)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFlag(neutralized, add = ac)
    checkmate::assertChoice(set, sets(obj), null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (!is.null(set))
        x <- x[, set = set]
    
    if (neutralized)
    {
        ret <- callNextMethod(x)
        ret[, set := rep.int(sets(x), times = sapply(x@setObjects, length))]
        setcolorder(ret, "set")
        return(ret)
    }
    
    # return original ionized features
    return(rbindlist(lapply(lapply(x@setObjects, featureTable), rbindlist, idcol = "analysis", fill = TRUE),
                     idcol = "set", fill = TRUE))
})

#' @describeIn featuresSet Subset on analyses.
#' @param \dots Ignored.
#' @export
setMethod("[", c("featuresSet", "ANY", "missing", "missing"), function(x, i, ..., set = NULL, drop = TRUE)
{
    checkmate::assertChoice(set, sets(obj), null.ok = TRUE)
    if (!is.null(set))
    {
        setAna <- x@analysisInfo$analysis[x@analysisInfo$set == set]
        i <- if (missing(i)) setAna else union(i, setAna)
    }
    
    x <- callNextMethod(x, i, ...)
    if (!missing(i))
    {
        subSets <- unique(x@analysisInfo$set)
        x@adducts <- x@adducts[subSets]
        x@setObjects <- x@setObjects[subSets]
        x@ionizedFeatures <- x@ionizedFeatures[x@analysisInfo$analysis]
    }
    return(x)
})

#' @describeIn features Extract a feature table for an analysis.
#' @export
setMethod("[[", c("featuresSet", "ANY", "missing"), function(x, i, neutralized)
{
    ac <- checkmate::makeAssertCollection()
    assertExtractArg(i, add = ac)
    checkmate::assertFlag(neutralized, add = ac)
    checkmate::reportAssertions(ac)
    
    return(if (neutralized) x@features[[i]] else x@ionizedFeatures[[i]])
})

# UNDONE: re-use docs from this and other methods somehow?
# UNDONE: support ionized masses for mass filters? or just clarify it doesn't.
#' @describeIn features Performs common rule based filtering of features. Note
#'   that this (and much more) functionality is also provided by the
#'   \code{filter} method defined for \code{\link{featureGroups}}. However,
#'   filtering a \code{features} object may be useful to avoid grouping large
#'   amounts of features.
#' @templateVar feat TRUE
#' @template feat-filter-args
#' @export
setMethod("filter", "featuresSet", function(obj, ..., set = NULL)
{
    checkmate::assertChoice(set, sets(obj), null.ok = TRUE)
    
    if (!is.null(sets))
        obj <- obj[, sets = sets]
    
    obj <- callNextMethod(obj, ...)
    
    # synchronize other features objects by remaining IDs
    cat("Synchronizing feature set objects...")
    remainingIDsPerAna <- sapply(obj@features, "[[", "ID", simplify = FALSE)
    obj@ionizedFeatures <- lapply(analyses(obj), function(ana) obj@ionizedFeatures[[ana]][ID %in% remainingIDsPerAna[[ana]]])
    obj@setObjects <- lapply(obj@setObjects, function(so) lapply(names(so),
                                                                 function(ana) featureTable(so)[[ana]][ID %in% remainingIDsPerAna]))
    cat("Done!\n")
    
    return(obj)
})

setMethod("c", "features", function(x, ..., adducts)
{
    # UNDONE: check anaInfos to be unique
    # UNDONE: cache
    
    featuresList <- list(x, ...)
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertList(featuresList, types = "features", any.missing = FALSE,
                          unique = TRUE, .var.name = "...", add = ac)
    checkmate::assert(checkmate::checkCharacter(adducts, any.missing = FALSE, min.len = 1,
                                                max.len = length(featuresList)),
                      checkmate::checkList(adducts, types = c("adduct", "character"), any.missing = FALSE,
                                           min.len = 1, max.len = length(featuresList)),
                      .var.name = "adducts")
    checkmate::reportAssertions(ac)
    
    n <- getArgNames(..., def = sapply(featuresList, algorithm))
    names(featuresList) <- make.unique(n)
    
    adductNamed <- checkmate::testNames(names(adducts))
    if (adductNamed)
        checkmate::assertNames(names(adducts), type = "unique", must.include = names(featuresList))
    adducts <- lapply(adducts, checkAndToAdduct, .var.name = "adducts")
    adducts <- rep(adducts, length.out = length(featuresList))
    
    if (!adductNamed)
        names(adducts) <- names(featuresList)
    else
        adducts <- adducts[names(featuresList)] # synchronize order
    
    # neutralize features
    neutralizedFeatures <- mapply(featuresList, adducts, FUN = neutralizeFeatures,
                                  SIMPLIFY = FALSE, USE.NAMES = TRUE)
    
    # combine anaInfo and tag
    combAnaInfo <- do.call(rbind, lapply(names(featuresList), function(set)
    {
        ret <- featuresList[[set]]@analysisInfo
        ret$set <- set
        return(ret)
    }))
    
    # combine (neutralized) features
    combFeatures <- Reduce(modifyList, lapply(neutralizedFeatures, featureTable))
    combFeaturesIon <- Reduce(modifyList, lapply(featuresList, featureTable))

    return(featuresSet(adducts = adducts, setObjects = featuresList, ionizedFeatures = combFeaturesIon,
                       features = combFeatures, analysisInfo = combAnaInfo))
})

# UNDONE: implement someday?
fSetNotYetImplemented <- function() stop("This function is not yet implemented for featuresSet.", call. = FALSE)

setMethod("screenSuspects", "featuresSet", function(obj, suspects, rtWindow, mzWindow, adduct,
                                                    skipInvalid) fSetNotYetImplemented())
setMethod("getXCMSSet", "featuresSet", function(obj, verbose) fSetNotYetImplemented())
setMethod("getXCMSnExp", "featuresSet", function(obj, verbose) fSetNotYetImplemented())
