#' @include main.R
#' @include workflow-step-set.R
NULL

neutralizeFeatures <- function(feat, adduct)
{
    if (!is.null(adduct))
    {
        adductChar <- as.character(adduct)
        adductMZ <- adductMZDelta(adduct)
    }
    else
    {
        allAdducts <- unique(unlist(lapply(feat@features, "[[", "adduct")))
        adductMZ <- sapply(allAdducts, function(a) adductMZDelta(as.adduct(a)))
    }
    
    feat@features <- lapply(feat@features, function(fTab)
    {
        fTab <- copy(fTab)
        
        if (!is.null(adduct))
        {
            mzd <- adductMZ
            fTab[, adduct := adductChar]
        }
        else
        {
            mzd <- adductMZ[fTab$adduct]
            fTab[, adduct := fTab$adduct]
        }
        
        fTab[, mz := mz - mzd]
        fTab[, mzmin := mzmin - mzd]
        fTab[, mzmax := mzmax - mzd]
        
        return(fTab)
    })
    return(feat)
}

doMakeFeaturesSet <- function(featuresList, adducts)
{
    if (!is.null(adducts))
        neutralizedFeatures <- mapply(featuresList, adducts, FUN = neutralizeFeatures,
                                      SIMPLIFY = FALSE, USE.NAMES = TRUE)
    else
        neutralizedFeatures <- sapply(featuresList, neutralizeFeatures, adduct = NULL,
                                      simplify = FALSE)
    
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
                       features = combFeatures, analysisInfo = combAnaInfo, algorithm = makeSetAlgorithm(featuresList)))
}

#' @export
featuresSet <- setClass("featuresSet",
                        slots = c(ionizedFeatures = "list"),
                        contains = c("features", "workflowStepSet"))


#' @describeIn featuresSet Shows summary information for this object.
#' @export
setMethod("show", "featuresSet", function(object)
{
    callAllNextMethods(object, show, firstClass = "features", startFrom = "featuresSet")
})

#' @describeIn featuresSet Get table with feature information
#'
#' @return \code{featureTable}: A \code{list} containing a
#'   \code{\link{data.table}} for each analysis with feature data
#'
#' @export
setMethod("featureTable", "featuresSet", function(obj, neutralized = TRUE)
{
    checkmate::assertFlag(neutralized)
    
    if (neutralized)
        return(callNextMethod(obj))
    
    return(obj@ionizedFeatures)
})


#' @describeIn featuresSet Returns all feature data in a table.
#' @export
setMethod("as.data.table", "featuresSet", function(x, neutralized = TRUE)
{
    checkmate::assertFlag(neutralized)
    
    if (neutralized)
    {
        ret <- callNextMethod(x)
        ret[, set := rep.int(sets(x), times = sapply(x@setObjects, length))]
        setcolorder(ret, "set")
        return(ret[])
    }
    
    # return original ionized features
    return(rbindlist(lapply(lapply(x@setObjects, featureTable), rbindlist, idcol = "analysis", fill = TRUE),
                     idcol = "set", fill = TRUE))
})

#' @describeIn featuresSet Subset on analyses.
#' @param \dots Ignored.
#' @export
setMethod("[", c("featuresSet", "ANY", "missing", "missing"), function(x, i, ..., sets = NULL, drop = TRUE)
{
    assertSets(x, sets, TRUE)
    
    if (!is.null(sets) && length(sets) > 0)
        i <- mergeAnaSubsetArgWithSets(i, sets, analysisInfo(x))
        
    x <- callNextMethod(x, i, ...)
    if (!missing(i))
    {
        subSets <- unique(x@analysisInfo$set)
        x@adducts <- x@adducts[subSets]
        # NOTE: assume that subsetting with non-existing analyses will not result in errors
        x@setObjects <- lapply(x@setObjects[subSets], "[", i = analyses(x))
        x@ionizedFeatures <- x@ionizedFeatures[x@analysisInfo$analysis]
    }
    return(x)
})

#' @describeIn features Extract a feature table for an analysis.
#' @export
setMethod("[[", c("featuresSet", "ANY", "missing"), function(x, i, neutralized = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    assertExtractArg(i, add = ac)
    checkmate::assertFlag(neutralized, add = ac)
    checkmate::reportAssertions(ac)
    
    return(if (neutralized) x@features[[i]] else x@ionizedFeatures[[i]])
})

#' @describeIn features Performs common rule based filtering of features. Note
#'   that this (and much more) functionality is also provided by the
#'   \code{filter} method defined for \code{\link{featureGroups}}. However,
#'   filtering a \code{features} object may be useful to avoid grouping large
#'   amounts of features.
#' @templateVar feat TRUE
#' @template feat-filter-args
#' @export
setMethod("filter", "featuresSet", function(obj, ..., negate = FALSE, sets = NULL)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFlag(negate, add = ac)
    assertSets(obj, sets, TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (!is.null(sets) && length(sets) > 0)
    {
        if (negate)
            sets <- setdiff(obj@sets, sets)
        obj <- obj[, sets = sets]
    }
    
    if (length(list(...)) > 0)
    {
        obj <- callNextMethod(obj, ..., negate = negate)
        
        # synchronize other features objects by remaining IDs
        cat("Synchronizing feature set objects...")
        remainingIDsPerAna <- sapply(obj@features, "[[", "ID", simplify = FALSE)
        obj@ionizedFeatures <- lapply(analyses(obj), function(ana) obj@ionizedFeatures[[ana]][ID %in% remainingIDsPerAna[[ana]]])
        obj@setObjects <- lapply(obj@setObjects, function(so) lapply(names(so),
                                                                     function(ana) featureTable(so)[[ana]][ID %in% remainingIDsPerAna]))
        cat("Done!\n")
    }
    
    return(obj)
})

# UNDONE: mention that object will be unset
#' @export
setMethod("getXCMSSet", "featuresSet", function(obj, ..., set) getXCMSSet(unset(obj, set), ...))

# UNDONE: mention that object will be unset
#' @export
setMethod("getXCMSnExp", "featuresSet", function(obj, ..., set) getXCMSnExp(unset(obj, set), ...))

#' @export
setMethod("makeSet", "features", function(obj, ..., adducts, labels)
{
    # UNDONE: check anaInfos to be unique
    # UNDONE: cache
    
    featuresList <- list(obj, ...)
    ac <- checkmate::makeAssertCollection()
    assertMakeSetArgs(featuresList, "features", adducts, FALSE, labels, ac)
    checkmate::reportAssertions(ac)

    adducts <- prepareMakeSetAdducts(featuresList, adducts, labels)
    names(featuresList) <- names(adducts)

    return(doMakeFeaturesSet(featuresList, adducts))
})

featuresUnset <- setClass("featuresUnset", contains = "features")
setMethod("unset", "featuresSet", function(obj, set)
{
    assertSets(obj, set, FALSE)
    obj <- obj[, sets = set]
    return(featuresUnset(features = obj@ionizedFeatures, analysisInfo = analysisInfo(obj),
                         algorithm = paste0(algorithm(obj), "_unset")))
})
