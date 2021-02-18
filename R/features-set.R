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
    
    return(featuresSet(features = combFeatures, analysisInfo = combAnaInfo,
                       algorithm = makeSetAlgorithm(featuresList)))
}

#' @export
featuresSet <- setClass("featuresSet", contains = "features")

#' @export
setMethod("sets", "featuresSet", function(obj) unique(analysisInfo(obj)$set))

#' @describeIn featuresSet Shows summary information for this object.
#' @export
setMethod("show", "featuresSet", function(object)
{
    callNextMethod()
    printf("Sets: %s\n", paste0(sets(object), collapse = ", "))
})

#' @describeIn featuresSet Returns all feature data in a table.
#' @export
setMethod("as.data.table", "featuresSet", function(x)
{
    ret <- callNextMethod(x)
    anaInfo <- analysisInfo(x)
    ret[, set := anaInfo$set[match(analysis, anaInfo$analysis)]]
    setcolorder(ret, "set")
    return(ret[])
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
    
    return(x)
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
    
    if (...length() > 0)
        obj <- callNextMethod(obj, ..., negate = negate)
    
    return(obj)
})

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
    
    ionizedFTable <- lapply(featureTable(obj), function(ft)
    {
        ft <- copy(ft)
        
        if (nrow(ft) > 0)
        {
            adducts <- sapply(unique(ft$adduct), as.adduct)
            addMZs <- sapply(adducts, adductMZDelta)
            addMZs <- addMZs[ft$adduct]
            
            set(ft, j = c("mz", "mzmin", "mzmax"),
                value = list(ft$mz + addMZs, ft$mzmin + addMZs, ft$mzmax + addMZs))
        }
        
        ft[, adduct := NULL] # UNDONE: keep?
        
        return(ft[])
    })
    
    return(featuresUnset(features = ionizedFTable, analysisInfo = unSetAnaInfo(analysisInfo(obj)),
                         algorithm = paste0(algorithm(obj), "_unset")))
})
