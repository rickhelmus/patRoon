#' @include main.R
#' @include workflow-step-set.R
#' @include features.R
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
        
        if (nrow(fTab) == 0)
            fTab[, adduct := character()]
        else
        {
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
        }
        
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
        ret$set <- if (nrow(ret) == 0) character() else set
        return(ret)
    }))
    
    # combine (neutralized) features
    combFeatures <- Reduce(modifyList, lapply(neutralizedFeatures, featureTable))
    
    return(featuresSet(features = combFeatures, analysisInfo = combAnaInfo,
                       algorithm = makeSetAlgorithm(featuresList)))
}

#' @param set \setsWF The name of the set.
#' @param sets \setsWF For \code{[} and \code{filter}: a \code{character} with name(s) of the sets to keep (or remove if
#'   \code{negate=TRUE}).
#'
#' @section Sets workflows: \setsWFClass{featuresSet}{features}
#'
#'   \setsWFNewMethodsFeat{featuresUnset}{The adduct annotations for the selected set (\emph{e.g.} as passed to
#'   \code{makeSet}) are used to convert all feature masses to ionic \emph{m/z} values. }
#'
#'   \setsWFChangedMethods{
#'
#'   \item \code{filter} and the subset operator (\code{[}) have specific arguments to choose/filter by (feature
#'   presence in) sets. See the \code{sets} argument description.
#'
#'   }
#'
#' @rdname features-class
#' @export
featuresSet <- setClass("featuresSet", contains = "features")

#' @rdname features-class
#' @export
setMethod("sets", "featuresSet", function(obj) unique(analysisInfo(obj)$set))

#' @rdname features-class
#' @export
setMethod("show", "featuresSet", function(object)
{
    callNextMethod()
    printf("Sets: %s\n", paste0(sets(object), collapse = ", "))
})

#' @rdname features-class
#' @export
setMethod("as.data.table", "featuresSet", function(x)
{
    ret <- callNextMethod(x)
    anaInfo <- analysisInfo(x)
    ret[, set := anaInfo$set[match(analysis, anaInfo$analysis)]]
    setcolorder(ret, "set")
    return(ret[])
})

#' @rdname features-class
#' @export
setMethod("[", c("featuresSet", "ANY", "missing", "missing"), function(x, i, ..., sets = NULL, drop = TRUE)
{
    assertSets(x, sets, TRUE)
    
    if (!is.null(sets))
        i <- mergeAnaSubsetArgWithSets(i, sets, analysisInfo(x))
        
    x <- callNextMethod(x, i, ...)
    
    return(x)
})

#' @rdname features-class
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
            sets <- setdiff(get("sets", pos = 2)(obj), sets)
        obj <- obj[, sets = sets]
    }
    
    if (...length() > 0)
        obj <- callNextMethod(obj, ..., negate = negate)
    
    return(obj)
})

#' @template makeSet
#' @rdname feature-finding
#' @return \code{makeSet} returns a \code{\link{featuresSet}} object.
#' @export
setMethod("makeSet", "features", function(obj, ..., adducts, labels = NULL)
{
    # UNDONE: cache
    
    featuresList <- list(obj, ...)
    ac <- checkmate::makeAssertCollection()
    assertMakeSetArgs(featuresList, "features", adducts, FALSE, labels, ac)
    checkmate::reportAssertions(ac)

    allAnas <- unlist(lapply(featuresList, analyses))
    if (anyDuplicated(allAnas))
        stop("Some objects have non-unique analyses: ", paste0(unique(allAnas[duplicated(allAnas)]), collapse = ","))
    
    adducts <- prepareMakeSetAdducts(featuresList, adducts, labels)
    names(featuresList) <- names(adducts)

    return(doMakeFeaturesSet(featuresList, adducts))
})

#' @note \code{makeSet} Currently does not support making sets from \code{\link{featuresSet}} objects.
#' @rdname feature-finding
#' @export
setMethod("makeSet", "featuresSet", function(obj, ...)
{
    stop("Making a set from set objects is not supported", call. = FALSE)
})

#' @rdname features-class
#' @export
featuresUnset <- setClass("featuresUnset", contains = "features")

#' @rdname features-class
#' @export
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
