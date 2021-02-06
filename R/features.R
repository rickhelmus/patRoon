#' @include main.R
#' @include workflow-step.R
NULL

printFeatStats <- function(fList)
{
    fCounts <- sapply(fList, nrow)
    fTotCount <- sum(fCounts)
    printf("Feature statistics:\n")
    printf("%s: %d (%.1f%%)\n", names(fList), fCounts, if (fTotCount == 0) 0 else fCounts * 100 / fTotCount)
    printf("Total: %d\n", fTotCount)
}

#' Base features class
#'
#' Holds information for all features present within a set of analysis.
#'
#' This class provides a way to store intensity, retention times, \emph{m/z} and
#' other data for all features in a set of analyses. The class is \code{virtual}
#' and derived objects are created by 'feature finders' such as
#' \code{findFeaturesOpenMS}, \code{findFeaturesXCMS} and
#' \code{findFeaturesBruker}.
#'
#' @param obj,x,object \code{features} object to be accessed
#'
#' @seealso \code{\link{feature-finding}}
#'
#' @slot features List of features per analysis file. Use the
#'   \code{featureTable} method for access.
#' @slot analysisInfo Analysis group information. Use the \code{analysisInfo} method
#'   for access.
#'
#' @templateVar seli analyses
#' @templateVar selOrderi analyses()
#' @templateVar dollarOpName analysis
#' @template sub_op-args
#'
#' @param \dots Ignored.
#'
#' @templateVar class features
#' @template class-hierarchy
#'
#' @export
features <- setClass("features",
                     slots = c(features = "list", analysisInfo = "data.frame"),
                     contains = c("VIRTUAL", "workflowStep"))

setMethod("initialize", "features", function(.Object, ...)
{
    .Object <- callNextMethod(.Object, ...)
    .Object@features <- makeEmptyListNamed(.Object@features)
    return(.Object)
})

#' @describeIn features Obtain total number of features.
#' @export
setMethod("length", "features", function(x) if (length(x@features) > 0) sum(sapply(x@features, nrow)) else 0)

#' @describeIn features Shows summary information for this object.
#' @export
setMethod("show", "features", function(object)
{
    callNextMethod(object)
    ftcounts <- if (length(object@features) > 0) sapply(object@features, nrow) else 0
    printf("Total feature count: %d\n", sum(ftcounts))
    printf("Average feature count/analysis: %.0f\n", if (length(object) > 0) sum(ftcounts) / nrow(analysisInfo(object)) else 0)
    printf("Least features: %s\n", names(object)[which.min(ftcounts)])
    printf("Most features: %s\n", names(object)[which.max(ftcounts)])
    showAnaInfo(analysisInfo(object))
})

#' @describeIn features Get table with feature information
#'
#' @return \code{featureTable}: A \code{list} containing a
#'   \code{\link{data.table}} for each analysis with feature data
#'
#' @export
setMethod("featureTable", "features", function(obj) obj@features)

setReplaceMethod("featureTable", "features", function(obj, value)
{
    # UNDONE: verify value
    obj@features <- value
    return(obj)
})

#' @describeIn features Get analysis information
#' @return \code{analysisInfo}: A \code{data.frame} containing a column with
#'   analysis name (\code{analysis}), its path (\code{path}), and other columns
#'   such as replicate group name (\code{group}) and blank reference
#'   (\code{blank}).
#' @export
setMethod("analysisInfo", "features", function(obj) obj@analysisInfo)

#' @templateVar class features
#' @templateVar what analyses
#' @template strmethod
#' @export
setMethod("analyses", "features", function(obj) analysisInfo(obj)$analysis)

#' @templateVar class features
#' @templateVar what replicate groups
#' @template strmethod
#' @export
setMethod("replicateGroups", "features", function(obj) unique(analysisInfo(obj)$group))

#' @describeIn features Returns all feature data in a table.
#' @export
setMethod("as.data.table", "features", function(x) rbindlist(featureTable(x), idcol = "analysis", fill = TRUE))

#' @describeIn features Performs common rule based filtering of features. Note
#'   that this (and much more) functionality is also provided by the
#'   \code{filter} method defined for \code{\link{featureGroups}}. However,
#'   filtering a \code{features} object may be useful to avoid grouping large
#'   amounts of features.
#' @templateVar feat TRUE
#' @template feat-filter-args
#' @export
setMethod("filter", "features", function(obj, absMinIntensity = NULL, relMinIntensity = NULL,
                                         retentionRange = NULL, mzRange = NULL, mzDefectRange = NULL,
                                         chromWidthRange = NULL, negate = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertNumber, . ~ absMinIntensity + relMinIntensity, lower = 0, finite = TRUE,
           null.ok = TRUE, fixed = list(add = ac))
    aapply(assertRange, . ~ retentionRange + mzRange + mzDefectRange + chromWidthRange, null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertFlag(negate, add = ac)
    checkmate::reportAssertions(ac)

    if (length(obj) == 0)
        return(obj)

    oldn <- length(obj)

    hash <- makeHash(obj, absMinIntensity, relMinIntensity, retentionRange, mzRange, mzDefectRange, chromWidthRange, negate)
    cache <- loadCacheData("filterFeatures", hash)
    if (!is.null(cache))
        obj <- cache
    else
    {
        anaInfo <- analysisInfo(obj)

        absIntPred <- if (!negate) function(x) x >= absMinIntensity else function(x) x < absMinIntensity
        relIntPred <- if (!negate) function(x, m) (x / m) >= relMinIntensity else function(x, m) (x / m) < relMinIntensity
        rangePred <- function(x, range) numGTE(x, range[1]) & numLTE(x, range[2])

        if (negate)
            rangePred <- Negate(rangePred)

        fList <- obj@features
        for (ana in analyses(obj))
        {
            if (!is.null(absMinIntensity))
                fList[[ana]] <- fList[[ana]][absIntPred(intensity)]

            if (!is.null(relMinIntensity))
            {
                maxInt <- max(fList[[ana]]$intensity)
                fList[[ana]] <- fList[[ana]][relIntPred(intensity, maxInt)]
            }

            if (!is.null(retentionRange))
                fList[[ana]] <- fList[[ana]][rangePred(ret, retentionRange)]

            if (!is.null(mzRange))
                fList[[ana]] <- fList[[ana]][rangePred(mz, mzRange)]

            if (!is.null(mzDefectRange))
                fList[[ana]] <- fList[[ana]][rangePred(mz - floor(mz), mzDefectRange)]

            if (!is.null(chromWidthRange))
                fList[[ana]] <- fList[[ana]][rangePred(retmax - retmin, chromWidthRange)]
        }

        featureTable(obj) <- fList
        
        saveCacheData("filterFeatures", obj, hash)
    }

    newn <- length(obj)
    printf("Done! Filtered %d (%.2f%%) features. Remaining: %d\n", oldn - newn, if (oldn == 0) 0 else (1-(newn/oldn))*100, newn)

    return(obj)
})

#' @describeIn features Subset on analyses.
#' @export
setMethod("[", c("features", "ANY", "missing", "missing"), function(x, i, ...)
{
    if (!missing(i))
    {
        i <- assertSubsetArgAndToChr(i, analyses(x))
        x@features <- x@features[i]
        x@analysisInfo <- x@analysisInfo[x@analysisInfo$analysis %in% i, ]
    }

    return(x)
})

setReplaceMethod("[", c("features", "ANY", "missing"), function(x, i, j, value)
{
    # UNDONE: verify value
    i <- assertSubsetArgAndToChr(i, analyses(x))
    x@features[i] <- value
    return(x)
})

#' @describeIn features Extract a feature table for an analysis.
#' @export
setMethod("[[", c("features", "ANY", "missing"), function(x, i)
{
    assertExtractArg(i)
    return(x@features[[i]])
})

setReplaceMethod("[[", c("features", "ANY", "missing"), function(x, i, j, value)
{
    # UNDONE: verify value
    assertExtractArg(i)
    x@features[[i]] <- value
    return(x)
})

#' @describeIn features Extract a feature table for an analysis.
#' @export
setMethod("$", "features", function(x, name)
{
    eval(substitute(x@features$NAME_ARG, list(NAME_ARG = name)))
})

setReplaceMethod("$", "features", function(x, name, value)
{
    eval(substitute(x@features$NAME_ARG <- value, list(NAME_ARG = name)))
    return(x)
})

setMethod("calculatePeakQualities", "features", function(obj, weights, flatnessFactor)
{
    featQualities <- featureQualities()
    featQualityNames <- featureQualityNames()
    featScoreNames <- featureScoreNames()
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertNumeric(weights, finite = TRUE, any.missing = FALSE, min.len = 1, names = "unique",
                             null.ok = TRUE, add = ac)
    if (!is.null(weights))
        checkmate::assertNames(names(weights), subset.of = featScoreNames, add = ac)
    checkmate::assertNumber(flatnessFactor, add = ac)
    checkmate::reportAssertions(ac)
    
    hash <- makeHash(obj, weights, flatnessFactor)
    cd <- loadCacheData("calculatePeakQualities", hash)
    if (!is.null(cd))
        return(cd)
    
    EICs <- getEICsForFeatures(obj)
    
    checkPackage("MetaClean")
    
    # HACK HACK HACK: MetaClean::calculateGaussianSimilarity needs to have
    # xcms::SSgauss attached
    # based on https://stackoverflow.com/a/36611896
    withr::local_environment(list(SSgauss = xcms::SSgauss))
    
    calcFeatQualities <- function(ret, retmin, retmax, intensity, EIC)
    {
        args <- list(c(rt = ret, rtmin = retmin, rtmax = retmax, maxo = intensity), as.matrix(EIC))
        return(sapply(featQualityNames, function(q)
        {
            a <- args
            if (q %in% c("Jaggedness", "Modality"))
                a <- c(a, flatnessFactor)
            return(do.call(featQualities[[q]]$func, a))
        }, simplify = FALSE))
    }
    
    printf("Calculating feature peak qualities and scores...\n")
    prog <- openProgBar(0, length(EICs))
    
    fTable <- Map(featureTable(obj)[names(EICs)], EICs, seq_along(EICs), f = function(ft, eic, i)
    {
        ft <- copy(ft)
        eic <- as.matrix(eic) # MetaClean expects matrices
        ft[, (featQualityNames) := rbindlist(Map(calcFeatQualities, ret, retmin, retmax, intensity, eic))]
        ft[, (featScoreNames) := Map(scoreFeatQuality, featQualities, .SD), .SDcols = featQualityNames]
        setTxtProgressBar(prog, i)
        return(ft)
    })
    
    if (!is.null(weights))
    {
        weights[setdiff(featScoreNames, names(weights))] <- 1
        weights <- weights[featScoreNames]
    }
    
    fTable <- lapply(fTable, function(ft)
    {
        wft <- ft[, featScoreNames, with = FALSE]
        if (!is.null(weights))
            wft[, names(wft) := Map("*", .SD, weights)]
        set(ft, j = "totalScore", value = rowSums(wft, na.rm = TRUE))
    })

    setTxtProgressBar(prog, length(EICs))
    
    featureTable(obj) <- fTable
    
    saveCacheData("calculatePeakQualities", obj, hash)
    
    return(obj)
})


#' @templateVar func findFeatures
#' @templateVar what find features
#' @templateVar ex1 findFeaturesOpenMS
#' @templateVar ex2 findFeaturesBruker
#' @templateVar algos bruker,openms,xcms,xcms3,envipick,kpic2
#' @template generic-algo
#'
#' @rdname feature-finding
#' @aliases findFeatures
#' @export
findFeatures <- function(analysisInfo, algorithm, ..., verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, add = ac)
    checkmate::assertChoice(algorithm, c("bruker", "openms", "xcms", "xcms3", "envipick", "kpic2"), add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)

    f <- switch(algorithm,
                bruker = findFeaturesBruker,
                openms = findFeaturesOpenMS,
                xcms = findFeaturesXCMS,
                xcms3 = findFeaturesXCMS3,
                envipick = findFeaturesEnviPick,
                kpic2 = findfeaturesKPIC2)

    f(analysisInfo, ..., verbose = verbose)
}

#' @details \code{importFeatures} is a generic function to import feature groups
#'   produced by other software. The actual functionality is provided by
#'   specific functions such as \code{importFeaturesXCMS} and
#'   \code{importFeaturesEnviMass}.
#' @param type What type of data should be imported: \code{"xcms"},
#'   \code{"xcms3"} or \code{"envimass"}.
#'
#' @rdname feature-finding
#' @export
importFeatures <- function(analysisInfo, type, ...)
{
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo)

    f <- switch(type,
                xcms = importFeaturesXCMS,
                xcms3 = importFeaturesXCMS3,
                envimass = importFeaturesEnviMass,
                stop("Invalid algorithm! Should be: xcms, xcms3 or envimass"))

    f(analysisInfo = analysisInfo, ...)
}
