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
#' @seealso \code{\link{findFeatures}}
#'
#' @slot features List of features per analysis file. Use the
#'   \code{featureTable} method for access.
#' @slot analysisInfo Analysis group information. Use the \code{analysisInfo} method
#'   for access.
#'
#' @templateVar seli analyses
#' @templateVar selOrderi analyses()
#' @templateVar del TRUE
#' @templateVar deli analyses
#' @templateVar delj features
#' @templateVar deljtype numeric index (row) of the feature
#' @templateVar delfwhat analysis
#' @templateVar delfa the feature table (a \code{data.table}), the analysis name
#' @templateVar delfr the feature indices (rows) to be removed (specified as an \code{integer} or \code{logical} vector)
#' @templateVar dollarOpName analysis
#' @template sub_sel_del-args
#'
#' @param \dots For \code{delete}: passed to the function specified as \code{j}.
#'   
#'   \setsPassedArgs1{features}
#'
#' @templateVar class features
#' @template class-hierarchy
#'
#' @export
features <- setClass("features",
                     slots = c(features = "list", analysisInfo = "data.table"),
                     contains = c("VIRTUAL", "workflowStep"))

setMethod("initialize", "features", function(.Object, ...)
{
    .Object <- callNextMethod(.Object, ...)
    .Object@features <- makeEmptyListNamed(.Object@features)
    return(.Object)
})

#' @describeIn features Obtain total number of features.
#' @export
setMethod("length", "features", function(x)
{
    if (length(x@features) == 0 || !any(lengths(x@features) != 0))
        return(0)
    return(sum(sapply(x@features[lengths(x@features) > 0], nrow)))
})

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
#' @param df If \code{TRUE} then the returned value is a \code{data.frame}, otherwise a \code{data.table}.
#' @return \code{analysisInfo}: A \code{data.table} containing a column with
#'   analysis name (\code{analysis}), its path (\code{path}), and other columns
#'   such as replicate group name (\code{group}) and blank reference
#'   (\code{blank}).
#' @export
setMethod("analysisInfo", "features", function(obj, df = FALSE)
{
    checkmate::assertFlag(df)
    return(if (df) as.data.frame(obj@analysisInfo) else obj@analysisInfo)
})

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
                                         chromWidthRange = NULL, qualityRange = NULL, negate = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertNumber, . ~ absMinIntensity + relMinIntensity, lower = 0, finite = TRUE,
           null.ok = TRUE, fixed = list(add = ac))
    aapply(assertRange, . ~ retentionRange + mzRange + mzDefectRange + chromWidthRange, null.ok = TRUE,
           fixed = list(add = ac))
    assertScoreRange(qualityRange, c(featureQualityNames(group = FALSE),
                                     featureQualityNames(group = FALSE, scores = TRUE)), add = ac)
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
        absIntPred <- if (!negate) function(x) x >= absMinIntensity else function(x) x < absMinIntensity
        relIntPred <- if (!negate) function(x, m) (x / m) >= relMinIntensity else function(x, m) (x / m) < relMinIntensity
        rangePred <- function(x, range) numGTE(x, range[1]) & numLTE(x, range[2])
        if (negate)
            rangePred <- Negate(rangePred)
        scorePred <- function(scTab, qr)
        {
            qualsOK <- scTab[, Map(.SD, qr, f = `%inrange%`), by = seq_len(nrow(scTab))]
            qualsOK <- qualsOK[, -1] # remove dummy by column
            if (negate)
                qualsOK[, keep := any(!unlist(.SD)), by = seq_len(nrow(qualsOK))]
            qualsOK[, keep := all(unlist(.SD)), by = seq_len(nrow(qualsOK))]
            return(qualsOK$keep)
        }

        obj <- delete(obj, j = function(ft, ...)
        {
            ft <- copy(ft)
            ft[, keep := TRUE]
            
            if (!is.null(absMinIntensity))
                ft[, keep := absIntPred(intensity)]
            
            if (!is.null(relMinIntensity))
            {
                maxInt <- max(ft$intensity)
                ft[keep == TRUE, keep := relIntPred(intensity, maxInt)]
            }
            
            if (!is.null(retentionRange))
                ft[keep == TRUE, keep := rangePred(ret, retentionRange)]
            
            if (!is.null(mzRange))
                ft[keep == TRUE, keep := rangePred(mz, mzRange)]
            
            if (!is.null(mzDefectRange))
                ft[keep == TRUE, keep := rangePred(mz - floor(mz), mzDefectRange)]
            
            if (!is.null(chromWidthRange))
                ft[keep == TRUE, keep := rangePred(retmax - retmin, chromWidthRange)]

            if (!is.null(qualityRange))
                ft[keep == TRUE, keep := scorePred(.SD, qualityRange), .SDcols = names(qualityRange)]
            
            return(!ft$keep)
        })
        
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
        x <- delete(x, setdiff(analyses(x), i))
    }

    return(x)
})

#' @describeIn features Extract a feature table for an analysis.
#' @export
setMethod("[[", c("features", "ANY", "missing"), function(x, i)
{
    assertExtractArg(i)
    return(x@features[[i]])
})

#' @describeIn features Extract a feature table for an analysis.
#' @export
setMethod("$", "features", function(x, name)
{
    eval(substitute(x@features$NAME_ARG, list(NAME_ARG = name)))
})

#' @templateVar where features
#' @templateVar what features
#' @template delete
#' @export
setMethod("delete", "features", function(obj, i = NULL, j = NULL, ...)
{
    ac <- checkmate::makeAssertCollection()
    i <- assertDeleteArgAndToChr(i, analyses(obj), add = ac)
    checkmate::assert(
        checkmate::checkIntegerish(j, any.missing = FALSE, null.ok = TRUE),
        checkmate::checkFunction(j, null.ok = TRUE),
        .var.name = "j"
    )
    checkmate::reportAssertions(ac)

    if (length(i) == 0 || (!is.null(j) && length(j) == 0))
      return(obj) # nothing to remove...
    
    # UNDONE: NULL for i and j will remove all?
    
    # i = NULL; j = vector: remove from all analyses
    # i = vector; j = NULL: remove specified analyses
    # j = function: remove specific features from given analyses (or all analyses if i=NULL)
    
    if (!is.function(j))
    {
        if (is.null(j))
            obj@features <- obj@features[setdiff(analyses(obj), i)]
        else
        {
            obj@features[i] <- lapply(obj@features[i], function(ft)
            {
                inds <- j[j <= nrow(ft)]
                return(if (length(inds) > 0) ft[-inds] else ft)
            })
        }
    }
    else
    {
        obj@features[i] <- Map(obj@features[i], i, f = function(ft, ana)
        {
            rm <- j(ft, ana, ...)
            if (is.logical(rm))
                return(ft[!rm])
            return(ft[setdiff(seq_len(nrow(ft)), rm)])
        })
    }
    obj@analysisInfo <- obj@analysisInfo[analysis %in% names(obj@features)]
    
    return(obj)
})

#' @describeIn features Calculates peak qualities for each feature. This uses
#'   \href{https://github.com/KelseyChetnik/MetaClean/}{MetaClean} \R package to calculate the following metrics:
#'   \verb{Apex-Boundary Ratio}, \verb{FWHM2Base}, \verb{Jaggedness}, \verb{Modality}, \verb{Symmetry}, \verb{Gaussian
#'   Similarity}, \verb{Sharpness}, \verb{Triangle Peak Area Similarity Ratio} and \verb{Zig-Zag index}. Please see the
#'   \pkg{MetaClean} publication (referenced below) for more details. For each metric, an additional score is calculated
#'   by normalizing all feature values (unless the quality metric definition has a fixed range) and scale from \samp{0}
#'   (worst) to \samp{1} (best). Then, a \verb{totalScore} for each feature is calculated by the (weighted) sum of all
#'   score values.
#'
#' @param weights A named \code{numeric} vector that defines the weight for each score to calculate the
#'   \verb{totalScore}. The names of the vector follow the score names. Unspecified weights are defaulted to \samp{1}.
#'   Example: \code{weights=c(ApexBoundaryRatioScore=0.5, GaussianSimilarityScore=2)}.
#' @param flatnessFactor Passed to \pkg{MetaClean} as the \code{flatness.factor} argument to
#'   \code{\link[MetaClean]{calculateJaggedness}} and \code{\link[MetaClean]{calculateModality}}.
#'
#' @template parallel-arg
#' 
#' @references \insertRef{Chetnik2020}{patRoon}
#'
#' @return \code{calculatePeakQualities} returns a modified object amended with peak qualities and scores.
#'
#' @note For \code{calculatePeakQualities}: sometimes \pkg{MetaClean} may return \code{NA} for the \verb{Gaussian
#'   Similarity} metric, in which case it will be set to \samp{0}.
#' @export
setMethod("calculatePeakQualities", "features", function(obj, weights, flatnessFactor, parallel = TRUE)
{
    checkPackage("MetaClean")
    
    if (length(obj) == 0)
        return(obj) # nothing to do...
    
    featQualities <- featureQualities()
    featQualityNames <- featureQualityNames(group = FALSE)
    featScoreNames <- featureQualityNames(group = FALSE, scores = TRUE, totScore = FALSE)
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertNumeric(weights, finite = TRUE, any.missing = FALSE, min.len = 1, names = "unique",
                             null.ok = TRUE, add = ac)
    if (!is.null(weights))
        checkmate::assertNames(names(weights), subset.of = featScoreNames, add = ac)
    checkmate::assertNumber(flatnessFactor, add = ac)
    checkmate::assertFlag(parallel, add = ac)
    checkmate::reportAssertions(ac)
    
    hash <- makeHash(obj, weights, flatnessFactor)
    cd <- loadCacheData("calculatePeakQualities", hash)
    if (!is.null(cd))
        return(cd)
    
    EICs <- getEICsForFeatures(obj)
    
    # HACK HACK HACK: MetaClean::calculateGaussianSimilarity needs to have
    # xcms::SSgauss attached
    # based on https://stackoverflow.com/a/36611896
    eg <- new.env()
    eg$SSgauss <- xcms::SSgauss
    withr::local_environment(eg)
    
    calcFeatQualities <- function(ret, retmin, retmax, intensity, EIC)
    {
        # NOTE: MetaClean expects matrices
        args <- list(c(rt = ret, rtmin = retmin, rtmax = retmax, maxo = intensity), as.matrix(EIC))
        return(sapply(featQualityNames, function(q)
        {
            a <- args
            if (q %in% c("Jaggedness", "Modality"))
                a <- c(a, flatnessFactor)
            qual <- do.call(featQualities[[q]]$func, a)
            if (q == "GaussianSimilarity" && is.na(qual))
                qual <- 0
            return(qual)
        }, simplify = FALSE))
    }
    
    doCalcs <- function(ft, eic)
    {
        ft <- copy(ft)
        
        if (nrow(ft) == 0)
            ft[, c(featQualityNames, featScoreNames) := numeric()]
        else
        {
            ft[, (featQualityNames) := rbindlist(Map(calcFeatQualities, ret, retmin, retmax, intensity, eic))]
            ft[, (featScoreNames) := Map(patRoon:::scoreFeatQuality, featQualities, .SD), .SDcols = featQualityNames]
        }
        patRoon:::doProgress()
        return(ft)
    }
    
    printf("Calculating feature peak qualities and scores...\n")

    if (parallel)
        fTable <- withProg(length(EICs), TRUE, future.apply::future_Map(doCalcs, featureTable(obj)[names(EICs)], EICs))
    else
        fTable <- withProg(length(EICs), FALSE, Map(doCalcs, featureTable(obj)[names(EICs)], EICs))
    
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
      
        if (parallel)
            setDT(ft)
        set(ft, j = "totalScore", value = rowSums(wft, na.rm = TRUE))
    })

    featureTable(obj) <- fTable
    
    saveCacheData("calculatePeakQualities", obj, hash)
    
    return(obj)
})


#' Finding features
#'
#' Automatically find features.
#'
#' Several functions exist to collect features (\emph{i.e.} retention and MS information that represent potential
#' compounds) from a set of analyses. All 'feature finders' return an object derived from the \code{\link{features}}
#' base class. The next step in a general workflow is to group and align these features across analyses with
#' \code{\link{groupFeatures}}. Note that some feature finders have a plethora of options which sometimes may have a
#' large effect on the quality of results. Fine-tuning parameters is therefore important, and the optimum is largely
#' dependent upon applied analysis methodology and instrumentation.
#'
#' @param verbose If set to \code{FALSE} then no text output is shown.
#' @param \dots Further parameters passed to the selected feature finding algorithms.
#'
#' @template analysisInfo-arg
#'
#' @templateVar func findFeatures
#' @templateVar what find features
#' @templateVar ex1 findFeaturesOpenMS
#' @templateVar ex2 findFeaturesXCMS
#' @templateVar algos bruker,openms,xcms,xcms3,envipick,sirius,kpic2,safd
#' @templateVar algosSuffix Bruker,OpenMS,XCMS,XCMS3,EnviPick,SIRIUS,KPIC2,SAFD
#' @templateVar ret features
#' @template generic-algo
#'
#' @template centroid_note
#'
#' @return An object of a class which is derived from \code{\link{features}}.
#'
#' @export
findFeatures <- function(analysisInfo, algorithm, ..., verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo, add = ac)
    checkmate::assertChoice(algorithm, c("bruker", "openms", "xcms", "xcms3", "envipick", "sirius", "kpic2", "safd"),
                            add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)

    f <- switch(algorithm,
                bruker = findFeaturesBruker,
                openms = findFeaturesOpenMS,
                xcms = findFeaturesXCMS,
                xcms3 = findFeaturesXCMS3,
                envipick = findFeaturesEnviPick,
                sirius = findFeaturesSIRIUS,
                kpic2 = findFeaturesKPIC2,
                safd = findFeaturesSAFD)

    f(analysisInfo, ..., verbose = verbose)
}

#' Import features
#'
#' Generic function to import features produced by other software.
#'
#' @templateVar func importFeatures
#' @templateVar what import features
#' @templateVar ex1 importFeaturesXCMS3
#' @templateVar ex2 importFeaturesKPIC2
#' @templateVar algosSuffix XCMS,XCMS3,KPIC2,EnviMass
#' @templateVar ret features
#' @templateVar noParam TRUE
#' @template generic-algo
#'
#' @template analysisInfo-arg
#' 
#' @param type What type of data should be imported: \code{"xcms"}, \code{"xcms3"}, \code{"kpic2"} or \code{"envimass"}.
#' @param \dots Further arguments passed to the selected import algorithm function.
#'
#' @inherit findFeatures return
#' 
#' @seealso \code{\link{findFeatures}} to find new features.
#'
#' @export
importFeatures <- function(analysisInfo, type, ...)
{
    analysisInfo <- assertAndPrepareAnaInfo(analysisInfo)

    f <- switch(type,
                xcms = importFeaturesXCMS,
                xcms3 = importFeaturesXCMS3,
                kpic2 = importFeaturesKPIC2,
                envimass = importFeaturesEnviMass,
                stop("Invalid algorithm! Should be: xcms, xcms3, kpic2 or envimass"))

    f(analysisInfo = analysisInfo, ...)
}
