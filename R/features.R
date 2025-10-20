# SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

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
#' This class provides a way to store intensity, retention times, \emph{m/z} and other data for all features in a set of
#' analyses. The class is \code{virtual} and derived objects are created by 'feature finders' such as
#' \code{findFeaturesOpenMS}, \code{findFeaturesXCMS} and \code{findFeaturesBruker}.
#'
#' @param obj,x,object \code{features} object to be accessed
#'
#' @seealso \code{\link{findFeatures}}
#'
#' @slot features List of features per analysis file. Use the \code{featureTable} method for access.
#' @slot analysisInfo A \code{data.table} with the \link[=analysis-information]{analysis information}. Use the
#'   \code{analysisInfo} method for access.
#' @slot hasMobilities A \code{logical} that is \code{TRUE} if the features object contain mobility information. Use the
#'   \code{hasMobilities} method for access.
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
#' @param retMin Plot retention time in minutes (instead of seconds).
#' @param title Character string used for title of the plot. If \code{NULL} a title will be automatically generated.
#' @param showLegend Plot a legend if TRUE.
#' @template plot-lim
#' @param \dots For \code{delete}: passed to the function specified as \code{j}.
#'
#'   For \code{plotTICs} and \code{plotBPCs}: further arguments passed to \code{\link[graphics]{plot}}.
#'
#'   \setsPassedArgs1{features}
#'
#' @author Rick Helmus <\email{r.helmus@@uva.nl}> and Ricardo Cunha <\email{cunha@@iuta.de}> (\code{getTICs},
#'   \code{getBPCs}, \code{plotTICs} and \code{plotBPCs} functions)
#'
#' @templateVar class features
#' @template class-hierarchy
#'
#' @export
features <- setClass("features",
                     slots = c(features = "list", analysisInfo = "data.table", hasMobilities = "logical"),
                     contains = c("VIRTUAL", "workflowStep"))

setMethod("initialize", "features", function(.Object, ...)
{
    .Object <- callNextMethod(.Object, ...)
    .Object@features <- makeEmptyListNamed(.Object@features)
    if (length(.Object@hasMobilities) == 0)
        .Object@hasMobilities <- FALSE # initialize
    return(.Object)
})

setMethod("clearMobilities", "features", function(obj)
{
    if (!hasMobilities(obj))
        return(obj)
    
    featureTable(obj) <- lapply(featureTable(obj), function(ft)
    {
        ft <- copy(ft)
        ft <- removeDTColumnsIfPresent(ft, c("mobility", "CCS", "ims_parent_ID", "mobmin", "mobmax"))
        mob_cols <- grep("^mob_", names(ft), value = TRUE)
        if (length(mob_cols) > 0)
            ft <- ft[, (mob_cols) := NULL]
        return(ft)
    })
    
    obj@hasMobilities <- FALSE
    return(obj)
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
    printf("Has IMS data: %s\n", if (hasMobilities(object)) "yes" else "no")
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

setMethod("reorderAnalyses", "features", function(obj, anas)
{
    anaInfo <- analysisInfo(obj)
    order <- match(anas, anaInfo$analysis)
    stopifnot(length(order) == nrow(anaInfo) && !anyNA(order))
    
    obj@analysisInfo <- anaInfo[order]
    featureTable(obj) <- featureTable(obj)[order]
    
    return(obj)
})

#' @describeIn features Get analysis information
#' @param df If \code{TRUE} then the returned value is a \code{data.frame}, otherwise a \code{data.table}.
#' @return \code{analysisInfo}: The \link[=analysis-information]{analysis information} of this \code{features} object.
#' @export
setMethod("analysisInfo", "features", function(obj, df = FALSE)
{
    checkmate::assertFlag(df)
    return(if (df) as.data.frame(obj@analysisInfo) else obj@analysisInfo)
})

#' @templateVar class features
#' @template analysisInfo-set
#' @param value A \code{data.frame} or \code{data.table} with the new analysis information.
#' @export
setReplaceMethod("analysisInfo", "features", function(obj, value)
{
    checkmate::assertDataFrame(value)
    
    oldAnaInfo <- analysisInfo(obj)
    
    if (nrow(oldAnaInfo) > nrow(value))
    {
        stop("Cannot remove analyses when changing the analysis information. Please subset the object instead.",
             call. = FALSE)
    }
    if (nrow(oldAnaInfo) < nrow(value))
        stop("Cannot add analyses.", call. = FALSE)
    if (!setequal(oldAnaInfo$analysis, value$analysis))
        stop("Cannot modify analysis column.", call. = FALSE)

    if (!checkmate::testNames(value$analysis, identical.to = oldAnaInfo$analysis)) # re-ordered?
        obj <- reorderAnalyses(obj, value$analysis)
    
    obj@analysisInfo <- assertAndPrepareAnaInfo(value)
    
    return(obj)
})

#' @templateVar class features
#' @templateVar what analyses
#' @template strmethod
#' @export
setMethod("analyses", "features", function(obj) analysisInfo(obj)$analysis)

#' @templateVar class features
#' @templateVar what replicates
#' @template strmethod
#' @export
setMethod("replicates", "features", function(obj) unique(analysisInfo(obj)$replicate))

#' @describeIn features Returns \code{TRUE} if the features object has mobility information.
#' @export
setMethod("hasMobilities", "features", function(obj) obj@hasMobilities)

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
#' @template IMSRangeParams-arg
#' @export
setMethod("filter", "features", function(obj, absMinIntensity = NULL, relMinIntensity = NULL,
                                         retentionRange = NULL, mzRange = NULL, mzDefectRange = NULL,
                                         chromWidthRange = NULL, IMSRangeParams = NULL, qualityRange = NULL,
                                         negate = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertNumber, . ~ absMinIntensity + relMinIntensity, lower = 0, finite = TRUE,
           null.ok = TRUE, fixed = list(add = ac))
    aapply(assertRange, . ~ retentionRange + mzRange + mzDefectRange + chromWidthRange, null.ok = TRUE,
           fixed = list(add = ac))
    assertIMSRangeParams(IMSRangeParams, null.ok = TRUE, add = ac)
    assertScoreRange(qualityRange, c(featureQualityNames(group = FALSE),
                                     featureQualityNames(group = FALSE, scores = TRUE)), add = ac)
    checkmate::assertFlag(negate, add = ac)
    checkmate::reportAssertions(ac)

    if (!is.null(IMSRangeParams) && !hasMobilities(obj))
        stop("Cannot apply IMS Range filter: no mobilities assigned", call. = FALSE)
    
    if (length(obj) == 0)
        return(obj)

    oldn <- length(obj)

    hash <- makeHash(obj, absMinIntensity, relMinIntensity, retentionRange, mzRange, mzDefectRange, chromWidthRange,
                     IMSRangeParams, qualityRange, negate)
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

            if (!is.null(IMSRangeParams))
            {
                if (is.null(ft[["CCS"]]) && IMSRangeParams$param == "CCS")
                    stop("Cannot apply IMS Range filter: no CCS values assigned", call. = FALSE)
                
                vals <- ft[[IMSRangeParams$param]]
                if (IMSRangeParams$mzRelative)
                    vals <- vals / ft$mz
                
                ft[keep == TRUE, keep := rangePred(vals, c(IMSRangeParams$lower, IMSRangeParams$upper))]
            }
            
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
#' @templateVar ex fList
#' @template feat-subset
#' @export
setMethod("[", c("features", "ANY", "missing", "missing"), function(x, i, j, ..., ni, reorder = FALSE, drop = TRUE)
{
    checkmate::assertFlag(reorder)
    env <- parent.frame()
    x <- doSubsetFeaturesByAna(x, i, ni, reorder = reorder, env = env)
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
#' @template EICParams-arg
#' @template parallel-arg
#'
#' @templateVar what \code{calculatePeakQualities} and TIC/BPC related functions
#' @template uses-msdata
#'
#' @references \insertRef{Chetnik2020}{patRoon}
#'
#' @return \code{calculatePeakQualities} returns a modified object amended with peak qualities and scores.
#'
#' @note For \code{calculatePeakQualities}: sometimes \pkg{MetaClean} may return \code{NA} for the \emph{Gaussian
#'   Similarity} metric, in which case it will be set to \samp{0}.
#' @export
setMethod("calculatePeakQualities", "features", function(obj, weights, flatnessFactor,
                                                         EICParams = getDefEICParams(window = 0), parallel = TRUE)
{
    checkPackage("MetaClean")
    
    if (length(obj) == 0)
        return(obj) # nothing to do...
    
    featScoreNames <- featureQualityNames(group = FALSE, scores = TRUE, totScore = FALSE)
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertNumeric(weights, finite = TRUE, any.missing = FALSE, min.len = 1, names = "unique",
                             null.ok = TRUE, add = ac)
    if (!is.null(weights))
        checkmate::assertNames(names(weights), subset.of = featScoreNames, add = ac)
    checkmate::assertNumber(flatnessFactor, add = ac)
    assertEICParams(EICParams, add = ac)
    checkmate::assertFlag(parallel, add = ac)
    checkmate::reportAssertions(ac)
    
    hash <- makeHash(obj, weights, flatnessFactor)
    cd <- loadCacheData("calculatePeakQualities", hash)
    if (!is.null(cd))
        return(cd)
    
    EICs <- getFeatureEIXs(obj, "EIC", EIXParams = EICParams)
    
    # HACK HACK HACK: MetaClean::calculateGaussianSimilarity needs to have
    # xcms::SSgauss attached
    # based on https://stackoverflow.com/a/36611896
    eg <- new.env()
    eg$SSgauss <- xcms::SSgauss
    withr::local_environment(eg)
    
    doCalcs <- function(ft, eic, flatf)
    {
        featQualities <- featureQualities()
        featQualityNames <- featureQualityNames(group = FALSE)
        featScoreNames <- featureQualityNames(group = FALSE, scores = TRUE, totScore = FALSE)
        
        ft <- copy(ft)
        
        if (nrow(ft) == 0)
            ft[, c(featQualityNames, featScoreNames) := numeric()]
        else
        {
            eicat <- attr(eic, "allXValues")
            ft[, (featQualityNames) := rbindlist(Map(patRoon:::calcFeatQualities, ret, retmin, retmax, intensity, eic,
                                                     MoreArgs = list(eicat, flatf)))]
            ft[, (featScoreNames) := Map(patRoon:::scoreFeatQuality, featQualities, .SD), .SDcols = featQualityNames]
        }
        return(ft)
    }
    
    printf("Calculating feature peak qualities and scores...\n")

    fTable <- doMap(parallel, featureTable(obj)[names(EICs)], EICs, f = doCalcs, MoreArgs = list(flatnessFactor))
    
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

#' @describeIn features Obtain the total ion chromatogram/s (TICs) of the analyses.
#' @inheritParams getTICs,data.frame-method
#' @export
setMethod("getTICs", "features", function(obj, retentionRange = NULL, MSLevel = 1)
{
    getTICs(analysisInfo(obj), retentionRange, MSLevel)
})


#' @describeIn features Obtain the base peak chromatogram/s (BPCs) of the analyses.
#' @export
setMethod("getBPCs", "features", function(obj, retentionRange = NULL, MSLevel = 1)
{
    getBPCs(analysisInfo(obj), retentionRange, MSLevel)
})


#' @describeIn features Plots the TICs of the analyses.
#' @export
setMethod("plotTICs", "features", function(obj, retentionRange = NULL, MSLevel = 1, retMin = FALSE, title = NULL, 
                                           groupBy = NULL, showLegend = TRUE, xlim = NULL,  ylim = NULL, ...)
{
    plotTICs(analysisInfo(obj), retentionRange, MSLevel, retMin, title, groupBy, showLegend, xlim, ylim, ...)
})


#' @describeIn features Plots the BPCs of the analyses.
#' @export
setMethod("plotBPCs", "features", function(obj, retentionRange = NULL, MSLevel = 1, retMin = FALSE, title = NULL,
                                           groupBy = NULL, showLegend = TRUE, xlim = NULL,  ylim = NULL, ...)
{
    plotBPCs(analysisInfo(obj), retentionRange, MSLevel, retMin, title, groupBy, showLegend, xlim, ylim, ...)
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
#' @templateVar algos bruker,openms,xcms,xcms3,envipick,sirius,kpic2,safd,piek
#' @templateVar algosSuffix Bruker,OpenMS,XCMS,XCMS3,EnviPick,SIRIUS,KPIC2,SAFD,Piek
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
    checkmate::assertChoice(algorithm, c("bruker", "openms", "xcms", "xcms3", "envipick", "sirius", "kpic2", "safd",
                                         "piek"), add = ac)
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
                safd = findFeaturesSAFD,
                piek = findFeaturesPiek)

    f(analysisInfo, ..., verbose = verbose)
}

#' Import features
#'
#' Generic function to import features produced by other software.
#'
#' @templateVar func importFeatures
#' @templateVar what import features
#' @templateVar ex1 importFeaturesXCMS3
#' @templateVar ex2 importFeaturesTable
#' @templateVar algosSuffix XCMS,XCMS3,KPIC2,Table,EnviMass
#' @templateVar ret features
#' @templateVar noParam TRUE
#' @template generic-algo
#'
#' @param input The input object or path that should be imported. See the algorithm specific functions for more details.
#' @param type What type of data should be imported: \code{"xcms"}, \code{"xcms3"}, \code{"kpic2"}, \code{"table"}, or
#'   \code{"envimass"}.
#' @param \dots Further arguments passed to the selected import algorithm function.
#'
#' @inherit findFeatures return
#'
#' @seealso \code{\link{findFeatures}} to find new features.
#'
#' @export
importFeatures <- function(input, type, ...)
{
    f <- switch(type,
                xcms = importFeaturesXCMS,
                xcms3 = importFeaturesXCMS3,
                kpic2 = importFeaturesKPIC2,
                table = importFeaturesTable,
                envimass = importFeaturesEnviMass,
                stop("Invalid algorithm! Should be: xcms, xcms3, kpic2, table or envimass"))

    f(input, ...)
}
