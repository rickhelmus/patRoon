#' @include main.R
NULL

getHighestAbsValue <- function(abs, rel, size)
{
    abs <- NULLToZero(abs); rel <- NULLToZero(rel)
    return(max(abs, rel * size))
}

#' Filtering of grouped features
#'
#' Basic rule based filtering of feature groups.
#'
#' @param fGroups,obj \code{\link{featureGroups}} object to which the filter is
#'   applied.
#' @param rGroups A character vector of replicate groups that should be kept
#'   (\code{filter}) or subtracted from (\code{replicateGroupSubtract}).
#'
#' @return A filtered \code{\link{featureGroups}} object. Feature groups that
#'   are filtered away have their intensity set to zero. In case a feature group
#'   is not present in any of the analyses anymore it will be removed
#'   completely.
#'
#' @name feature-filtering
#' @seealso \code{\link{featureGroups-class}}
#' @seealso \code{\link{feature-grouping}}
NULL

doFilter <- function(fGroups, what, hashParam, func, cacheCateg = what, verbose = TRUE)
{
    if (verbose)
    {
        printf("Applying %s filter... ", what)
        oldn <- ncol(fGroups@groups)
    }

    cacheName <- sprintf("filterFGroups_%s", cacheCateg)
    hash <- makeHash(fGroups, hashParam)
    ret <- loadCacheData(cacheName, hash)
    if (is.null(ret))
    {
        fGroups@groups <- copy(fGroups@groups)
        ret <- if (length(fGroups) > 0) func(fGroups) else fGroups
        saveCacheData(cacheName, ret, hash)
    }

    if (verbose)
    {
        newn <- ncol(ret@groups)
        printf("Done! Filtered %d (%.2f%%) groups. Remaining: %d.\n", oldn - newn,
               if (oldn > 0) (1-(newn/oldn))*100 else 0, newn)
    }

    return(ret)
}

intensityFilter <- function(fGroups, absThreshold, relThreshold, negate = FALSE)
{
    threshold <- getHighestAbsValue(absThreshold, relThreshold, max(groups(fGroups)))
    if (threshold == 0)
        return(fGroups)

    return(doFilter(fGroups, "intensity", c(threshold, negate), function(fGroups)
    {
        compF <- if (negate) function(x) x >= threshold else function(x) x < threshold

        # use set to speed stuff up: http://stackoverflow.com/a/20545629
        for (v in seq_along(fGroups@groups))
            set(fGroups@groups, which(compF(fGroups@groups[[v]])), v, 0)
        return(updateFeatIndex(removeEmptyGroups(fGroups)))
    }))
}

blankFilter <- function(fGroups, threshold, negate = FALSE)
{
    anaInfo <- analysisInfo(fGroups)
    gNames <- names(fGroups)
    rGroups <- unique(anaInfo$group)

    # multiple groups may be specified separated by comma
    refGroups <- sapply(anaInfo$ref, function(rg) strsplit(rg, ","), USE.NAMES = FALSE)
    allRefs <- unique(unlist(refGroups))
    allRefs <- allRefs[allRefs %in% rGroups]

    if (length(allRefs) == 0)
    {
        warning("No suitable reference analyses found, skipping blank filter...")
        return(fGroups)
    }

    return(doFilter(fGroups, "blank", c(threshold, negate), function(fGroups)
    {
        pred <- function(x, t) x >= t
        if (negate)
            pred <- Negate(pred)

        for (ref in allRefs)
        {
            refAnalyses <- which(anaInfo$group == ref)
            thr <- fGroups@groups[refAnalyses, lapply(.SD, function(x)
            {
                m <- mean(x[x > 0])
                if (is.na(m))
                    return(0)
                else
                    return(m * threshold)
            })]

            fGroups@groups[, (gNames) := lapply(seq_along(.SD), function(n) ifelse(pred(.SD[[n]], thr[[n]]), .SD[[n]], 0))]
        }

        return(updateFeatIndex(removeEmptyGroups(fGroups)))
    }))
}

minAnalysesFilter <- function(fGroups, absThreshold = 0, relThreshold = 0, negate = FALSE, verbose = TRUE)
{
    threshold <- getHighestAbsValue(absThreshold, relThreshold, length(analyses(fGroups)))
    if (threshold == 0)
        return(fGroups)
    return(doFilter(fGroups, "minimum analyses", c(threshold, negate), verbose = verbose, function(fGroups)
    {
        pred <- function(x) sum(x > 0) >= threshold
        if (negate)
            pred <- Negate(pred)
        return(fGroups[, sapply(groups(fGroups), pred, USE.NAMES = FALSE)])
    }, "minAnalyses"))
}

minReplicatesFilter <- function(fGroups, absThreshold = 0, relThreshold = 0, negate = FALSE, verbose = TRUE)
{
    threshold <- getHighestAbsValue(absThreshold, relThreshold, length(replicateGroups(fGroups)))
    if (threshold == 0)
        return(fGroups)

    rGroupsAna <- analysisInfo(fGroups)$group

    return(doFilter(fGroups, "minimum replicates", c(threshold, negate), function(fGroups)
    {
        pred <- function(x) length(unique(rGroupsAna[x > 0])) >= threshold
        if (negate)
            pred <- Negate(pred)

        return(fGroups[, sapply(groups(fGroups), pred, USE.NAMES = FALSE)])
    }, "minReplicates", verbose))
}

minFeaturesFilter <- function(fGroups, absThreshold = 0, relThreshold = 0, negate = FALSE, verbose = TRUE)
{
    threshold <- getHighestAbsValue(absThreshold, relThreshold, length(fGroups))
    if (threshold == 0)
        return(fGroups)

    return(doFilter(fGroups, "minimum features", c(threshold, negate), function(fGroups)
    {
        pred <- function(x) sum(x > 0) >= threshold
        if (negate)
            pred <- Negate(pred)

        return(fGroups[sapply(transpose(groups(fGroups)), pred, USE.NAMES = FALSE)])
    }, "minReplicates", verbose))
}

replicateAbundanceFilter <- function(fGroups, absThreshold, relThreshold, maxIntRSD, negate = FALSE)
{
    if (NULLToZero(absThreshold) == 0 && NULLToZero(relThreshold) == 0 && NULLToZero(maxIntRSD) == 0)
        return(fGroups) # all thresholds NULL/0

    gNames <- names(fGroups)
    rGroupsAna <- fGroups@analysisInfo$group

    doThr <- !is.null(absThreshold) || !is.null(relThreshold)
    if (doThr)
    {
        if (!is.null(relThreshold))
            thresholds <- sapply(replicateGroups(fGroups),
                                 function(rg) getHighestAbsValue(absThreshold, relThreshold, sum(rGroupsAna == rg)))
        else
            thresholds <- setNames(rep(absThreshold, length(replicateGroups(fGroups))), replicateGroups(fGroups))
    }

    return(doFilter(fGroups, "replicate abundance", c(absThreshold, relThreshold, maxIntRSD, negate), function(fGroups)
    {
        # add replicate groups temporarily
        fGroups@groups[, group := rGroupsAna]

        pred <- function(x, n, rg)
        {
            ret <- TRUE
            if (doThr)
                ret <- sum(x > 0) >= thresholds[[rg]]
            if (ret && NULLToZero(maxIntRSD) != 0 && any(x > 0))
                ret <- (sd(x) / mean(x)) < maxIntRSD # UNDONE: remove zero's?
            return(ret)
        }
        if (negate)
            pred <- Negate(pred)

        fGroups@groups[, (gNames) := lapply(.SD, function(x) if (pred(x, .N, group)) x else 0),
                       by = group, .SDcols = gNames]
        fGroups@groups[, group := NULL]

        return(updateFeatIndex(removeEmptyGroups(fGroups)))
    }, "replicateAbundance"))
}

retentionMzFilter <- function(fGroups, range, negate, what)
{
    return(doFilter(fGroups, what, c(range, negate), function(fGroups)
    {
        pred <- function(x) numGTE(x, range[1]) & numLTE(x, range[2])

        if (negate)
            pred <- Negate(pred)

        checkVals <- switch(what,
                            retention = fGroups@groupInfo$rts,
                            mz = fGroups@groupInfo$mzs,
                            mzDefect = fGroups@groupInfo$mzs - floor(fGroups@groupInfo$mzs))

        return(fGroups[, pred(checkVals)])
    }))
}

chromWidthFilter <- function(fGroups, range, negate)
{
    ftindex <- groupFeatIndex(fGroups)
    fTable <- featureTable(fGroups)
    anaInfo <- analysisInfo(fGroups)
    gNames <- names(fGroups)

    return(doFilter(fGroups, "chromwidth", c(range, negate), function(fGroups)
    {
        pred <- function(finds)
        {
            cwidths <- sapply(seq_along(finds), function(i)
            {
                if (finds[i] == 0)
                    0
                else
                    diff(unlist(fTable[[anaInfo$analysis[i]]][finds[i], c("retmin", "retmax")]))
            }, USE.NAMES = FALSE)
            return(numGTE(cwidths, range[1]) & numLTE(cwidths, range[2]))
        }

        if (negate)
            pred <- Negate(pred)

        fGroups@groups[, (gNames) := lapply(seq_along(.SD), function(n) ifelse(pred(ftindex[[n]]), .SD[[n]], 0))]

        return(updateFeatIndex(removeEmptyGroups(fGroups)))
    }))
}

replicateGroupFilter <- function(fGroups, rGroups, negate = FALSE, verbose = TRUE)
{
    return(doFilter(fGroups, "replicate group", c(rGroups, negate), function(fGroups)
    {
        pred <- function(g) g %in% rGroups
        if (negate)
            pred <- Negate(pred)

        fGroups <- removeAnalyses(fGroups, which(!pred(fGroups@analysisInfo$group)))
        return(removeEmptyGroups(fGroups))
    }, "replicate_group", verbose))
}

#' @details \code{filter} performs common rule based filtering of feature groups
#'   such as blank subtraction, minimum intensity and minimal replicate
#'   abundance.
#'
#' @param absMinIntensity Minimum intensity for a features to be kept . Set to 0
#'   or \code{NULL} to skip this step.
#' @param
#' relAbundance,absAbundance,interRelRGroupAbundance,interAbsRGroupAbundance
#' Minimum overall relative/absolute abundance for a feature group to be present
#' within all analyses/replicate groups.
#' @param intraRGroupAbundance Minimum relative abundance (0-1) that a feature
#'   group should be present within a replicate group. Set to \code{NULL} to
#'   skip this step.
#' @param minBlankThreshold Feature groups that are also present in reference
#'   analyses (see \link[=analysis-information]{analysis info}) are filtered out
#'   unless their relative intensity is above this threshold. Set to \code{NULL}
#'   to skip this step.
#' @param retentionRange,mzRange,chromWidthRange Range of retention time (in
#'   seconds), \emph{m/z} or chromatographic peak width (in seconds),
#'   respectively. Features outside this range will be removed. Should be a
#'   numeric vector with length of two containing the min/max values. \code{Inf}
#'   can be specified for no maximum range. Set to \code{NULL} to skip this
#'   step.
#' @param removeRefAnalyses Set to \code{TRUE} to remove all analyses that
#'   belong to replicate groups that are specified as a reference in the
#'   \link{analysis-information}. This is useful to simplify the analyses in the
#'   specified \code{\link{featureGroups}} object after blank subtraction.
#'   Consequently, when both \code{minBlankThreshold} and this argument are set,
#'   blank subtraction is performed prior to removing any analyses.
#' @param repetitions Sometimes more feature groups may be removed by repeating
#'   filtering steps after another. Usually a value of two is enough to filter
#'   the maximum amount of feature groups.
#' @param negate If set to \code{TRUE} then filtering operations are performed
#'   in opposite manner.
#'
#' @rdname feature-filtering
#' @export
setMethod("filter", "featureGroups", function(obj, absMinIntensity = NULL, relMinIntensity = NULL,
                                              preAbsMinIntensity = NULL, preRelMinIntensity = NULL,
                                              absMinAnalyses = NULL, relMinAnalyses = NULL,
                                              absMinReplicates = NULL, relMinReplicates = NULL,
                                              absMinFeatures = NULL, relMinFeatures = NULL,
                                              absMinReplicateAbundance = NULL, relMinReplicateAbundance = NULL,
                                              maxReplicateIntRSD = NULL, blankThreshold = NULL,
                                              retentionRange = NULL, mzRange = NULL, mzDefectRange = NULL,
                                              chromWidthRange = NULL, rGroups = NULL, removeBlanks = FALSE,
                                              negate = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertNumber, . ~ absMinIntensity + relMinIntensity + preAbsMinIntensity + preRelMinIntensity +
               absMinAnalyses + relMinAnalyses + absMinReplicates + relMinReplicates + absMinFeatures + relMinFeatures +
               absMinReplicateAbundance + relMinReplicateAbundance + maxReplicateIntRSD +
               blankThreshold,
           lower = 0, finite = TRUE, null.ok = TRUE, fixed = list(add = ac))
    aapply(assertRange, . ~ retentionRange + mzRange + mzDefectRange + chromWidthRange, null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertCharacter(rGroups, min.chars = 1, min.len = 1, any.missing = FALSE, null.ok = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ removeBlanks + negate, fixed = list(add = ac))
    checkmate::reportAssertions(ac)

    if (length(obj) == 0)
        return(obj)

    maybeDoFilter <- function(func, arg1, ..., otherArgs = list())
    {
        args <- c(list(arg1), ...)
        if (any(!sapply(args, is.null)))
            return(do.call(func, c(list(obj, arg1, ..., negate = negate), otherArgs)))
        return(obj)
    }

    obj <- maybeDoFilter(intensityFilter, preAbsMinIntensity, preRelMinIntensity)

    obj <- maybeDoFilter(retentionMzFilter, retentionRange, otherArgs = list(what = "retention"))
    obj <- maybeDoFilter(retentionMzFilter, mzRange, otherArgs = list(what = "mz"))
    obj <- maybeDoFilter(retentionMzFilter, mzDefectRange, otherArgs = list(what = "mzDefect"))
    obj <- maybeDoFilter(chromWidthFilter, chromWidthRange)

    # replicate round #1
    obj <- maybeDoFilter(replicateAbundanceFilter, absMinReplicateAbundance, relMinReplicateAbundance, maxReplicateIntRSD)
    lenAfter <- length(obj)

    obj <- maybeDoFilter(blankFilter, blankThreshold)
    obj <- maybeDoFilter(intensityFilter, absMinIntensity, relMinIntensity)

    # replicate round #2 (only do if previous filters affected results)
    if (length(obj) != lenAfter)
        obj <- maybeDoFilter(replicateAbundanceFilter, absMinReplicateAbundance, relMinReplicateAbundance, maxReplicateIntRSD)


    obj <- maybeDoFilter(minAnalysesFilter, absMinAnalyses, relMinAnalyses)
    obj <- maybeDoFilter(minReplicatesFilter, absMinReplicates, relMinReplicates)
    obj <- maybeDoFilter(minFeaturesFilter, absMinFeatures, relMinFeatures)

    obj <- maybeDoFilter(replicateGroupFilter, rGroups)
    if (removeBlanks)
        obj <- replicateGroupFilter(obj, unique(analysisInfo(obj)$ref), negate = !negate)

    return(obj)
})

#' @details \code{replicateGroupSubtract} removes feature groups present in a
#'   given set of replicate groups (unless intensities are above a given
#'   threshold). The replicate groups that are subtracted will be removed.
#'
#' @param threshold Minimum relative threshold (compared to mean intensity of
#'   replicate group being subtracted) for a feature group to be \emph{not}
#'   removed. When \samp{0} a feature group is always removed when present in
#'   the given replicate groups.
#'
#' @rdname feature-filtering
#' @aliases replicateGroupSubtract
#' @export
setMethod("replicateGroupSubtract", "featureGroups", function(fGroups, rGroups, threshold)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertCharacter(rGroups, min.chars = 1, add = ac)
    checkmate::assertNumber(threshold, lower = 0, finite = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    if (length(fGroups) == 0)
        return(fGroups)

    checkIntensities <- threshold > 0
    gNames <- names(fGroups)
    fGroups@groups <- copy(fGroups@groups)

    filteredGroups <- replicateGroupFilter(fGroups, rGroups, verbose = FALSE)
    sharedGroups <- gNames[gNames %in% names(filteredGroups)]

    if (length(sharedGroups) == 0)
        return(fGroups)

    if (checkIntensities)
    {
        avgGroups <- averageGroups(filteredGroups)
        thrs <- sapply(avgGroups, max) * threshold
    }

    for (b in sharedGroups)
    {
        if (checkIntensities)
            set(fGroups@groups, which(fGroups@groups[[b]] < thrs[b]), b, 0)
        else
            fGroups@groups[, (b) := 0] # no threshold, zero-out everything
    }

    return(replicateGroupFilter(updateFeatIndex(removeEmptyGroups(fGroups)), rGroups, negate = TRUE, verbose = FALSE))
})
