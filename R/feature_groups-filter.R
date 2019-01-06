#' @include main.R
NULL

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
        ret <- func(fGroups)
        saveCacheData(cacheName, ret, hash)
    }

    if (verbose)
    {
        newn <- ncol(ret@groups)
        printf("Done! Filtered %d (%.2f%%) groups. Remaining: %d.\n", oldn - newn, (1-(newn/oldn))*100, newn)
    }

    return(ret)
}

intensityFilter <- function(fGroups, threshold, negate = FALSE)
{
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
            samplesWithRef <- which(sapply(refGroups, function(refs) ref %in% refs))
            refSamples <- which(anaInfo$group == ref)
            thr <- fGroups@groups[refSamples, lapply(.SD, function(x)
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

abundanceFilter <- function(fGroups, relThreshold = 0, absThreshold = 0, negate = FALSE)
{
    absThreshold <- max(absThreshold, relThreshold * nrow(analysisInfo(fGroups)))
    return(doFilter(fGroups, "abundance", c(absThreshold, negate), function(fGroups)
    {
        gTable <- groups(fGroups)
        pred <- function(x) sum(x > 0) >= absThreshold
        if (negate)
            pred <- Negate(pred)
        abundantGroups <- sapply(gTable, pred)
        return(fGroups[, abundantGroups])
    }))
}

interReplicateAbundanceFilter <- function(fGroups, relThreshold = 0, absThreshold = 0, negate = FALSE, verbose = TRUE)
{
    anaInfo <- analysisInfo(fGroups)
    rGroups <- unique(anaInfo$group)
    absThreshold <- max(absThreshold, relThreshold * length(rGroups))
    groupNames <- eval(colnames(fGroups@groups))

    return(doFilter(fGroups, "inter replicate abundance", c(absThreshold, negate), function(fGroups)
    {
        fGroups@groups[, group := fGroups@analysisInfo$group]

        pred <- function(x) length(unique(anaInfo$group[which(x > 0)])) >= absThreshold
        if (negate)
            pred <- Negate(pred)

        fGroups@groups[, (groupNames) := lapply(.SD, function(x) if (pred(x)) x else 0),
                       .SDcols = groupNames]

        fGroups@groups[, group := NULL]

        return(updateFeatIndex(removeEmptyGroups(fGroups)))
    }, "inter_rgroup_abundance", verbose))
}

intraReplicateFilter <- function(fGroups, threshold, negate = FALSE)
{
    groupNames <- colnames(fGroups@groups)

    return(doFilter(fGroups, "replicate", c(threshold, negate), function(fGroups)
    {
        # add groups temporarily
        fGroups@groups[, group := fGroups@analysisInfo$group]

        pred <- function(x, n) (sum(x > 0) / n) >= threshold
        if (negate)
            pred <- Negate(pred)

        fGroups@groups[, (groupNames) := lapply(.SD, function(x) if (pred(x, .N)) x else 0),
                       by = group, .SDcols = groupNames]
        fGroups@groups[, group := NULL]

        return(updateFeatIndex(removeEmptyGroups(fGroups)))
    }))
}

retentionMzFilter <- function(fGroups, range, negate, doRet)
{
    return(doFilter(fGroups, if (doRet) "retention" else "mz", c(range, negate), function(fGroups)
    {
        if (range[2] < 0)
            pred <- function(x) numGTE(x, range[1])
        else
            pred <- function(x) numGTE(x, range[1]) & numLTE(x, range[2])

        if (negate)
            pred <- Negate(pred)

        return(fGroups[, pred(if (doRet) fGroups@groupInfo$rts else fGroups@groupInfo$mzs)])
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
            if (range[2] < 0)
                return(numGTE(cwidths, range[1]))
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
#' @param intensityThreshold Minimum intensity. Set to 0 or \code{NULL} to skip
#'   this step.
#' @param
#' relAbundance,absAbundance,interRelRGroupAbundance,interAbsRGroupAbundance
#' Minimal overall relative/absolute abundance for a feature group to be present
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
#'   numeric vector with length of two containing the min/max values. If the max
#'   value is set to a value below 0 then no maximum is assumed. Set to
#'   \code{NULL} to skip this step.
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
setMethod("filter", "featureGroups", function(obj, intensityThreshold = NULL, relAbundance = NULL,
                                              absAbundance = NULL, interRelRGroupAbundance = NULL,
                                              interAbsRGroupAbundance = NULL, intraRGroupAbundance = NULL,
                                              minBlankThreshold = NULL, retentionRange = NULL, mzRange = NULL,
                                              chromWidthRange = NULL, rGroups = NULL, removeRefAnalyses = FALSE,
                                              repetitions = 2, negate = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertNumber, . ~ intensityThreshold + relAbundance + absAbundance +
               interRelRGroupAbundance + interAbsRGroupAbundance + intraRGroupAbundance +
               minBlankThreshold,
           lower = 0, finite = TRUE, null.ok = TRUE, fixed = list(add = ac))
    aapply(assertRange, . ~ retentionRange + mzRange + chromWidthRange, null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertCharacter(rGroups, min.chars = 1, min.len = 1, any.missing = FALSE, null.ok = TRUE, add = ac)
    checkmate::assertCount(repetitions, positive = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~ removeRefAnalyses + negate, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    if (length(obj) == 0)
        return(obj)
    
    if (!is.null(intensityThreshold) && intensityThreshold > 0)
        obj <- intensityFilter(obj, intensityThreshold)

    if (!is.null(retentionRange))
        obj <- retentionMzFilter(obj, retentionRange, negate, doRet = TRUE)

    if (!is.null(mzRange))
        obj <- retentionMzFilter(obj, mzRange, negate, doRet = FALSE)

    if (!is.null(chromWidthRange))
        obj <- chromWidthFilter(obj, chromWidthRange, negate)
    
    if (!is.null(rGroups))
        obj <- replicateGroupFilter(obj, rGroups, negate)

    # make sure that both absolute and relative abundances are set when only one of them is specified
    if (!is.null(relAbundance) && is.null(absAbundance))
        absAbundance <- 0
    if (is.null(relAbundance) && !is.null(absAbundance))
        relAbundance <- 0
    if (!is.null(interRelRGroupAbundance) && is.null(interAbsRGroupAbundance))
        interAbsRGroupAbundance <- 0
    if (is.null(interRelRGroupAbundance) && !is.null(interAbsRGroupAbundance))
        interRelRGroupAbundance <- 0

    if (!is.null(relAbundance) || !is.null(interRelRGroupAbundance) ||
        !is.null(intraRGroupAbundance) || !is.null(minBlankThreshold))
    {
        for (n in seq_len(repetitions))
        {
            if (repetitions > 1)
                printf("filter repetition %d/%d\n", n, repetitions)
            if (!is.null(relAbundance))
                obj <- abundanceFilter(obj, relAbundance, absAbundance, negate)
            if (!is.null(interRelRGroupAbundance))
                obj <- interReplicateAbundanceFilter(obj, interRelRGroupAbundance, interAbsRGroupAbundance, negate)
            if (!is.null(intraRGroupAbundance))
                obj <- intraReplicateFilter(obj, intraRGroupAbundance, negate)
            if (!is.null(minBlankThreshold))
                obj <- blankFilter(obj, minBlankThreshold, negate)
        }
    }

    if (removeRefAnalyses)
        obj <- replicateGroupFilter(obj, unique(analysisInfo(obj)$ref), negate = !negate)
    
    return(obj)
})

#' @details \code{replicateGroupSubtract} removes feature groups present in a
#'   given set of replicate groups (unless intensities are above a given
#'   threshold). The replicate groups that are subtracted will be removed.
#'
#' @param threshold Minimal relative threshold (compared to mean intensity of
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
