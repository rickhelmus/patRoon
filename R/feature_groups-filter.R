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
    if (length(fGroups) == 0)
        return(fGroups)
    
    threshold <- getHighestAbsValue(absThreshold, relThreshold, max(sapply(groupTable(fGroups), max)))
    if (threshold == 0)
        return(fGroups)

    return(doFilter(fGroups, "intensity", c(threshold, negate), function(fGroups)
    {
        compF <- if (negate) function(x) x >= threshold else function(x) x < threshold
        delGroups <- setnames(as.data.table(matrix(FALSE, length(analyses(fGroups)), length(fGroups))),
                              names(fGroups))
        delGroups[, (names(delGroups)) := lapply(fGroups@groups, compF), by = rep(1, nrow(delGroups))]
        return(delete(fGroups, j = delGroups))
    }))
}

blankFilter <- function(fGroups, threshold, negate = FALSE)
{
    anaInfo <- analysisInfo(fGroups)
    gNames <- names(fGroups)
    rGroups <- unique(anaInfo$group)

    # multiple groups may be specified separated by comma
    blankGroups <- sapply(anaInfo$blank, function(rg) strsplit(rg, ","), USE.NAMES = FALSE)
    allBlanks <- unique(unlist(blankGroups))
    allBlanks <- allBlanks[allBlanks %in% rGroups]
    blAnaInds <- anaInfo$group %chin% allBlanks
    
    if (length(allBlanks) == 0)
    {
        warning("No suitable blank analyses found, skipping blank filter...")
        return(fGroups)
    }

    return(doFilter(fGroups, "blank", c(threshold, negate), function(fGroups)
    {
        pred <- function(x, t) x < t
        if (negate)
            pred <- Negate(pred)

        avgBls <- lapply(allBlanks, function(bl)
        {
            avg <- vapply(fGroups@groups[anaInfo$group == bl], function(x) mean(x[x > 0]), FUN.VALUE = numeric(1),
                          USE.NAMES = FALSE)
            avg[is.na(avg)] <- 0
            return(avg)
        })
        avgBls <- transpose(avgBls)
        minInts <- sapply(avgBls, max) * threshold

        delGroups <- copy(fGroups@groups)
        
        for (j in seq_along(delGroups))
            set(delGroups, j = j, value = fifelse(pred(fGroups@groups[[j]], minInts[[j]]), 1, 0))
        return(delete(fGroups, j = delGroups))
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
        return(fGroups[, sapply(groupTable(fGroups), pred, USE.NAMES = FALSE)])
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

        return(fGroups[, sapply(groupTable(fGroups), pred, USE.NAMES = FALSE)])
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

        return(fGroups[sapply(transpose(groupTable(fGroups)), pred, USE.NAMES = FALSE)])
    }, "minReplicates", verbose))
}

replicateAbundanceFilter <- function(fGroups, absThreshold, relThreshold, maxIntRSD, negate = FALSE)
{
    if (NULLToZero(absThreshold) == 0 && NULLToZero(relThreshold) == 0 && NULLToZero(maxIntRSD) == 0)
        return(fGroups) # all thresholds NULL/0

    gNames <- names(fGroups)
    rGroupsAna <- fGroups@analysisInfo$group
    rGroups <- replicateGroups(fGroups)
    rGroupLens <- table(rGroupsAna)
    rGroupInds <- sapply(rGroups, function(rg) which(rGroupsAna == rg), simplify = FALSE)

    doThr <- !is.null(absThreshold) || !is.null(relThreshold)
    if (doThr)
    {
        if (!is.null(relThreshold))
            thresholds <- sapply(replicateGroups(fGroups),
                                 function(rg) getHighestAbsValue(absThreshold, relThreshold, sum(rGroupsAna == rg)))
        else
            thresholds <- setNames(rep(absThreshold, length(replicateGroups(fGroups))), replicateGroups(fGroups))
    }
    
    maxIntRSD <- NULLToZero(maxIntRSD)

    return(doFilter(fGroups, "replicate abundance", c(absThreshold, relThreshold, maxIntRSD, negate), function(fGroups)
    {
        pred <- function(x, n, rg)
        {
            if (doThr && sum(x > 0) < thresholds[[rg]])
                return(TRUE)
            return(maxIntRSD != 0 && length(x) > 1 && any(x > 0) && (sd(x) / mean(x)) > maxIntRSD) # UNDONE: remove zeros?
        }
        
        if (negate)
            pred <- Negate(pred)

        delGroups <- copy(fGroups@groups)
        set(delGroups, j = "group", value = rGroupsAna)
        delGroups[, (gNames) := lapply(.SD, function(x) if (pred(x, .N, group)) 1 else 0), by = group, .SDcols = gNames]
        return(delete(fGroups, j = delGroups[, -"group"]))
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
    anas <- analyses(fGroups)
    gNames <- names(fGroups)

    return(doFilter(fGroups, "chromwidth", c(range, negate), function(fGroups)
    {
        pred <- function(finds)
        {
            cwidths <- sapply(seq_along(finds), function(i)
            {
                if (finds[i] == 0)
                    return(0)
                else
                    return(fTable[[anas[i]]][["retmax"]][finds[i]] - fTable[[anas[i]]][["retmin"]][finds[i]])
            }, USE.NAMES = FALSE)
            return(cwidths < range[1] | cwidths > range[2])
        }
        
        if (negate)
            pred <- Negate(pred)
        
        delGroups <- setnames(as.data.table(matrix(FALSE, length(analyses(fGroups)), length(fGroups))),
                              names(fGroups))
        delGroups[, (names(delGroups)) := lapply(ftindex, pred), by = rep(1, nrow(delGroups))]
        return(delete(fGroups, j = delGroups))
    }))
}

replicateGroupFilter <- function(fGroups, rGroups, negate = FALSE, verbose = TRUE)
{
    return(doFilter(fGroups, "replicate group", c(rGroups, negate), function(fGroups)
    {
        pred <- function(g) !g %chin% rGroups
        if (negate)
            pred <- Negate(pred)
        return(delete(fGroups, pred(analysisInfo(fGroups)$group)))
    }, "replicate_group", verbose))
}

checkFeaturesFilter <- function(fGroups, checkFeaturesSession, negate)
{
    sessionPath <- getCheckFeaturesSessionPath(checkFeaturesSession)
    return(doFilter(fGroups, "checked features session", c(makeFileHash(sessionPath), negate), function(fGroups)
    {
        session <- readRDS(sessionPath)
        if (negate)
            fGroups <- fGroups[, setdiff(names(fGroups), session$enabledFGroups)]
        else
            fGroups <- fGroups[, session$enabledFGroups]
        
        gNames <- names(fGroups)
        pred <- if (negate) function(ef) !ef else function(ef) ef
        fGroups@groups[, (gNames) := Map(.SD, session$enabledFeatures[gNames],
                                         f = function(fg, ef) ifelse(pred(ef), fg, 0))]
        return(cleanGroups(fGroups, TRUE))
    }, "checkedFeatures"))
}

#' @details \code{filter} performs common rule based filtering of feature groups
#'   such as blank subtraction, minimum intensity and minimum replicate
#'   abundance. Removing of features occurs by zeroing their intensity values.
#'   Furthermore, feature groups that are left completely empty (\emph{i.e.} all
#'   intensities are zero) will be automatically removed.
#'
#' @param preAbsMinIntensity,preRelMinIntensity As
#'   \code{absMinIntensity}/\code{relMinIntensity}, but applied \emph{before}
#'   any other filters. This is typically used to speed-up subsequent filter
#'   steps. However, care must be taken that a sufficiently low value is choosen
#'   that is not expected to affect subsequent filtering steps. See below why
#'   this may be important.
#' @param absMinAnalyses,relMinAnalyses Feature groups are only kept when they
#'   contain data for at least this (absolute or relative) amount of analyses.
#'   Set to \code{NULL} to ignore.
#' @param absMinReplicates,relMinReplicates Feature groups are only kept when
#'   they contain data for at least this (absolute or relative) amount of
#'   replicates. Set to \code{NULL} to ignore.
#' @param absMinFeatures,relMinFeatures Analyses are only kept when they contain
#'   at least this (absolute or relative) amount of features. Set to \code{NULL}
#'   to ignore.
#' @param absMinReplicateAbundance,relMinReplicateAbundance Minimum
#'   absolute/relative abundance that a grouped feature should be present within
#'   a replicate group. If this mimimum is not met all features within the
#'   replicate group are removed. Set to \code{NULL} to skip this step.
#' @param maxReplicateIntRSD Maximum relative standard deviation (RSD) of
#'   intensity values for features within a replicate group. If the RSD is above
#'   this value all features within the replicate group are removed. Set to
#'   \code{NULL} to ignore.
#' @param blankThreshold Feature groups that are also present in blank analyses
#'   (see \link[=analysis-information]{analysis info}) are filtered out unless
#'   their relative intensity is above this threshold. For instance, a value of
#'   \samp{5} means that only features with an intensity five times higher than
#'   that of the blank are kept. The relative intensity values between blanks
#'   and non-blanks are determined from the mean of all non-zero blank
#'   intensities. Set to \code{NULL} to skip this step.
#' @param removeBlanks Set to \code{TRUE} to remove all analyses that belong to
#'   replicate groups that are specified as a blank in the
#'   \link{analysis-information}. This is useful to simplify the analyses in the
#'   specified \code{\link{featureGroups}} object after blank subtraction. When
#'   both \code{blankThreshold} and this argument are set, blank subtraction is
#'   performed prior to removing any analyses.
#'
#' @templateVar feat FALSE
#' @template feat-filter-args
#'
#' @section Filter order: When multiple arguments are specified to
#'   \code{filter}, multiple filters are applied in sequence. Since some of
#'   these filters may affect each other, choosing their order correctly may be
#'   important for effective data filtering. For instance, when an intensity
#'   filter removes features from blank analyses, a subsequent blank filter may
#'   not adequately perform blank subtraction. Similarly, when intensity and
#'   blank filters are executed after the replicate abundance filter it may be
#'   necessary to ensure minimum replicate abundance again as the intensity and
#'   blank filters may have removed some features within a replicate group.
#'
#'   With this in mind, filters (if specified) occur in the following order:
#'
#'   \enumerate{
#'
#'   \item Pre-Intensity filters (\emph{i.e.} \code{preAbsMinIntensity} and
#'   \code{preRelMinIntensity}).
#'
#'   \item Chromatography and mass filters (\emph{i.e} \code{retentionRange},
#'   \code{mzRange}, \code{mzDefectRange} and \code{chromWidthRange}).
#'
#'   \item Replicate abundance filters (\emph{i.e.}
#'   \code{absMinReplicateAbundance}, \code{relMinReplicateAbundance} and
#'   \code{maxReplicateIntRSD}).
#'
#'   \item Blank filter (\emph{i.e.} blankThreshold).
#'
#'   \item Intensity filters (\emph{i.e.} \code{absMinIntensity} and
#'   \code{relMinIntensity}).
#'
#'   \item Replicate abundance filters (2nd time, only if previous filters
#'   affected results).
#'
#'   \item General abundance filters (\emph{i.e.} \code{absMinAnalyses},
#'   \code{relMinAnalyses}, \code{absMinReplicates}, \code{relMinReplicates},
#'   \code{absMinFeatures} and \code{relMinFeatures}).
#'
#'   \item Replicate group filter (\emph{i.e.} \code{rGroups}) and blank
#'   analyses removal (\emph{i.e.} if \code{removeBlanks=TRUE}).
#'
#'   }
#'
#'   If another filtering order is desired then \code{filter} should be called
#'   multiple times with only one filter argument at a time.
#'
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
                                              checkFeaturesSession = NULL, negate = FALSE)
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
    assertCheckFeaturesSession(checkFeaturesSession, obj, mustExist = TRUE, canClearSession = FALSE, didClearSession = FALSE,
                               null.ok = TRUE, add = ac)
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

    obj <- maybeDoFilter(checkFeaturesFilter, checkFeaturesSession)
    
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
        obj <- replicateGroupFilter(obj, unique(analysisInfo(obj)$blank), negate = !negate)

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
    sharedGroups <- intersect(gNames, names(filteredGroups))

    if (length(sharedGroups) == 0)
        return(fGroups)

    if (checkIntensities)
    {
        avgGroups <- averageGroups(filteredGroups)
        thrs <- sapply(avgGroups, max) * threshold
    }

    if (!checkIntensities)
        fGroups <- delete(fGroups, j = sharedGroups)
    else
    {
        delGroups <- setnames(as.data.table(matrix(FALSE, length(analyses(fGroups)), length(fGroups))),
                              names(fGroups))
        delGroups[, (sharedGroups) := Map(fGroups@groups[, sharedGroups, with = FALSE], sharedGroups,
                                          f = function(x, grp) x < thrs[grp]),
                  by = rep(1, nrow(delGroups))]
        # fGroups <- delete(fGroups, j = function(x, grp) grp %chin% sharedGroups & x < thrs[grp])
        fGroups <- delete(fGroups, j = delGroups)
    }
    
    return(replicateGroupFilter(fGroups, rGroups, negate = TRUE, verbose = FALSE))
})
