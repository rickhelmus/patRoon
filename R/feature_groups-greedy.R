# SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
NULL

#' @rdname featureGroups-class
#' @export
featureGroupsGreedy <- setClass("featureGroupsGreedy", contains = "featureGroups")

setMethod("initialize", "featureGroupsGreedy",
          function(.Object, ...) callNextMethod(.Object, algorithm = "greedy", ...))

#' Group features using greedy algorithm
#'
#' Group features using a greedy algorithm that maximizes group scores based on retention time, m/z, mobility, and
#' intensity similarities.
#'
#' @templateVar algo greedy
#' @templateVar do group features
#' @templateVar generic groupFeatures
#' @templateVar algoParam greedy
#' @template algo_generator
#'
#' @details The \code{greedy} algorithm is a simple feature grouping algorithm that can work with both HRMS and IMS-HRMS
#'   data. The algorithm groups features by iteratively building the best possible groups. Features are processed in
#'   order of decreasing intensity. For each feature, candidate groups are formed from all other (ungrouped) features
#'   within the specified retention time, \emph{m/z} and mobility windows. Each candidate group only contains a maximum
#'   of one feature per analysis. The candidates are then scored and the group with the lowest overall variations in
#'   retention time, \emph{m/z}, mobility and replicate intensity is then selected. This process is repeated until all
#'   features have been assigned to a group. The weights for each of the scoring terms can be configured.
#'
#' @template feat-arg
#'
#' @param rtalign Not yet supported. Provided for consistency with other grouping methods.
#' @param rtWindow,mzWindow,IMSWindow Numeric tolerances for retention time (seconds), \emph{m/z}, and mobility,
#'   respectively. The scoring terms are normalized to these values. Defaults to \code{defaultLim("retention",
#'   "medium")}, \code{defaultLim("mz", "medium")}, and \code{defaultLim("mobility", "medium")}, respectively (see \link{limits}).
#' @param scoreWeights Numeric vector specifying the scoring weights. Should contain the following named elements:
#'   \code{"retention"}, \code{"mz"}, \code{"mobility"}, and \code{"intensity"}.
#'
#' @inheritParams groupFeatures
#'
#' @inherit groupFeatures return
#'
#' @note Any links between IMS parents and mobility features are removed. This can occur \emph{e.g.} when \code{greedy}
#' is used to generate a \link[=featureGroups-compare]{feature consensus} from a \link[=assignMobilities_feat]{post
#' mobility assignment} workflow.
#'
#' @templateVar what groupFeaturesGreedy
#' @templateVar cl features
#' @template main-rd-method
#' @export
setMethod("groupFeaturesGreedy", "features", function(feat, rtalign = FALSE,
                                                      rtWindow = defaultLim("retention", "medium"),
                                                      mzWindow = defaultLim("mz", "medium"),
                                                      IMSWindow = defaultLim("mobility", "medium"),
                                                      scoreWeights = c(retention = 1, mz = 1, mobility = 1, intensity = 1),
                                                      verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFlag(rtalign, add = ac)
    aapply(checkmate::assertNumber, . ~ rtWindow + mzWindow + IMSWindow, lower = 0, finite = TRUE,
           fixed = list(add = ac))
    checkmate::assertNumeric(scoreWeights, len = 4, lower = 0, finite = TRUE, any.missing = FALSE, add = ac)
    assertHasNames(scoreWeights, c("retention", "mz", "mobility", "intensity"), add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    if (rtalign)
        stop("Retention time alignment (rtalign=TRUE) is not yet supported for greedy grouping!", call. = FALSE)
    
    return(doGroupFeatures(feat, doGroupFeaturesGreedy, "greedy", rtWindow, mzWindow, IMSWindow, scoreWeights,
                           verbose = verbose))
})

doGroupFeaturesGreedy <- function(feat, rtWindow, mzWindow, IMSWindow, scoreWeights, verbose)
{
    hash <- makeHash(feat, rtWindow, mzWindow, IMSWindow, scoreWeights)
    cd <- loadCacheData("groupFeaturesGreedy", hash)
    if (!is.null(cd))
        return(cd)
    
    anaInfo <- analysisInfo(feat)
    hasMob <- hasMobilities(feat)
    sqeps <- sqrt(.Machine$double.eps) # cache for numLTE below
    
    fTable <- as.data.table(feat)
    fTable[, groupID := NA_integer_]
    fTable[, anaRow := seq_len(.N), by = "analysis"]
    fTable[, row := seq_len(.N)]
    
    # HACK: add dummy mobilities so grouping code can handle features without mobilities
    if (!hasMob)
    {
        fTable[, mobility := -100] # some random value so mobilities will not influence grouping
        scoreWeights["mobility"] <- 0 # do not use mobility in scoring
    }
    else
    {
        fTable[, mobility_orig := mobility] # we need the original data later
        
        if (any(is.na(fTable$mobility)))
        {
            unreasonableMob <- if (all(is.na(fTable$mobility)))
                -100
            else
                max(fTable$mobility, na.rm = TRUE) + (10 * IMSWindow)
            
            fTable[is.na(mobility), mobility := unreasonableMob]
        }
    }

    if (verbose)
        printf("Grouping %d features... ", nrow(fTable))
    
    fTable[, replicate := anaInfo$replicate[match(analysis, anaInfo$analysis)]]
    reps <- replicates(feat)
    fTable[, groupID := getGroupIDs(ret, mz, mobility, intensity, match(analysis, anaInfo$analysis),
                                    match(replicate, reps), rtWindow, mzWindow, IMSWindow, scoreWeights)]

    if (verbose)
        printf("Done! Found %d groups.\n", if (nrow(fTable) > 0) max(fTable$groupID) else 0L)
    
    gInfo <- fTable[, .(ret = mean(ret), mz = mean(mz)), by = "groupID"]
    if (hasMob)
    {
        gInfo[, mobility := fTable[, mean(mobility_orig), by = "groupID"][[2]]]
        gInfo[, ims_parent_group := NA_character_]
        if (!is.null(fTable[["CCS"]]))
            gInfo[, CCS := fTable[, mean(CCS), by = "groupID"][[2]]]
        setorderv(gInfo, c("mz", "ret", "mobility"), na.last = FALSE)
        gInfo[, group := fifelse(!is.na(mobility), makeIMSFGroupName(.I, ret, mz, mobility), makeFGroupName(.I, ret, mz))]
    }
    else
    {
        setorderv(gInfo, c("mz", "ret"))
        gInfo[, group := makeFGroupName(.I, ret, mz)]
    }
    setcolorder(gInfo, "group")
    
    gTable <- data.table(matrix(0, nrow = nrow(anaInfo), ncol = nrow(gInfo)))
    setnames(gTable, gInfo$group)
    
    ftindex <- data.table(matrix(0L, nrow = nrow(anaInfo), ncol = nrow(gInfo)))
    setnames(ftindex, gInfo$group)

    fTable[, group := gInfo$group[match(groupID, gInfo$groupID)]]
    fTable[, anaInfoInd := match(analysis, anaInfo$analysis)]
    for (row in seq_len(nrow(fTable)))
    {
        set(gTable, fTable$anaInfoInd[row], fTable$group[row], fTable$intensity[row])
        set(ftindex, fTable$anaInfoInd[row], fTable$group[row], fTable$anaRow[row])
    }
        
    gInfo[, groupID := NULL]
    
    if (hasMob && any(!is.na(fTable$ims_parent_ID)))
    {
        fTableNoIMSPar <- lapply(featureTable(feat), function(ft)
        {
            ft <- copy(ft)
            ft[, ims_parent_ID := NA_character_]
            return(ft)
        })
        featureTable(feat) <- fTableNoIMSPar
        warning("Any links between IMS parents and mobility features are removed!", call. = FALSE)
    }
    
    ret <- featureGroupsGreedy(groups = gTable, groupInfo = gInfo, features = feat, ftindex = ftindex)
    saveCacheData("groupFeaturesGreedy", ret, hash)
    return(ret)
}
