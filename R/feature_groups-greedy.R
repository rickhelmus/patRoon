#' @include main.R
NULL

#' @rdname featureGroups-class
#' @export
featureGroupsGreedy <- setClass("featureGroupsGreedy", contains = "featureGroups")

setMethod("initialize", "featureGroupsGreedy",
          function(.Object, ...) callNextMethod(.Object, algorithm = "greedy", ...))

calcFeatDist <- function(ret, mz, mob, retRef, mzRef, mobRef, rtWindow, mzWindow, IMSWindow)
{
    retd <- abs(ret - retRef); mzd <- abs(mz - mzRef); mobd <- abs(mob - mobRef)
    return(sqrt((retd / rtWindow)^2 + (mzd / mzWindow)^2 + (mobd / IMSWindow)^2))
}

calcGroupScore <- function(fTableGrp, curAna, rtWindow, mzWindow, IMSWindow)
{
    # add eps: avoid zeros
    getDim <- \(x, win) ((max(x) - min(x)) / win) + .Machine$double.eps
    
    wh <- fTableGrp[analysis == curAna, which = TRUE]
    
    ret <- mapply(fTableGrp$ret[wh], fTableGrp$mz[wh], fTableGrp$mobility[wh], FUN = function(refRT, refMZ, refMob)
    {
        # select potential group
        fTableGrpOtherAna <- fTableGrp[analysis != curAna]
        
        # check if this feature would lead to an invalid group
        grpRetD <- max(fTableGrpOtherAna$ret, refRT) - min(fTableGrpOtherAna$ret, refRT)
        grpMZD <- max(fTableGrpOtherAna$mz, refMZ) - min(fTableGrpOtherAna$mz, refMZ)
        grpMobD <- max(fTableGrpOtherAna$mobility, refMob) - min(fTableGrpOtherAna$mobility, refMob)
        if (grpRetD > (rtWindow * 2) || grpMZD > (mzWindow * 2) || grpMobD > (IMSWindow * 2))
            return(-1)

        fTableGrpOtherAna[, keep := {
            if (.N == 1L)
                TRUE
            else
            {
                dists <- calcFeatDist(ret, mz, mobility, refRT, refMZ, refMob, rtWindow, mzWindow, IMSWindow)
                seq_len(.N) == which.min(dists)
            }
        }, by = "analysis"]
        fTableGrpOtherAna <- fTableGrpOtherAna[keep == TRUE]
        
        # calc score
        
        if (nrow(fTableGrpOtherAna) == 0)
            return(0)
        
        volume <- getDim(c(refRT, fTableGrpOtherAna$ret), rtWindow) *
            getDim(c(refMZ, fTableGrpOtherAna$mz), mzWindow) *
            getDim(c(refMob, fTableGrpOtherAna$mobility), IMSWindow)

        # higher score for less volume and more members
        return(1/volume + nrow(fTableGrpOtherAna) * 0.1)
    })
}

#' @export
groupFeaturesGreedy <- function(features, rtWindow = defaultLim("retention", "medium"),
                                mzWindow = defaultLim("mz", "medium"),
                                IMSWindow = defaultLim("mobility", "medium"),
                                verbose = TRUE, useCPP = TRUE)
{
    # UNDONE: rtalign arg?
    # UNDONE: remove useCPP and R code
    
    # add checkmates
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertNumber, . ~ rtWindow + mzWindow + IMSWindow, lower = 0, finite = TRUE,
           fixed = list(add = ac))
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    hash <- makeHash(features, rtWindow, mzWindow, IMSWindow)
    cd <- loadCacheData("groupFeaturesGreedy", hash)
    if (!is.null(cd))
        return(cd)
    
    anaInfo <- analysisInfo(features)
    hasMob <- hasMobilities(features)
    sqeps <- sqrt(.Machine$double.eps) # cache for numLTE below
    
    fTable <- as.data.table(features)
    fTable[, groupID := NA_integer_]
    fTable[, anaRow := seq_len(.N), by = "analysis"]
    fTable[, row := seq_len(.N)]
    
    curGroup <- 0L
    
    # UNDONE: remove
    doGroup <- function(withMob)
    {
        for (ftRow in order(fTable$intensity, decreasing = TRUE))
        {
            if (!is.na(fTable$groupID[ftRow]))
                next # already assigned
            
            # NOTE: withMob is only TRUE if hasMob is
            if ((withMob && is.na(fTable$mobility[ftRow])) || (!withMob && hasMob && !is.na(fTable$mobility[ftRow])))
                next
            
            fTableGrp <- fTable[is.na(groupID) & numLTE(abs(ret - ret[ftRow]), rtWindow, sqeps) &
                                    numLTE(abs(mz - mz[ftRow]), mzWindow, sqeps)]
            if (withMob)
                fTableGrp <- fTableGrp[numLTE(abs(mobility - fTable$mobility[ftRow]), IMSWindow, sqeps)]
            
            if (nrow(fTableGrp) > 1 && anyDuplicated(fTableGrp$analysis))
            {
                fTable[match(fTableGrp$row, row), groupID := {
                    scores <- calcGroupScore(fTableGrp, analysis, rtWindow, mzWindow, IMSWindow, withMob)
                    gid <- rep(NA_integer_, .N)
                    gid[seq_along(gid) == which.max(scores) & scores >= 0] <- curGroup # select highest score if duplicate feature in one analysis
                    gid
                }, by = "analysis"]
            }
            else
                fTable[match(fTableGrp$row, row), groupID := curGroup]
            
            # if (curGroup == 117) browser()
            
            if (nrow(fTableGrp) > 0)
                curGroup <<- curGroup + 1L
        }
    }
    
    # HACK: add dummy mobilities so grouping code can handle features without mobilities
    if (!hasMob || all(is.na(fTable$mobility)))
        fTable[, mobility := -100] # some random value so mobilities will not influence grouping
    else if (any(is.na(fTable$mobility)))
    {
        unreasonableMob <- max(fTable$mobility, na.rm = TRUE) + (10 * IMSWindow)
        fTable[, mobility_orig := mobility] # we need the original data later
        fTable[is.na(mobility), mobility := unreasonableMob]
    }

    if (useCPP)
    {
        fTable[, groupID := getGroupIDs(ret, mz, mobility, intensity, match(analysis, anaInfo$analysis), rtWindow,
                                        mzWindow, IMSWindow, verbose)]
        curGroup <- max(fTable$groupID) + 1L
    }
    else
    {
        for (ftRow in order(fTable$intensity, decreasing = TRUE))
        {
            if (!is.na(fTable$groupID[ftRow]))
                next # already assigned
            
            fTableGrp <- fTable[is.na(groupID) & numLTE(abs(ret - ret[ftRow]), rtWindow, sqeps) &
                                    numLTE(abs(mz - mz[ftRow]), mzWindow, sqeps) &
                                    numLTE(abs(mobility - fTable$mobility[ftRow]), IMSWindow, sqeps)]
            
            if (nrow(fTableGrp) > 1 && anyDuplicated(fTableGrp$analysis))
            {
                fTable[match(fTableGrp$row, row), groupID := {
                    scores <- calcGroupScore(fTableGrp, analysis, rtWindow, mzWindow, IMSWindow)
                    gid <- rep(NA_integer_, .N)
                    gid[seq_along(gid) == which.max(scores) & scores >= 0] <- curGroup # select highest score if duplicate feature in one analysis
                    gid
                }, by = "analysis"]
            }
            else
                fTable[match(fTableGrp$row, row), groupID := curGroup]
            
            if (nrow(fTableGrp) > 0)
                curGroup <- curGroup + 1L
        }
    }
    
    # some duplicate features in a single analysis may be unassigned --> fix here
    fTable[is.na(groupID) | groupID < 0, groupID := seq_len(.N) + curGroup]
    
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
        fTableNoIMSPar <- lapply(featureTable(features), function(ft)
        {
            ft <- copy(ft)
            ft[, ims_parent_ID := NA_character_]
            return(ft)
        })
        featureTable(features) <- fTableNoIMSPar
        warning("Any links between IMS parents and mobility features are removed!", call. = FALSE)
    }
    
    ret <- featureGroupsGreedy(groups = gTable, groupInfo = gInfo, features = features, ftindex = ftindex)
    saveCacheData("groupFeaturesGreedy", ret, hash)
    return(ret)
}
