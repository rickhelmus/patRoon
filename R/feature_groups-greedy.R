#' @include main.R
NULL

#' @rdname featureGroups-class
#' @export
featureGroupsGreedy <- setClass("featureGroupsGreedy", contains = "featureGroups")

setMethod("initialize", "featureGroupsGreedy",
          function(.Object, ...) callNextMethod(.Object, algorithm = "greedy", ...))

calcFeatDist <- function(ret, mz, mob, retRef, mzRef, mobRef, rtWindow, mzWindow, IMSWindow)
{
    retd <- abs(ret - retRef); mzd <- abs(mz - mzRef)
    if (is.null(mob))
        return(sqrt((retd / rtWindow)^2 + (mzd / mzWindow)^2))
    mobd <- abs(mob - mobRef)
    return(sqrt((retd / rtWindow)^2 + (mzd / mzWindow)^2 + (mobd / IMSWindow)^2))
}

calcGroupScore <- function(fTableGrp, curAna, rtWindow, mzWindow, IMSWindow, withMob)
{
    # add eps: avoid zeros
    getDim <- \(x, win) ((max(x) - min(x)) / win) + .Machine$double.eps
    
    wh <- fTableGrp[analysis == curAna, which = TRUE]
    
    refMobs <- if (withMob) fTableGrp$mobility[wh] else rep(NA_real_, length(wh))
    ret <- mapply(fTableGrp$ret[wh], fTableGrp$mz[wh], refMobs, FUN = function(refRT, refMZ, refMob)
    {
        # select potential group
        fTableGrpOtherAna <- fTableGrp[analysis != curAna]
        
        # check if this feature would lead to an invalid group
        grpRetD <- max(fTableGrpOtherAna$ret, refRT) - min(fTableGrpOtherAna$ret, refRT)
        grpMZD <- max(fTableGrpOtherAna$mz, refMZ) - min(fTableGrpOtherAna$mz, refMZ)
        if (grpRetD > (rtWindow * 2) || grpMZD > (mzWindow * 2))
            return(-1)
        if (withMob)
        {
            grpMobD <- max(fTableGrpOtherAna$mobility, refMob) - min(fTableGrpOtherAna$mobility, refMob)
            if (grpMobD > (IMSWindow * 2))
                return(-1)
        }
        
        fTableGrpOtherAna[, keep := {
            if (.N == 1L)
                TRUE
            else
            {
                dists <- calcFeatDist(ret, mz, if (withMob) mobility else NULL, refRT, refMZ, refMob, rtWindow, mzWindow, IMSWindow)
                seq_len(.N) == which.min(dists)
            }
        }, by = "analysis"]
        fTableGrpOtherAna <- fTableGrpOtherAna[keep == TRUE]
        
        # calc score
        
        if (nrow(fTableGrpOtherAna) == 0)
            return(0)
        
        volume <- getDim(c(refRT, fTableGrpOtherAna$ret), rtWindow) *
            getDim(c(refMZ, fTableGrpOtherAna$mz), mzWindow)
        if (withMob)
            volume <- volume * getDim(c(refMob, fTableGrpOtherAna$mobility), IMSWindow)
        
        # higher score for less volume and more members
        return(1/volume + nrow(fTableGrpOtherAna) * 0.1)
    })
}

#' @export
groupFeaturesGreedy <- function(features, rtWindow = defaultLim("retention", "medium"),
                                mzWindow = defaultLim("mz", "medium"),
                                IMSWindow = defaultLim("mobility", "medium"),
                                verbose = TRUE, useCPP)
{
    # UNDONE: rtalign arg?
    
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
    
    doGroupCPP <- function(withMob)
    {
        mobs <- if (!withMob)
            rep(-100, nrow(fTable)) # UNDONE!!
        else
            fTable$mobility
        fTable[, groupID := getGroupIDs(ret, mz, mobs, intensity, match(analysis, anaInfo$analysis), rtWindow,
                                        mzWindow, IMSWindow)]
        curGroup <<- max(fTable$groupID) + 1L
    }
    
    if (hasMob)
    {
        if (verbose)
            message("Grouping features with mobilities...")
        if (useCPP)
            doGroupCPP(TRUE)
        else
            doGroup(TRUE)
    }
    if (!hasMob || any(is.na(fTable$mobility)))
    {
        if (verbose && hasMob)
            message("Grouping features without mobilities...")
        if (useCPP)
            doGroupCPP(FALSE)
        else
            doGroup(FALSE)
    }
    
    # some duplicate features in a single analysis may be unassigned --> fix here
    fTable[is.na(groupID) | groupID < 0, groupID := seq_len(.N) + curGroup]
    
    gInfo <- fTable[, .(ret = mean(ret), mz = mean(mz)), by = "groupID"]
    if (hasMob)
    {
        gInfo[, mobility := fTable[, mean(mobility), by = "groupID"][[2]]]
        gInfo[, ims_parent_group := NA_character_]
        if (!is.null(fTable[["CCS"]]))
            gInfo[, CCS := fTable[, mean(CCS), by = "groupID"][[2]]]
        setorderv(gInfo, c("ret", "mz", "mobility"), na.last = FALSE)
        gInfo[, group := fifelse(!is.na(mobility), makeIMSFGroupName(.I, ret, mz, mobility), makeFGroupName(.I, ret, mz))]
    }
    else
    {
        setorderv(gInfo, c("ret", "mz"))
        gInfo[, group := makeFGroupName(.I, ret, mz)]
    }
    setcolorder(gInfo, "group")
    
    gTable <- data.table(matrix(0, nrow = nrow(anaInfo), ncol = nrow(gInfo)))
    setnames(gTable, gInfo$group)
    
    ftindex <- data.table(matrix(0L, nrow = nrow(anaInfo), ncol = nrow(gInfo)))
    setnames(ftindex, gInfo$group)
    
    for (grpi in seq_len(nrow(gInfo)))
    {
        ft <- fTable[groupID == gInfo$groupID[grpi]]
        anai <- match(ft$analysis, anaInfo$analysis) # align analysis order
        set(gTable, anai, gInfo$group[grpi], ft$intensity)
        set(ftindex, anai, gInfo$group[grpi], ft$anaRow)
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
