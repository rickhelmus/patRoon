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

calcGroupScore <- function(fTableGrp, curAna, rtWindow, mzWindow, IMSWindow)
{
    # add eps: avoid zeros
    getDim <- \(x, win) ((max(x) - min(x)) / win) + .Machine$double.eps
    
    wh <- fTableGrp[analysis == curAna, which = TRUE]
    
    ret <- mapply(fTableGrp$ret[wh], fTableGrp$mz[wh], FUN = function(refRT, refMZ)
    {
        # select potential group
        # UNDONE: also keep track if max/min of group stays within tolerance?
        
        fTableGrpOtherAna <- fTableGrp[analysis != curAna]
        
        grpRetD <- max(fTableGrpOtherAna$ret, refRT) - min(fTableGrpOtherAna$ret, refRT)
        grpMZD <- max(fTableGrpOtherAna$mz, refMZ) - min(fTableGrpOtherAna$mz, refMZ)
        if (grpRetD > rtWindow || grpMZD > mzWindow)
            return(-1) # this feature would lead to an invalid group
        
        fTableGrpOtherAna[, keep := {
            if (.N == 1L)
                TRUE
            else
            {
                dists <- calcFeatDist(ret, mz, NULL, refRT, refMZ, NULL, rtWindow, mzWindow, IMSWindow)
                seq_len(.N) == which.min(dists)
            }
        }, by = "analysis"]
        fTableGrpOtherAna <- fTableGrpOtherAna[keep == TRUE]
        
        # calc score
        
        if (nrow(fTableGrpOtherAna) == 0)
            return(0)
        
        volume <- getDim(c(refRT, fTableGrpOtherAna$ret), rtWindow) *
            getDim(c(refMZ, fTableGrpOtherAna$mz), mzWindow)
        
        # higher score for less volume and more members
        return(1/volume + nrow(fTableGrpOtherAna) * 0.1)
    })
}

#' @export
groupFeaturesGreedy <- function(features, rtWindow = defaultLim("retention", "medium"),
                                mzWindow = defaultLim("mz", "medium"),
                                IMSWindow = defaultLim("mobility", "medium"),
                                verbose = TRUE)
{
    # UNDONE: rtalign arg?
    # UNDONE: IMS support
    
    anaInfo <- analysisInfo(features)
    sqeps <- sqrt(.Machine$double.eps) # cache for numLTE below
    
    fTable <- as.data.table(features)
    fTable[, groupID := NA_integer_]
    fTable[, anaRow := seq_len(.N), by = "analysis"]
    
    curGroup <- 1L
    for (ftRow in order(fTable$intensity, decreasing = TRUE))
    {
        if (!is.na(fTable$groupID[ftRow]))
            next # already assigned
        
        grpInds <- fTable[is.na(groupID) & numLTE(abs(ret - ret[ftRow]), rtWindow, sqeps) & numLTE(abs(mz - mz[ftRow]), mzWindow, sqeps),
                          which = TRUE]
        if (length(grpInds) > 1 && anyDuplicated(fTable$analysis[grpInds]))
        {
            fTable[grpInds, groupID := {
                scores <- calcGroupScore(fTable[grpInds], analysis, rtWindow, mzWindow, IMSWindow)
                gid <- rep(NA_integer_, .N)
                gid[seq_along(gid) == which.max(scores) & scores >= 0] <- curGroup # select highest score if duplicate feature in one analysis
                gid
            }, by = "analysis"]
        }
        else
            fTable[grpInds, groupID := curGroup]
        
        curGroup <- curGroup + 1L
    }
    
    # some duplicate features in a single analysis may be unassigned --> fix here
    fTable[is.na(groupID), groupID := seq_len(.N) + curGroup]
    
    gInfo <- fTable[, .(ret = mean(ret), mz = mean(mz)), by = "groupID"]
    gInfo[, group := makeFGroupName(seq_len(nrow(gInfo)), ret, mz)]
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
    
    return(featureGroupsGreedy(groups = gTable, groupInfo = gInfo, features = features, ftindex = ftindex))
}
