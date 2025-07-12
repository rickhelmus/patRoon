#' @include main.R
NULL

#' @rdname featureGroups-class
#' @export
featureGroupsGreedy <- setClass("featureGroupsGreedy", contains = "featureGroups")

setMethod("initialize", "featureGroupsGreedy",
          function(.Object, ...) callNextMethod(.Object, algorithm = "greedy", ...))

#' @export
setMethod("groupFeaturesGreedy", "features", function(feat, rtalign = FALSE,
                                                      rtWindow = defaultLim("retention", "medium"),
                                                      mzWindow = defaultLim("mz", "medium"),
                                                      IMSWindow = defaultLim("mobility", "medium"),
                                                      scoreWeights = c(retention = 1, mz = 1, mobility = 1),
                                                      verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertFlag(rtalign, add = ac)
    aapply(checkmate::assertNumber, . ~ rtWindow + mzWindow + IMSWindow, lower = 0, finite = TRUE,
           fixed = list(add = ac))
    checkmate::assertNumeric(scoreWeights, len = 3, lower = 0, finite = TRUE, any.missing = FALSE, add = ac)
    assertHasNames(scoreWeights, c("retention", "mz", "mobility"), add = ac)
    checkmate::assertFlag(verbose, add = ac)
    checkmate::reportAssertions(ac)
    
    if (rtalign)
        stop("Retention time alignment (rtalign=TRUE) is not yet supported for greedy grouping!", call. = FALSE)
    
    hash <- makeHash(feat, rtWindow, mzWindow, IMSWindow)
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
    
    fTable[, groupID := getGroupIDs(ret, mz, mobility, intensity, match(analysis, anaInfo$analysis), rtWindow,
                                    mzWindow, IMSWindow, scoreWeights)]
    # some duplicate features in a single analysis may be unassigned --> fix here
    fTable[groupID < 0, groupID := seq_len(.N) + max(fTable$groupID)]

    if (verbose)
        printf("Done! Found %d groups.\n", max(fTable$groupID))
    
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
})
