# UNDONE: this should probably be deleted someday...

# nocov start

#' @include features.R
#' @include feature_groups.R
NULL

#' @rdname featureGroups-class
#' @export
featureGroupsEnviMass <- setClass("featureGroupsEnviMass", contains = "featureGroups")

setMethod("initialize", "featureGroupsEnviMass",
          function(.Object, ...) callNextMethod(.Object, algorithm = "envimass", ...))


#' Imports feature groups from enviMass
#'
#' Imports a 'profiles' produced by \pkg{enviMass}.
#'
#' @templateVar algo enviMass
#' @templateVar generic importFeatureGroups
#' @templateVar algoParam envimass
#' @template algo_importer
#'
#' @details This function \emph{only} imports 'raw' profiles, \emph{not} any results from further componentization steps
#'   performed in \pkg{enviMass}. Furthermore, this functionality has only been tested with older versions of
#'   \pkg{enviMass}. Finally, please note that this function only supports features imported by
#'   \code{\link{importFeaturesEnviMass}} (obviously, the same project should be used for both importing functions).
#'
#' @param path The path of the enviMass project.
#' @param feat The \code{\link{features}} object obtained with \code{\link{importFeaturesEnviMass}}.
#' @param positive Whether data from positive (\code{TRUE}) or negative (\code{FALSE}) should be loaded.
#'
#' @inherit importFeatureGroups return
#' 
#' @export
importFeatureGroupsEnviMass <- function(path, feat, positive)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertDirectoryExists(path, "r", add = ac)
    checkmate::assertClass(feat, "featuresEnviPick", add = ac)
    checkmate::assertFlag(positive, add = ac)
    
    if (ac$isEmpty())
    {
        resPath <- file.path(path, "results")
        profPeaksPath <- file.path(resPath, sprintf("profpeaks_%s", if (positive) "pos" else "neg"))
        profListPath <- file.path(resPath, sprintf("profileList_%s", if (positive) "pos" else "neg"))
        checkmate::assertFileExists(profPeaksPath, "r", .var.name = "path", add = ac)
        checkmate::assertFileExists(profListPath, "r", .var.name = "path", add = ac)
    }
        
    checkmate::reportAssertions(ac)
    
    hash <- makeHash(feat, path, positive)
    cachefg <- loadCacheData("featureGroupsEnviMass", hash)
    if (!is.null(cachefg))
        return(cachefg)

    cat("Importing enviMass feature groups (profiles)...\n")

    anaInfo <- analysisInfo(feat)
    fTable <- featureTable(feat)

    profPeaks <- as.data.table(loadRData(profPeaksPath))
    setorder(profPeaks, profileID)
    gInfo <- data.frame(mzs = profPeaks$`mean_m/z`, rts = profPeaks$mean_RT)
    rownames(gInfo) <- makeFGroupName(seq_len(nrow(gInfo)), gInfo$rts, gInfo$mzs)

    profList <- loadRData(profListPath)
    peaks <- as.data.table(profList$peaks)

    groups <- as.data.table(matrix(0, nrow(anaInfo), nrow(gInfo)))
    setnames(groups, rownames(gInfo))
    ftind <- copy(groups)

    gcount <- ncol(groups)
    prog <- openProgBar(0, gcount)

    for (grpi in seq_along(groups))
    {
        p <- peaks[profileIDs == grpi]
        sids <- match(as.character(p$sampleIDs), anaInfo$analysis)

        set(groups, sids, grpi, p$intensity)

        for (si in sids)
            set(ftind, si, grpi, match(p[sampleIDs == as.numeric(anaInfo$analysis[si]), peakIDs],
                                       fTable[[anaInfo$analysis[si]]]$ID))

        setTxtProgressBar(prog, grpi)
    }

    setTxtProgressBar(prog, gcount)
    close(prog)

    ret <- featureGroupsEnviMass(groups = groups, groupInfo = gInfo, features = feat, ftindex = ftind)

    saveCacheData("featureGroupsEnviMass", ret, hash)

    cat("Done!\n")

    return(ret)
}

# nocov end
