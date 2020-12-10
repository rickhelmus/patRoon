#' @include features.R
#' @include feature_groups.R
NULL

#' @rdname featureGroups-class
#' @export
featureGroupsOpenMS <- setClass("featureGroupsOpenMS", contains = "featureGroups")

setMethod("initialize", "featureGroupsOpenMS",
          function(.Object, ...) callNextMethod(.Object, algorithm = "openms", ...))


#' @details \code{groupFeaturesOpenMS} uses the OpenMS tools for grouping of
#'   features (see \url{http://www.openms.de}). Retention times may be aligned
#'   by the
#'   \href{https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_MapAlignerPoseClustering.html}{MapAlignerPoseClustering}
#'    TOPP tool. Grouping is achieved by either the
#'   \href{https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_FeatureLinkerUnlabeled.html}{FeatureLinkerUnlabeled}
#'    or
#'   \href{https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_FeatureLinkerUnlabeledQT.html}{FeatureLinkerUnlabeledQT}
#'    TOPP tools.
#'
#' @param QT If enabled, use \command{FeatureLinkerUnlabeledQT} instead of
#'   \command{FeatureLinkerUnlabeled} for feature grouping.
#' @param maxAlignRT,maxAlignMZ Used for retention alignment. Maximum retention
#'   time or m/z difference (seconds/Dalton) for feature pairing. Sets
#'   \code{-algorithm:pairfinder:distance_RT:max_difference} and
#'   \code{-algorithm:pairfinder:distance_MZ:max_difference} otpions,
#'   respectively.
#' @param maxGroupRT,maxGroupMZ as \code{maxAlignRT} and \code{maxAlignMZ}, but
#'   for grouping of features. Sets \code{-algorithm:distance_RT:max_difference}
#'   and \code{-algorithm:distance_MZ:max_difference} options, respectively.
#' @param extraOptsRT,extraOptsGroup Named \code{list} containing extra options
#'   that will be passed to \command{MapAlignerPoseClustering} or
#'   \command{FeatureLinkerUnlabeledQT/FeatureLinkerUnlabeled}, respectively.
#'   Any options specified here will override any of the above. Example:
#'   \code{extraOptsGroup=list("-algorithm:distance_RT:max_difference"=12)}
#'   (corresponds to setting \code{maxGroupRT=12}). Set to \code{NULL} to
#'   ignore.
#'
#' @template refs-openms
#'
#' @rdname feature-grouping
#' @export
setMethod("groupFeaturesOpenMS", "features", function(feat, rtalign = TRUE, QT = FALSE, maxAlignRT = 30,
                                                      maxAlignMZ = 0.005, maxGroupRT = 12,
                                                      maxGroupMZ = 0.005, extraOptsRT = NULL,
                                                      extraOptsGroup = NULL, verbose = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(feat, "features", add = ac)
    aapply(checkmate::assertFlag, . ~ rtalign + QT + verbose, fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ maxAlignRT + maxAlignMZ + maxGroupRT + maxGroupMZ,
           finite = TRUE, lower = 0, fixed = list(add = ac))
    aapply(checkmate::assertList, . ~ extraOptsRT + extraOptsGroup, any.missing = FALSE,
           names = "unique", null.ok = TRUE, fixed = list(add = ac))
    checkmate::reportAssertions(ac)

    if (length(feat) == 0)
        return(featureGroupsOpenMS(analysisInfo = analysisInfo(feat), features = feat))

    hash <- makeHash(feat, rtalign, QT, maxAlignRT, maxAlignMZ, maxGroupRT, maxGroupMZ, extraOptsRT, extraOptsGroup)
    cachefg <- loadCacheData("featureGroupsOpenMS", hash)
    if (!is.null(cachefg))
        return(cachefg)

    if (verbose)
        cat("Grouping features with OpenMS...\n===========\n")

    cfile <- tempfile("cons", fileext = ".consensusXML")
    generateConsensusXML(feat, cfile, rtalign, QT, maxAlignRT, maxAlignMZ, maxGroupRT, maxGroupMZ,
                         extraOptsRT, extraOptsGroup, verbose)
    fgimp <- importConsensusXML(feat, cfile, verbose)

    ret <- featureGroupsOpenMS(groups = fgimp$groups, groupInfo = fgimp$gInfo, analysisInfo = analysisInfo(feat),
                               features = feat, ftindex = fgimp$ftindex)

    saveCacheData("featureGroupsOpenMS", ret, hash)

    if (verbose)
        cat("\n===========\nDone!\n")

    return(ret)
})

generateConsensusXML <- function(feat, out, rtalign, QT, maxAlignRT, maxAlignMZ, maxGroupRT,
                                 maxGroupMZ, extraOptsRT, extraOptsGroup, verbose)
{
    sGroup <- analysisInfo(feat)
    fts <- featureTable(feat)

    edtaFiles <- list()
    featFiles <- list()
    for (datafile in sGroup$analysis)
    {
        featFiles[[datafile]] <- tempfile(datafile, fileext = ".featureXML")
        writeFeatureXML(fts[[datafile]], featFiles[[datafile]])
    }

    if (rtalign)
    {
        settings <- list("-algorithm:max_num_peaks_considered" = -1,
                         "-algorithm:superimposer:mz_pair_max_distance" = maxAlignMZ,
                         "-algorithm:superimposer:num_used_points" = 10000,
                         "-algorithm:pairfinder:distance_RT:max_difference" = maxAlignRT,
                         "-algorithm:pairfinder:distance_MZ:max_difference" = maxAlignMZ,
                         "-algorithm:pairfinder:distance_MZ:unit" = "Da")
        if (!is.null(extraOptsRT))
            settings <- modifyList(settings, extraOptsRT)
        
        executeCommand(getCommandWithOptPath("MapAlignerPoseClustering", "OpenMS"),
                       c(OpenMSArgListToOpts(settings), "-in", featFiles, "-out", featFiles),
                       stdout = if (verbose) "" else FALSE,
                       stderr = if (verbose) "" else FALSE)
    }

    settings <- list("-algorithm:distance_RT:max_difference" = maxGroupRT,
                     "-algorithm:distance_MZ:max_difference" = maxGroupMZ,
                     "-algorithm:distance_MZ:unit" = "Da")
    if (!is.null(extraOptsGroup))
        settings <- modifyList(settings, extraOptsGroup)
    
    executeCommand(getCommandWithOptPath(if (QT) "FeatureLinkerUnlabeledQT" else "FeatureLinkerUnlabeled", "OpenMS"),
                   c(OpenMSArgListToOpts(settings), "-in", featFiles, "-out", out),
                   stdout = if (verbose) "" else FALSE,
                   stderr = if (verbose) "" else FALSE)
}

importConsensusXML <- function(feat, cfile, verbose)
{
    if (verbose)
        cat("Importing consensus XML...")

    fTable <- featureTable(feat)
    anaCount <- nrow(analysisInfo(feat))

    consXML <- parseFeatConsXMLFile(cfile, anaCount)

    gCount <- nrow(consXML$gInfo)
    if (gCount > 0)
    {
        # generate group intensity table
        gtab <- data.table(matrix(0, nrow = anaCount, ncol = gCount))
        ftindex <- as.data.table(consXML$ftindex)

        for (gi in seq_len(gCount))
        {
            for (ai in seq_len(anaCount))
            {
                fti <- ftindex[[gi]][ai]
                if (fti != 0)
                    set(gtab, ai, gi, fTable[[ai]][["intensity"]][fti])
            }
        }

        gNames <- sapply(seq_len(gCount), function(grpi) makeFGroupName(grpi, consXML$gInfo$rts[grpi],
                                                                        consXML$gInfo$mzs[grpi]))
        rownames(consXML$gInfo) <- gNames
        setnames(gtab, gNames)
        setnames(ftindex, gNames)

        ret <- list(groups = gtab, gInfo = consXML$gInfo, ftindex = ftindex)
    }
    else
        ret <- list(groups = data.table(), gInfo = data.frame(), ftindex = data.table())

    if (verbose)
        cat("Done!\n")

    return(ret)
}
