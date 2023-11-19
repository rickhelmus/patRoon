#' @include features.R
#' @include feature_groups.R
NULL

#' @rdname featureGroups-class
#' @export
featureGroupsOpenMS <- setClass("featureGroupsOpenMS", contains = "featureGroups")

setMethod("initialize", "featureGroupsOpenMS",
          function(.Object, ...) callNextMethod(.Object, algorithm = "openms", ...))


#' Group features using OpenMS
#'
#' Group and align features with OpenMS tools
#'
#' @templateVar algo OpenMS
#' @templateVar do group features
#' @templateVar generic groupFeatures
#' @templateVar algoParam openms
#' @template algo_generator
#'
#' @details Retention times may be aligned by the
#'   \href{https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_MapAlignerPoseClustering.html}{MapAlignerPoseClustering}
#'    TOPP tool. Grouping is achieved by either the
#'   \href{https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_FeatureLinkerUnlabeled.html}{FeatureLinkerUnlabeled}
#'    or
#'   \href{https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_FeatureLinkerUnlabeledQT.html}{FeatureLinkerUnlabeledQT}
#'    TOPP tools.
#'
#' @template feat-arg
#' @template rtalign-arg
#'
#' @param QT If enabled, use \command{FeatureLinkerUnlabeledQT} instead of \command{FeatureLinkerUnlabeled} for feature
#'   grouping.
#' @param maxAlignRT,maxAlignMZ Used for retention alignment. Maximum retention time or m/z difference (seconds/Dalton)
#'   for feature pairing. Sets \code{-algorithm:pairfinder:distance_RT:max_difference} and
#'   \code{-algorithm:pairfinder:distance_MZ:max_difference} otpions, respectively.
#' @param maxGroupRT,maxGroupMZ as \code{maxAlignRT} and \code{maxAlignMZ}, but for grouping of features. Sets
#'   \code{-algorithm:distance_RT:max_difference} and \code{-algorithm:distance_MZ:max_difference} options,
#'   respectively.
#' @param extraOptsRT,extraOptsGroup Named \code{list} containing extra options that will be passed to
#'   \command{MapAlignerPoseClustering} or \command{FeatureLinkerUnlabeledQT/FeatureLinkerUnlabeled}, respectively. Any
#'   options specified here will override any of the above. Example:
#'   \code{extraOptsGroup=list("-algorithm:distance_RT:max_difference"=12)} (corresponds to setting
#'   \code{maxGroupRT=12}). Set to \code{NULL} to ignore.
#'
#' @inheritParams groupFeatures
#'
#' @inherit groupFeatures return
#'
#' @template refs-openms
#'
#' @templateVar what groupFeaturesOpenMS
#' @templateVar cl features
#' @template main-rd-method
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
        return(featureGroupsOpenMS(features = feat))

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

    ret <- featureGroupsOpenMS(groups = fgimp$groups, groupInfo = fgimp$gInfo, features = feat, ftindex = fgimp$ftindex)

    saveCacheData("featureGroupsOpenMS", ret, hash)

    if (verbose)
        cat("\n===========\nDone!\n")

    return(ret)
})

generateConsensusXML <- function(feat, out, rtalign, QT, maxAlignRT, maxAlignMZ, maxGroupRT,
                                 maxGroupMZ, extraOptsRT, extraOptsGroup, verbose)
{
    anas <- analysisInfo(feat)$analysis
    fts <- featureTable(feat)

    getIni <- function(cmd, inFiles, outFiles)
    {
        iniFile <- tempfile(fileext = ".ini")
        executeCommand(cmd, c("-write_ini", iniFile), stdout = FALSE, stderr = FALSE)
        addFilesToOpenMSIni(iniFile, inFiles, outFiles)
        return(iniFile)
    }
    
    featFiles <- character()
    for (datafile in anas)
    {
        featFiles[datafile] <- tempfile(datafile, fileext = ".featureXML")
        writeFeatureXML(fts[[datafile]], datafile, featFiles[datafile], FALSE)
    }
    
    # NOTE: we use ini files to specify the input/output files below to avoid troubles on Windows with large number of
    # analyses, see https://github.com/OpenMS/OpenMS/issues/6845

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
        
        ftNonEmpty <- sapply(fts, nrow) > 0
        
        cmd <- getExtDepPath("openms", "MapAlignerPoseClustering")
        iniFile <- getIni(cmd, featFiles[ftNonEmpty], featFiles[ftNonEmpty])
        executeCommand(cmd, c(OpenMSArgListToOpts(settings), "-ini", iniFile),
                       stdout = if (verbose) "" else FALSE,
                       stderr = if (verbose) "" else FALSE)
    }

    settings <- list("-algorithm:distance_RT:max_difference" = maxGroupRT,
                     "-algorithm:distance_MZ:max_difference" = maxGroupMZ,
                     "-algorithm:distance_MZ:unit" = "Da")
    if (!is.null(extraOptsGroup))
        settings <- modifyList(settings, extraOptsGroup)
    
    cmd <- getExtDepPath("openms", if (QT) "FeatureLinkerUnlabeledQT" else "FeatureLinkerUnlabeled")
    iniFile <- getIni(cmd, featFiles, character())
    executeCommand(cmd, c(OpenMSArgListToOpts(settings), "-ini", iniFile, "-out", out),
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
