# SPDX-FileCopyrightText: 2016 - 2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
#' @include features-xcms.R
#' @include feature_groups.R
#' @include feature_groups-xcms.R
NULL

# internal functions copied from XCMS (https://github.com/sneumann/xcms)
# UNDONE: would be nice to not need this, e.g. by improving xcms interop
XCMSInternal <- setRefClass("XCMSInternal", methods = list(

    .PROCSTEP.PEAK.DETECTION = function() "Peak detection",
    .PROCSTEP.PEAK.GROUPING = function() "Peak grouping",
    .PROCSTEP.RTIME.CORRECTION = function() "Retention time correction",

    .featureIDs = function(x, prefix = "FT") {
        sprintf(paste0(prefix, "%0", ceiling(log10(x + 1L)), "d"), 1:x)
    },

    .update_feature_definitions = function(x, original_names, subset_names) {
        x$peakidx <- lapply(x$peakidx, function(z) {
            idx <- base::match(original_names[z], subset_names)
            idx[!is.na(idx)]
        })
        x[lengths(x$peakidx) > 0, ]
    }

))()

XCMSFeatCols <- function() c("mz", "mzmin", "mzmax", "ret", "retmin", "retmax", "intensity", "area", "sn", "ID")

isXCMSClass <- function(cl, name) is(cl, getClass(name, where = "xcms"))

readMSDataForXCMS3 <- function(anaInfo)
{
    anaFiles <- mapply(anaInfo$analysis, anaInfo$path, FUN = getMzMLOrMzXMLAnalysisPath, MoreArgs = list(mustExist = TRUE))
    return(MSnbase::readMSData(files = anaFiles,
                               pdata = Biobase::AnnotatedDataFrame(data.frame(sample_name = anaInfo$analysis,
                                                                              sample_group = anaInfo$group,
                                                                              stringsAsFactors = FALSE)),
                               mode = "onDisk"))
}

makeXCMSGroups <- function(fGroups, verbose = TRUE)
{
    fTable <- featureTable(fGroups)

    # add unique feat IDs (corresponding to ftindex)
    combfts <- lapply(copy(fTable), function(f) f[, ftID := seq_len(.N)])
    names(combfts) <- NULL # get rid of file names ...
    combfts <- rbindlist(combfts, idcol = "sind") # ... so that rbindlist will generate numeric IDs
    # combfts should now be aligned with peak table of xcms set

    combfts[, xsID := seq_len(.N)] # ID corresponding to XCMS peaks

    gTable <- groupTable(fGroups)
    ftind <- copy(groupFeatIndex(fGroups))
    anaInfo <- analysisInfo(fGroups)

    if (verbose)
        cat("Making group index table\n")

    setkeyv(combfts, c("sind", "ftID"))
    idx <- lapply(seq_along(ftind), function(g)
    {
        fti <- data.table(ftID = ftind[[g]], sind = seq_along(ftind[[g]]))
        setkeyv(fti, c("sind", "ftID"))
        # join by sample ID and by feat ID so that we end up with a dt having only
        # relevant rows for the IDs given in ftind
        combfts[fti, nomatch = 0][["xsID"]] # xsID is the XCMS peak ID
    })

    if (verbose)
        cat("Making groups...\n")

    repGroups <- replicateGroups(fGroups)
    grps <- rbindlist(lapply(seq_along(ftind), function(g)
    {
        fprops <- rbindlist(lapply(seq_len(nrow(anaInfo)), function(s)
        {
            fi <- ftind[[g]][s]
            if (fi == 0)
                data.table(mz = NA, ret = NA)
            else
                fTable[[anaInfo$analysis[s]]][fi, c("mz", "ret")]
        }))

        ret <- data.table(mzmed = median(fprops[["mz"]], na.rm = TRUE), mzmin = min(fprops[["mz"]], na.rm = TRUE),
                          mzmax = max(fprops[["mz"]], na.rm = TRUE), rtmed = median(fprops[["ret"]], na.rm = TRUE),
                          rtmin = min(fprops[["ret"]], na.rm = TRUE), rtmax = max(fprops[["ret"]], na.rm = TRUE),
                          npeaks = sum(ftind[[g]] != 0))

        # add counts per sample group
        ret[, (repGroups) := lapply(repGroups, function(rg) sum(ftind[[g]][rg == anaInfo$group] != 0))]

        return(ret)
    }))

    return(list(groups = grps, idx = idx))
}


#' @rdname xcms-conv
#' @export
setMethod("getXCMSSet", "features", function(obj, verbose, loadRawData)
{
    # generate dummy XCMS set, based on https://groups.google.com/forum/m/#!topic/xcms/CGC0SKMVhAQ

    checkmate::assertFlag(loadRawData)
    checkmate::assertFlag(verbose)

    xs <- new(getClassDef("xcmsSet", package = "xcms"))
    anaInfo <- analysisInfo(obj)
    
    if (loadRawData)
        verifyDataCentroided(anaInfo)
    
    xcms::phenoData(xs) <- data.frame(class = anaInfo$group, row.names = anaInfo$analysis)

    if (loadRawData)
        xcms::filepaths(xs) <- sapply(seq_len(nrow(anaInfo)),
                                      function(i) getMzMLOrMzXMLAnalysisPath(anaInfo$analysis[i], anaInfo$path[i],
                                                                             mustExist = TRUE),
                                      USE.NAMES = FALSE)
    else
        xcms::filepaths(xs) <- anaInfo$analysis # dummy paths

    fts <- featureTable(obj)
    plist <- list()
    rlist <- list(raw = vector("list", length(anaInfo$analysis)), corrected = vector("list", length(anaInfo$analysis)))

    for (i in seq_len(nrow(anaInfo)))
    {
        ft <- fts[[anaInfo$analysis[i]]]

        if (nrow(ft) > 0)
        {
            plist[[i]] <- data.frame(mz = ft$mz, mzmin = ft$mzmin, mzmax = ft$mzmax, rt = ft$ret,
                                     rtmin = ft$retmin, rtmax = ft$retmax, maxo = ft$intensity, into = ft$area,
                                     sample = i, sn = if (!is.null(ft[["sn"]])) ft$sn else NA_real_)
        }
        else
            plist[[i]] <- data.frame(mz = numeric(), mzmin = numeric(), mzmax = numeric(), rt = numeric(),
                                     rtmin = numeric(), rtmax = numeric(), maxo = numeric(), into = numeric(),
                                     sample = numeric(), sn = numeric())
        
        if (loadRawData)
        {
            xr <- loadXCMSRaw(anaInfo$analysis[i], anaInfo$path[i], verbose = verbose)[[1]]
            rlist$raw[[i]] <- xr@scantime
            rlist$corrected[[i]] <- xr@scantime
        }
    }

    xcms::peaks(xs) <- as.matrix(do.call(function(...) rbind(..., make.row.names = FALSE), plist))
    xs@rt <- rlist
    xcms::profinfo(xs) <- list(method = "bin", step = 0.1)

    return(xs)
})

#' @rdname xcms-conv
#' @export
setMethod("getXCMSSet", "featuresXCMS", function(obj, ...)
{
    return(obj@xs)
})

#' @rdname xcms-conv
#' @export
setMethod("getXCMSSet", "featureGroups", function(obj, verbose, loadRawData)
{
    checkmate::assertFlag(loadRawData)
    checkmate::assertFlag(verbose)

    if (verbose)
        cat("Getting ungrouped xcmsSet...\n")
    xs <- getXCMSSet(getFeatures(obj), verbose = verbose, loadRawData = loadRawData)

    xsgrps <- makeXCMSGroups(obj, verbose)
    xcms::groupidx(xs) <- xsgrps$idx
    xcms::groups(xs) <- as.matrix(xsgrps$groups)

    return(xs)
})

#' @rdname xcms-conv
#' @export
setMethod("getXCMSSet", "featureGroupsXCMS", function(obj, verbose, loadRawData)
{
    # first see if we can just return the xcmsSet used during grouping

    anaInfo <- analysisInfo(obj)

    if (length(xcms::filepaths(obj@xs)) != length(anaInfo$analysis) ||
        !all(simplifyAnalysisNames(xcms::filepaths(obj@xs)) == anaInfo$analysis))
    {
        # files changed, need to update group statistics which is rather complex so just fallback
        return(callNextMethod(obj, verbose = verbose, loadRawData = loadRawData))
    }

    return(obj@xs)
})

#' @rdname xcms-conv
#' @export
setMethod("getXCMSSet", "featuresSet", function(obj, ..., set) getXCMSSet(unset(obj, set), ...))

#' @rdname xcms-conv
#' @export
setMethod("getXCMSSet", "featureGroupsSet", function(obj, ..., set) getXCMSSet(unset(obj, set), ...))

#' @rdname xcms-conv
#' @export
setMethod("getXCMSnExp", "features", function(obj, verbose, loadRawData)
{
    checkmate::assertFlag(loadRawData)
    
    rawData <- NULL
    if (loadRawData)
    {
        verifyDataCentroided(analysisInfo(obj))
        rawData <- readMSDataForXCMS3(analysisInfo(obj))
    }
    else
    {
        # create a dummy MSnExp object
        
        # dummy spectra: 1 peak/analysis
        # UNDONE: more columns needed?
        dummySpecs <- data.frame(fileIdx = integer(), spIdx = integer(), seqNum = integer(), msLevel = integer(),
                                 acquisitionNum = integer(), retentionTime = numeric(), precursorScanNum = integer())
        rownames(dummySpecs) <- character()
        
        anaInfo <- analysisInfo(obj)
        rawData <- new("OnDiskMSnExp", processingData = new("MSnProcess",
                                                            files = anaInfo$analysis),
                       featureData = Biobase::AnnotatedDataFrame(data = dummySpecs))
        xcms::phenoData(rawData) <- Biobase::AnnotatedDataFrame(data.frame(sample_name = anaInfo$analysis,
                                                                           sample_group = anaInfo$group,
                                                                           stringsAsFactors = FALSE))
    }

    msLevel = 1L # UNDONE?

    # this is mainly based on xcms:::.peaks_to_result()
    xph <- new("XProcessHistory", param = NULL, date = date(),
               type = XCMSInternal$.PROCSTEP.PEAK.DETECTION(),
               fileIndex = seq_along(analyses(obj)),
               msLevel = msLevel)
    ret <- as(rawData, "XCMSnExp")
    ret@.processHistory <- c(xcms::processHistory(ret), list(xph))

    if (length(obj) > 0)
    {
        fTable <- featureTable(obj)

        # use unname: get numeric id column
        allFeats <- rbindlist(lapply(unname(fTable), function(ft)
        {
            if (nrow(ft) > 0)
            {
                ret <- data.table(mz = ft$mz, mzmin = ft$mzmin, mzmax = ft$mzmax, rt = ft$ret,
                                  rtmin = ft$retmin, rtmax = ft$retmax, maxo = ft$intensity, into = ft$area,
                                  sn = if (!is.null(ft[["sn"]])) ft$sn else NA_real_)
            }
            else
                ret <- data.table(mz = numeric(), mzmin = numeric(), mzmax = numeric(), rt = numeric(),
                                  rtmin = numeric(), rtmax = numeric(), maxo = numeric(), into = numeric(),
                                  sn = numeric())
            return(ret)
        }), idcol = "sample")
        setcolorder(allFeats, setdiff(names(allFeats), "sample")) # move sample col to end

        allFeats <- as.matrix(allFeats)
        rownames(allFeats) <- XCMSInternal$.featureIDs(nrow(allFeats), "CP")

        xcms::chromPeaks(ret) <- allFeats
        xcms::chromPeakData(ret)$ms_level <- msLevel
        xcms::chromPeakData(ret)$is_filled <- FALSE
    }

    return(ret)
})

#' @rdname xcms-conv
#' @export
setMethod("getXCMSnExp", "featuresXCMS3", function(obj, ...)
{
    return(obj@xdata)
})

#' @rdname xcms-conv
#' @export
setMethod("getXCMSnExp", "featureGroups", function(obj, verbose, loadRawData)
{
    checkmate::assertFlag(verbose)
    
    if (verbose)
        cat("Getting ungrouped XCMSnExp...\n")

    msLevel = 1L # UNDONE?

    xdata <- getXCMSnExp(getFeatures(obj), verbose = verbose, loadRawData = loadRawData)

    xsgrps <- makeXCMSGroups(obj, verbose)
    xsgrps$groups[, peakidx := list(xsgrps$idx)]

    # this is mainly based on xcms::groupChromPeaks()

    grps <- as(xsgrps$groups, "DataFrame")
    if (!all(xcms::chromPeakData(xdata)$ms_level %in% msLevel))
        df <- XCMSInternal$.update_feature_definitions(grps, rownames(xcms::chromPeaks(xdata, msLevel = msLevel)),
                                                 rownames(xcms::chromPeaks(xdata)))
    rownames(grps) <- XCMSInternal$.featureIDs(nrow(grps))
    xcms::featureDefinitions(xdata) <- grps

    xph <- new("XProcessHistory", param = NULL, date = date(),
               type = XCMSInternal$.PROCSTEP.PEAK.GROUPING(),
               fileIndex = seq_along(analyses(obj)),
               msLevel = msLevel)
    xdata@.processHistory <- c(xcms::processHistory(xdata), list(xph))

    return(xdata)
})

#' @rdname xcms-conv
#' @export
setMethod("getXCMSnExp", "featureGroupsXCMS3", function(obj, verbose, loadRawData)
{
    # first see if we can just return the embedded xcms object
    # NOTE: we can't do this if analyses have been subset

    anaInfo <- analysisInfo(obj)

    if (nrow(Biobase::pData(obj@xdata)) != length(anaInfo$analysis) ||
        !all(simplifyAnalysisNames(Biobase::pData(obj@xdata)$sample_name) == anaInfo$analysis))
        return(callNextMethod(obj, verbose = verbose, loadRawData = loadRawData))

    return(obj@xdata)
})

#' @rdname xcms-conv
#' @export
setMethod("getXCMSnExp", "featuresSet", function(obj, ..., set) getXCMSnExp(unset(obj, set), ...))

#' @rdname xcms-conv
#' @export
setMethod("getXCMSnExp", "featureGroupsSet", function(obj, ..., set) getXCMSnExp(unset(obj, set), ...))

loadXCMSRaw <- function(analyses, paths, cacheDB = NULL, verbose = TRUE)
{
    ret <- sapply(seq_along(analyses), function(anai)
    {
        p <- getMzMLOrMzXMLAnalysisPath(analyses[anai], paths[anai])
        hash <- makeFileHash(p)
        xr <- loadCacheData("EICData", hash, cacheDB)
        if (is.null(xr))
        {
            if (verbose)
                printf("Loading raw data from '%s'...\n", analyses[anai])
            xr <- xcms::xcmsRaw(p, profstep = 0)
            saveCacheData("EICData", xr, hash, cacheDB)
        }
        return(xr)
    })
    names(ret) <- analyses
    return(ret)
}

importXCMSPeaks <- function(peaks, analysisInfo)
{
    plist <- as.data.table(peaks)

    feat <- lapply(seq_len(nrow(analysisInfo)), function(sind)
    {
        ret <- plist[sample == sind]
        ret[, ID := seq_len(nrow(ret))]
        setnames(ret, c("rt", "rtmin", "rtmax", "maxo", "into"), c("ret", "retmin", "retmax", "intensity", "area"))
        return(ret[, intersect(XCMSFeatCols(), names(ret)), with = FALSE])
    })
    names(feat) <- analysisInfo$analysis

    return(feat)
}

XCMSFeatTblEqual <- function(tbl1, tbl2)
{
    return(isTRUE(all.equal(tbl1[, intersect(XCMSFeatCols(), names(tbl1)), with = FALSE],
                            tbl2[, intersect(XCMSFeatCols(), names(tbl2)), with = FALSE])))
}

# used by delete()
getKeptXCMSPeakInds <- function(old, new)
{
    newft <- featureTable(new)
    oldXCMSInds <- rbindlist(Map(featureTable(old), analyses(old), f = function(ft, ana)
    {
        if (nrow(ft) == 0)
            return(NULL)
        ret <- data.table(row = seq_len(nrow(ft)), analysis = ana)
        set(ret, j = "keep", value = if (!is.null(newft[[ana]]) && nrow(newft[[ana]]) > 0) ft$ID %in% newft[[ana]]$ID else FALSE)
        return(ret)
    }))
    oldXCMSInds[, inds := seq_len(nrow(oldXCMSInds))]
    return(oldXCMSInds[keep == TRUE]$inds)
}

# UNDONE: use with feature group plot EIC plotting
# UNDONE: update, some day...
# plotXCMSEIC <- function(EIC, rtRange = NULL, intRange = NULL, fillRange = NULL, retMin = FALSE,
#                         title = NULL, add = FALSE, col = "black", fillAlpha = 0.35, ...)
# {
#     fillCol <- adjustcolor(col, alpha.f = fillAlpha)
#
#     # plot limits
#     if (is.null(rtRange))
#         rtRange <- range(EIC$time) + c(-30, 30)
#     if (is.null(intRange))
#         intRange <- c(0, max(EIC$intensity))
#
#     if (is.null(title))
#         title <- sprintf("EIC plot m/z %s", paste0(attr(EIC, "mzRange"), collapse = " - "))
#
#     if (!add)
#         plot(0, type = "n", main = title, xlab = if (retMin) "Minutes" else "Seconds", ylab = "Intensity",
#              xlim = rtRange, ylim = intRange, ...)
#
#     points(if (retMin) EIC$time / 60 else EIC$time, EIC$intensity, type = "l", col = col)
#
#     if (!is.null(fillRange))
#     {
#         EICFill <- EIC[numGTE(EIC$time, fillRange[1]) & numLTE(EIC$time, fillRange[2]), ]
#         if (retMin)
#             EICFill$time <- EICFill$time / 60
#         polygon(c(EICFill$time, rev(EICFill$time)), c(EICFill$intensity, rep(0, length(EICFill$intensity))),
#                 col = fillCol, border = NA)
#     }
# }
