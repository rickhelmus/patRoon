#' @include main.R
#' @include features-openms.R
#' @include features-xcms.R
#' @include feature_groups.R
#' @include feature_groups-xcms.R

#' @rdname getXCMSSet
#' @export
setMethod("getXCMSSet", "features", function(obj, exportedData, verbose)
{
    # generate dummy XCMS set, based on https://groups.google.com/forum/m/#!topic/xcms/CGC0SKMVhAQ

    checkmate::assertFlag(exportedData)
    checkmate::assertFlag(verbose)

    xs <- new(getClassDef("xcmsSet", package = "xcms"))
    anaInfo <- analysisInfo(obj)
    phenoData(xs) <- data.frame(class = anaInfo$group, row.names = anaInfo$analysis)

    if (exportedData)
        filepaths(xs) <- sapply(seq_len(nrow(anaInfo)), function(i) getMzMLOrMzXMLAnalysisPath(anaInfo$analysis[i], anaInfo$path[i]), USE.NAMES = F)
    else
        filepaths(xs) <- anaInfo$analysis # dummy paths

    fts <- featureTable(obj)
    plist <- list()
    rlist <- list(raw = vector("list", length(anaInfo$analysis)), corrected = vector("list", length(anaInfo$analysis)))

    for (i in seq_len(nrow(anaInfo)))
    {
        ft <- fts[[anaInfo$analysis[i]]]

        if (nrow(ft) > 0)
            plist[[i]] <- data.frame(mz = ft$mz, mzmin = ft$mzmin, mzmax = ft$mzmax, rt = ft$ret,
                                     rtmin = ft$retmin, rtmax = ft$retmax, maxo = ft$intensity, into = ft$area,
                                     sample = i, stringsAsFactors = F)
        else
            plist[[i]] <- data.frame(mz = numeric(), mzmin = numeric(), mzmax = numeric(), rt = numeric(),
                                     rtmin = numeric(), rtmax = numeric(), maxo = numeric(), into = numeric(),
                                     sample = numeric())

        if (exportedData)
        {
            xr <- loadXCMSRaw(anaInfo$analysis[i], anaInfo$path[i], verbose = verbose)[[1]]
            rlist$raw[[i]] <- xr@scantime
            rlist$corrected[[i]] <- xr@scantime
        }
    }

    peaks(xs) <- as.matrix(do.call(function(...) rbind(..., make.row.names=F), plist))
    xs@rt <- rlist
    profinfo(xs) <- list(method = "bin", step = 0.1)

    return(xs)
})

#' @rdname getXCMSSet
#' @export
setMethod("getXCMSSet", "featuresOpenMS", function(obj, exportedData, verbose)
{
    return(callNextMethod(obj, TRUE, verbose = verbose))
})

#' @rdname getXCMSSet
#' @export
setMethod("getXCMSSet", "featuresXCMS", function(obj, exportedData, verbose)
{
    return(obj@xs)
})

#' @rdname getXCMSSet
#' @export
setMethod("getXCMSSet", "featureGroups", function(obj, exportedData, verbose = TRUE)
{
    checkmate::assertFlag(exportedData)
    checkmate::assertFlag(verbose)

    fTable <- featureTable(obj)

    # add unique feat IDs (corresponding to ftindex)
    combfts <- lapply(copy(fTable), function(f) f[, ftID := seq_len(.N)])
    names(combfts) <- NULL # get rid of file names ...
    combfts <- rbindlist(combfts, idcol = "sind") # ... so that rbindlist will generate numeric IDs
    # combfts should now be aligned with peak table of xcms set

    combfts[, xsID := seq_len(.N)] # ID corresponding to XCMS peaks

    gTable <- groups(obj)
    ftind <- copy(groupFeatIndex(obj))
    anaInfo <- analysisInfo(obj)

    if (verbose)
        cat("Getting ungrouped xcmsSet...\n")
    xs <- getXCMSSet(getFeatures(obj), exportedData, verbose = verbose)

    if (verbose)
        cat("Making group index table\n")

    setkeyv(combfts, c("sind", "ftID"))
    groupidx(xs) <- lapply(seq_along(ftind), function(g)
    {
        fti <- data.table(ftID = ftind[[g]], sind = seq_along(ftind[[g]]))
        setkeyv(fti, c("sind", "ftID"))
        # join by sample ID and by feat ID so that we end up with a dt having only
        # relevant rows for the IDs given in ftind
        combfts[fti, nomatch = 0][["xsID"]] # xsID is the XCMS peak ID
    })

    if (verbose)
        cat("Making groups...\n")

    # generate xcmsSet group matrix
    grps <- matrix(unlist(sapply(seq_along(ftind), function(g)
    {
        fprops <- matrix(unlist(sapply(seq_len(nrow(anaInfo)), function(s)
        {
            fi <- ftind[[g]][s]
            if (fi == 0)
                list(NA, NA)
            else
                fTable[[anaInfo$analysis[s]]][fi, c("mz", "ret")]
        }), use.names = FALSE), ncol = 2, byrow = TRUE, dimnames = list(NULL, c("mz", "ret")))

        list(mzmed = median(fprops[, "mz"], na.rm = TRUE), mzmin = min(fprops[, "mz"], na.rm = TRUE),
             mzmax = max(fprops[, "mz"], na.rm = TRUE), rtmed = median(fprops[, "ret"], na.rm = TRUE),
             rtmin = min(fprops[, "ret"], na.rm = TRUE), rtmax = max(fprops[, "ret"], na.rm = TRUE),
             npeaks = sum(ftind[[g]] != 0))
    }, simplify = FALSE), use.names = FALSE), ncol = 7, byrow = TRUE, dimnames = list(NULL, c("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax", "npeaks")))

    # add counts per sample group
    xcms::groups(xs) <- cbind(grps, sapply(unique(anaInfo$group), function(sg)
    {
        sapply(seq_along(ftind), function(g) sum(ftind[[g]][sg == anaInfo$group] != 0))
    }))

    return(xs)
})

#' @rdname getXCMSSet
#' @export
setMethod("getXCMSSet", "featureGroupsXCMS", function(obj, exportedData, verbose)
{
    # first see if we can just return the xcmsSet used during grouping

    anaInfo <- analysisInfo(obj)

    if (length(filepaths(obj@xs)) != length(anaInfo$analysis) ||
        !all(simplifyAnalysisNames(filepaths(obj@xs)) == anaInfo$analysis))
        return(callNextMethod(obj, exportedData, verbose = verbose)) # files changed, need to update group statistics which is rather complex so just fallback

    return(obj@xs)
})

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
            xr <- xcmsRaw(p, profstep = 0)
            saveCacheData("EICData", xr, hash, cacheDB)
        }
        return(xr)
    })
    names(ret) <- analyses
    return(ret)
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
