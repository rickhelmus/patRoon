#' @include main.R
#' @include features-openms.R
#' @include features-xcms.R
#' @include feature_groups.R
#' @include feature_groups-xcms.R

#' @rdname getXcmsSet
#' @export
setMethod("getXcmsSet", "features", function(obj, exportedData)
{
    # generate dummy XCMS set, based on https://groups.google.com/forum/m/#!topic/xcms/CGC0SKMVhAQ

    checkmate::assertFlag(exportedData)
    
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
            xr <- loadXCMSRaw(anaInfo$analysis[i], anaInfo$path[i])[[1]]
            rlist$raw[[i]] <- xr@scantime
            rlist$corrected[[i]] <- xr@scantime
        }
    }

    peaks(xs) <- as.matrix(do.call(function(...) rbind(..., make.row.names=F), plist))
    xs@rt <- rlist
    profinfo(xs) <- list(method = "bin", step = 0.1)

    return(xs)
})

#' @rdname getXcmsSet
#' @export
setMethod("getXcmsSet", "featuresOpenMS", function(obj, exportedData)
{
    return(callNextMethod(obj, TRUE))
})

#' @rdname getXcmsSet
#' @export
setMethod("getXcmsSet", "featuresXCMS", function(obj, exportedData)
{
    return(obj@xs)
})

#' @rdname getXcmsSet
#' @export
setMethod("getXcmsSet", "featureGroups", function(obj, exportedData)
{
    checkmate::assertFlag(exportedData)
    
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

    cat("Getting ungrouped xcmsSet...\n")
    xs <- getXcmsSet(getFeatures(obj), exportedData)

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

#' @rdname getXcmsSet
#' @export
setMethod("getXcmsSet", "featureGroupsXCMS", function(obj, exportedData)
{
    # first see if we can just return the xcmsSet used during grouping

    anaInfo <- analysisInfo(obj)

    if (length(filepaths(obj@xs)) != length(anaInfo$analysis) ||
        !all(simplifyAnalysisNames(filepaths(obj@xs)) == anaInfo$analysis))
        return(callNextMethod(obj, exportedData)) # files changed, need to update group statistics which is rather complex so just fallback

    return(obj@xs)
})

# export
# setMethod("checkEICs", "featureGroups", function(fGroups, EICWidth, prevState)
# {
#     fTable <- featureTable(fGroups)
#     anaInfo <- analysisInfo(fGroups)
#     gTable <- groups(fGroups)
#     gInfo <- groupInfo(fGroups)
#     ftindex <- groupFeatIndex(fGroups)
#
#     if (is.null(prevState))
#         selectInfo <- list(currentGroup = 1, currentFile = 1, enabledGroups = rep(TRUE, nrow(gInfo)),
#                            enabledFiles = rep(TRUE, nrow(anaInfo)), visibleGroups = rep(TRUE, nrow(gInfo)))
#     else
#         selectInfo <- prevState
#
#     # based on example from getGraphicsEvent()
#     interactivePlot <- function()
#     {
#         xrs <- sapply(seq_len(nrow(anaInfo)), function(i)
#         {
#             p <- getMzMLOrMzXMLAnalysisPath(anaInfo$analysis[i], anaInfo$path[i])
#             hash <- makeFileHash(p)
#             ret <- loadCacheData("checkEICs", hash)
#             if (is.null(ret))
#             {
#                 printf("Loading raw data of '%s'...\n", anaInfo$analysis[i])
#                 ret <- xcmsRaw(p, profstep = 0)
#                 saveCacheData("checkEICs", ret, hash)
#             }
#             return(ret)
#         })
#         names(xrs) <- anaInfo$analysis
#
#         unsgrps <- unique(anaInfo$group)
#         if (length(unsgrps) > 1)
#         {
#             # unique colour per sample group
#             gcols <- qualpal(length(unsgrps), "pretty_dark")$hex
#             cols <- gcols[match(anaInfo$group, unsgrps)]
#         }
#         else if (nrow(anaInfo) == 1) # qualpal needs n>1, just assign a color
#             cols <- "black"
#         else
#             cols <- qualpal(nrow(anaInfo), "pretty_dark")$hex # unique colour per analysis
#
#         dragInfo <- list(startx = 0, starty = 0, prevx = 0, prevy = 0, dx = 0, dy = 0, usr = NULL)
#
#
#         doPlot <- function(rtRange = NULL, intRange = NULL)
#         {
#             fts <- sapply(anaInfo$analysis, function(s)
#             {
#                 sind <- match(s, anaInfo$analysis)
#                 ftind <- ftindex[sind, selectInfo$currentGroup, with = F][[1]]
#                 if (ftind == 0)
#                     return(data.table(ret=NA, retmin=NA, retmax=NA, mz=NA))
#                 return(fTable[[s]][ftind])
#             }, simplify = F)
#
#             fts <- rbindlist(fts, fill=T)
#
#             rtPeakRange <- c(min(fts[, retmin], na.rm = T), max(fts[, retmax], na.rm = T))
#             rt <- mean(fts[, ret], na.rm = T)
#
#             if (is.null(rtRange))
#                 rtRange <- c(rtPeakRange[1] - 60, rtPeakRange[2] + 60)
#             if (is.null(intRange))
#                 intRange <- c(0, max(gTable[, selectInfo$currentGroup, with = F]) * 1.2)
#
#             par(oma = c(nrow(anaInfo)+2, 0, 0, 0))
#
#             if (any(selectInfo$enabledFiles))
#             {
#                 add <- FALSE
#                 for(i in seq_along(xrs))
#                 {
#                     if (selectInfo$enabledFiles[i])
#                     {
#                         plotFeatGroupEIC(xrs[[i]], rownames(gInfo)[selectInfo$currentGroup], gInfo[selectInfo$currentGroup, "mzs"], rt,
#                                          rtRange, intRange, mzwin = EICWidth, add=add, col=cols[i],
#                                          lwd = if (i == selectInfo$currentFile) 3 else 0.5)
#                         add <- TRUE
#                     }
#                 }
#             }
#             else
#                 plot(1, type = "n", xlab = "Seconds", ylab = "Intensity")
#
#             if (!selectInfo$enabledGroups[selectInfo$currentGroup])
#                 text(par("usr")[1] + diff(par("usr")[1:2])/2, par("usr")[3] + diff(par("usr")[3:4])/2,
#                      labels = "X", adj = c(0.5, 0.5), cex = 20, col = "red")
#
#             # if (plotMinutes)
#             #     rtPeakRange <- rtPeakRange / 60
#
#             # mark complete time range
#             # abline(v = rtPeakRange, col="red")
#
#             # mark range of currently selected file
#             ffind <- ftindex[selectInfo$currentFile, selectInfo$currentGroup, with=F][[1]]
#             if (ffind != 0)
#                 abline(v = fTable[[anaInfo$analysis[selectInfo$currentFile]]][ffind, c(retmin, retmax)], col="red", lty="dashed")
#
#             font <- rep(1, nrow(anaInfo))
#             font[selectInfo$currentFile] <- 2 # selection has bold font
#             nc <- ceiling(nrow(anaInfo) / 10)
#             cx <- min(1, 6.5 / (nrow(anaInfo) / nc))
#             txts <- abbreviate(anaInfo$analysis, 20, strict = T, method = "both.sides", use.classes = F)
#             leg <- sapply(seq_along(anaInfo$analysis), function(i)
#             {
#                 sel <- if (i == selectInfo$currentFile) "->" else "  "
#                 if (ftindex[i, selectInfo$currentGroup, with = F] != 0)
#                     sprintf("%s %s: rt: %.1f (%.1f-%.1f) m/z: %.4f int: %.2e",
#                             sel, txts[[i]], fts[i, ret], fts[i, retmin], fts[i, retmax], fts[i, mz], fts[i, intensity])
#                 else
#                     sprintf("%s %s: NA", sel, txts[[i]])
#             })
#
#             op <- par(no.readonly = TRUE)
#
#             par(family = "mono")
#             lc <- cols
#             lc[!selectInfo$enabledFiles] <- "grey"
#
#             # first call legend (without plotting) to get X coordinate, then use line2user to place legend in margin area
#             lx <- legend("bottom", legend=leg, col=cols, text.col=lc, lty=1, text.font=font, xpd = NA,
#                          cex = cx, ncol = nc, plot = FALSE)$rect$left
#             legend(lx, line2user(3.3, 1), legend=leg, col=cols, text.col=lc, lty=1, text.font=font, xpd = NA, cex = cx, ncol = nc)
#
#             # ktxt <- c("[left/right]  browse groups         [a/A] add DA EIC (no/with bg subtr)                [+/-]     zoom           [h] hide+disable group",
#             #           "[up/down]     select analysis       [e/E] add DA EICs all enabled (no/with bg subtr)   [left MB] drag chrom.    [p] hide disabled grps",
#             #           "[space/enter] toggle analysis/grp   [t]   toggle all analyses                          [r]       reset groups   [q] quit")
#
#             ktxt <- c("[left/right]  browse groups         [a/A] add DA EIC (no/with bg subtr)                [+/-]     zoom",
#                       "[up/down]     select analysis       [e/E] add DA EICs all enabled (no/with bg subtr)   [left MB] pan chrom",
#                       "[space/enter] toggle analysis/grp   [h]   hide+disable group                           [r]       reset groups",
#                       "[t]           toggle all analyses   [p]   hide disabled grps                           [q]       quit")
#
#             cx <- 0.63
#             at <- 0.5 - (max(strwidth(ktxt, "figure", cex = cx))/2)
#             # at <- 0.5 - (grconvertX(max(strwidth(ktxt, cex = cx)), to="npc")/2)
#             par(mex = cx)
#             for (i in seq_along(ktxt))
#                 mtext(ktxt[i], 1, line = par("oma")[1]-1.5-length(ktxt)+i, TRUE, col = "blue", cex = cx,
#                       at = at, adj = 0)
#
#             par(op)
#         }
#
#         devset <- function()
#         {
#             if (dev.cur() != eventEnv$which)
#                 dev.set(eventEnv$which)
#         }
#
#         setNextGroup <- function(rev = FALSE)
#         {
#             if (!any(selectInfo$visibleGroups))
#                 return()
#
#             repeat
#             {
#                 if (rev)
#                     selectInfo$currentGroup <<- if (selectInfo$currentGroup == 1) nrow(gInfo) else selectInfo$currentGroup - 1
#                 else
#                     selectInfo$currentGroup <<- if (selectInfo$currentGroup == nrow(gInfo)) 1 else selectInfo$currentGroup + 1
#
#                 if (selectInfo$visibleGroups[selectInfo$currentGroup])
#                     break
#             }
#         }
#
#         mouseDown <- function(buttons, x, y)
#         {
#             dragInfo$startx <<- x
#             dragInfo$starty <<- y
#             dragInfo$prevx <<- dragInfo$prevy <- 0
#             dragInfo$usr <<- par("usr")
#             devset()
#             eventEnv$onMouseMove <- mouseMove
#             return(NULL)
#         }
#
#         mouseUp <- function(buttons, x, y)
#         {
#             eventEnv$onMouseMove <- NULL
#             return(NULL)
#         }
#
#         mouseMove <- function(buttons, x, y)
#         {
#             devset()
#             deltax <- diff(grconvertX(c(dragInfo$startx, x), "ndc", "user"))
#             deltay <- diff(grconvertY(c(dragInfo$starty, y), "ndc", "user"))
#             if (abs(deltax-dragInfo$prevx) + abs(deltay-dragInfo$prevy) > 0)
#             {
#                 dragInfo$prevx <<- deltax
#                 dragInfo$prevy <<- deltay
#                 doPlot(dragInfo$usr[1:2] - deltax, dragInfo$usr[3:4] - deltay)
#             }
#
#             return(NULL)
#         }
#
#         keyDown <- function(key)
#         {
#             rtr <- intr <- NULL
#             if (is.null(key) || key == "q" || key == "Q") # quit
#                 return(1)
#             else if (key == "Left") # browse group
#                 setNextGroup(rev = TRUE)
#             else if (key == "Right") # browse group
#                 setNextGroup()
#             else if (key == "Up") # select file
#                 selectInfo$currentFile <<- if (selectInfo$currentFile == 1) nrow(anaInfo) else selectInfo$currentFile - 1
#             else if (key == "Down") # select file
#                 selectInfo$currentFile <<- if (selectInfo$currentFile == nrow(anaInfo)) 1 else selectInfo$currentFile + 1
#             else if (key == "+" || key == "-") # zoom
#             {
#                 p <- par("usr")
#
#                 w <- diff(p[1:2])
#                 sw <- w * if (key == "+") 0.85 else 1.1
#                 offset <- (sw - w) / 2
#                 rtr <- p[1:2] + c(-offset, offset)
#
#                 w <- diff(p[3:4])
#                 sw <- w * if (key == "+") 0.85 else 1.1
#                 offset <- (sw - w) / 2
#                 intr <- p[3:4] + c(-offset, offset)
#             }
#             else if (key == " ") # toggle analysis
#                 selectInfo$enabledFiles[selectInfo$currentFile] <<- !selectInfo$enabledFiles[selectInfo$currentFile]
#             else if (key == "ctrl-J") # enter, toggle group
#                 selectInfo$enabledGroups[selectInfo$currentGroup] <<- !selectInfo$enabledGroups[selectInfo$currentGroup]
#             else if (key == "a" || key == "A") # add DA trace
#                 addDAEIC(anaInfo$analysis[selectInfo$currentFile], anaInfo$path[selectInfo$currentFile],
#                          gInfo[selectInfo$currentGroup, "mzs"], EICWidth, bgsubtr = (key == "A"))
#             else if (key == "e" || key == "E") # add DA trace enabled files
#             {
#                 for (f in anaInfo$analysis[selectInfo$enabledFiles])
#                     addDAEIC(f, dpath, gInfo[selectInfo$currentGroup, "mzs"], EICWidth, bgsubtr = (key == "E"))
#             }
#             else if (key == "t") # toggle all analysis
#                 selectInfo$enabledFiles <<- rep(mean(selectInfo$enabledFiles) < 0.5, nrow(anaInfo))
#             else if (key == "r") # reset groups
#             {
#                 selectInfo$enabledGroups <<- rep(TRUE, nrow(gInfo))
#                 selectInfo$visibleGroups <<- rep(TRUE, nrow(gInfo))
#             }
#             else if (key == "h") # hide group
#             {
#                 selectInfo$enabledGroups[selectInfo$currentGroup] <<- FALSE
#                 selectInfo$visibleGroups[selectInfo$currentGroup] <<- FALSE
#                 setNextGroup()
#             }
#             else if (key == "p") # purge groups
#             {
#                 selectInfo$visibleGroups <<- selectInfo$visibleGroups & selectInfo$enabledGroups
#                 if (!selectInfo$visibleGroups[selectInfo$currentGroup])
#                     setNextGroup()
#             }
#             else
#                 return(NULL)
#
#             devset()
#             doPlot(rtr, intr)
#
#             return(NULL)
#         }
#
#         setGraphicsEventHandlers(onMouseDown = mouseDown, onMouseUp = mouseUp, onKeybd = keyDown)
#         eventEnv <- getGraphicsEventEnv()
#         doPlot()
#     }
#
#     oldpar <- par(ask = FALSE)
#     x11(width = 9)
#     wdev <- dev.cur()
#     interactivePlot()
#     getGraphicsEvent()
#     dev.off(wdev)
#     par(oldpar)
#
#     invisible(selectInfo)
# })

# based on plotEIC from xcms
plotFeatGroupEIC <- function(xr, groupName, mz, rt, rtRange, intRange, mzwin=0.005, plotMinutes=F, add=F, ...)
{
    if (rtRange[1] < xr@scantime[length(xr@scantime)-1] && rtRange[2] > xr@scantime[2])
    {
        # only plot something when retention range is within analysis range
        # note: clamp between second and second last datapoint so that at least 1 point can be drawn

        mzr <- c(mz - mzwin, mz + mzwin)

        EIC <-  rawEIC(xr, mzrange = mzr, rtrange = rtRange)
        points <- cbind(xr@scantime[EIC$scan], EIC$intensity)

        if (plotMinutes)
            points[1] <- points[1] / 60
    }
    else
        points <- NULL # nothing to draw...

    if(add)
        points(points, type="l", ...)
    else
    {
        plot(points, type = "l", main = sprintf("Group '%s' - rt: %.1f - m/z: %.4f", groupName, rt, mz),
             xlab = "", ylab = "Intensity", xlim = if (plotMinutes) rtRange / 60 else rtRange, ylim = intRange, ...)
        # draw X axis later to position it closer
        title(xlab = if (plotMinutes) "Minutes" else "Seconds", line = 1.9)
    }


    return()
}

loadXCMSRaw <- function(analyses, paths, cacheDB = NULL)
{
    ret <- sapply(seq_along(analyses), function(anai)
    {
        p <- getMzMLOrMzXMLAnalysisPath(analyses[anai], paths[anai])
        hash <- makeFileHash(p)
        xr <- loadCacheData("EICData", hash, cacheDB)
        if (is.null(xr))
        {
            printf("Loading raw data from '%s'...\n", analyses[anai])
            xr <- xcmsRaw(p, profstep = 0)
            saveCacheData("EICData", xr, hash, cacheDB)
        }
        return(xr)
    })
    names(ret) <- analyses
    return(ret)
}

getXCMSEiCHashFromPath <- function(path, rtRange, mzRange) makeHash(makeFileHash(path), rtRange, mzRange)

loadXCMSEICFromCache <- function(hash, cacheDB = NULL)
{
    ret <- loadCacheData("XCMSEIC", hash, cacheDB)
    return(ret)
}

loadXCMSEICFromRaw <- function(xr, rtRange, mzRange, hash = NULL, cacheDB = NULL)
{
    if (is.null(hash))
        hash <- getXCMSEiCHashFromPath(xr@filepath, rtRange, mzRange)

    ret <- loadXCMSEICFromCache(hash, cacheDB)

    if (is.null(ret))
    {
        rtRange <- c(max(0, rtRange[1]), min(xr@scantime[length(xr@scantime)], rtRange[2]))
        rEIC <- rawEIC(xr, mzrange = mzRange, rtrange = rtRange)
        ret <- data.frame(time = xr@scantime[rEIC$scan], intensity = rEIC$intensity)
        saveCacheData("XCMSEIC", ret, hash, cacheDB)

        # UNDONE: keep this?
        attr(ret, "rtRange") <- rtRange
        attr(ret, "mzRange") <- mzRange
    }

    return(ret)
}

# loads EICs for all features within a feature groups object.
loadXCMSEICForFGroups <- function(fGroups, rtWindow, mzWindow, topMost, onlyPresent)
{
    gNames <- names(fGroups)
    gInfo <- groupInfo(fGroups)
    ftind <- groupFeatIndex(fGroups)
    gTable <- groups(fGroups)
    fTable <- featureTable(fGroups)
    anaInfo <- analysisInfo(fGroups)

    anaHashes <- list()
    xrs <- list()

    if (!is.null(topMost))
        topMost <- min(topMost, nrow(anaInfo))

    cacheDB <- openCacheDB()

    EICs <- lapply(seq_len(length(fGroups)), function(grpi)
    {
        if (!is.null(topMost)) #&& sum(!is.na(fts$intensity)) > topMost)
        {
            oint <- order(-gTable[[grpi]])
            topAnalysesInd <- oint[seq_len(topMost)]
        }
        else
            topAnalysesInd <- seq_len(nrow(anaInfo))

        if (onlyPresent)
            topAnalysesInd <- topAnalysesInd[gTable[[grpi]][topAnalysesInd] != 0]

        # rtRange <- gInfo[grpi, "rts"] + c(-rtWindow, rtWindow)
        mzRange <- gInfo[grpi, "mzs"] + c(-mzWindow, mzWindow)
        
        getFTCol <- function(anai, col)
        {
            ana <- anaInfo$analysis[anai]
            if (ftind[[grpi]][anai] == 0)
                NA
            else
                fTable[[anaInfo$analysis[anai]]][ftind[[grpi]][anai]][[col]]
        }

        rtRange <- c(min(sapply(topAnalysesInd,
                                function(anai) getFTCol(anai, "retmin")),
                         na.rm = TRUE),
                     max(sapply(topAnalysesInd,
                                function(anai) getFTCol(anai, "retmax")),
                         na.rm = TRUE))
        if (any(is.na(rtRange)))
            rtRange <- gInfo[grpi, "rts"]
        rtRange <- rtRange + c(-rtWindow, rtWindow)
        
        ret <- lapply(topAnalysesInd, function(anai)
        {
            ana <- anaInfo$analysis[anai]

            # feat <- fTable[[ana]][ftind[[grpi]][anai]]
            # rtRange <- c(feat$retmin - rtWindow, feat$retmax + rtWindow)

            if (is.null(anaHashes[[ana]]))
                anaHashes[[ana]] <<- makeFileHash(getMzMLOrMzXMLAnalysisPath(ana, anaInfo$path[anai]))

            hash <- makeHash(anaHashes[[ana]], rtRange, mzRange)

            eic <- loadXCMSEICFromCache(hash, cacheDB)
            if (is.null(eic))
            {
                if (is.null(xrs[[ana]]))
                    xrs[[ana]] <<- loadXCMSRaw(ana, anaInfo$path[anai], cacheDB)[[1]]

                eic <- loadXCMSEICFromRaw(xrs[[ana]], rtRange, mzRange, hash, cacheDB)
            }
            return(eic)
        })
        names(ret) <- anaInfo$analysis[topAnalysesInd]
        return(ret)
    })
    names(EICs) <- gNames

    closeCacheDB(cacheDB)

    return(EICs)
}

# UNDONE: use with feature group plot EIC plotting
plotXCMSEIC <- function(EIC, rtRange = NULL, intRange = NULL, fillRange = NULL, retMin = FALSE,
                        title = NULL, add = FALSE, col = "black", fillAlpha = 0.35, ...)
{
    fillCol <- adjustcolor(col, alpha.f = fillAlpha)

    # plot limits
    if (is.null(rtRange))
        rtRange <- range(EIC$time) + c(-30, 30)
    if (is.null(intRange))
        intRange <- c(0, max(EIC$intensity))

    if (is.null(title))
        title <- sprintf("EIC plot m/z %s", paste0(attr(EIC, "mzRange"), collapse = " - "))

    if (!add)
        plot(0, type = "n", main = title, xlab = if (retMin) "Minutes" else "Seconds", ylab = "Intensity",
             xlim = rtRange, ylim = intRange, ...)

    points(if (retMin) EIC$time / 60 else EIC$time, EIC$intensity, type = "l", col = col)

    if (!is.null(fillRange))
    {
        EICFill <- EIC[EIC$time >= fillRange[1] & EIC$time <= fillRange[2], ]
        if (retMin)
            EICFill$time <- EICFill$time / 60
        polygon(c(EICFill$time, rev(EICFill$time)), c(EICFill$intensity, rep(0, length(EICFill$intensity))),
                col = fillCol, border = NA)
    }
}
