#' @include main.R
NULL

#' Plotting of grouped features
#'
#' Various plotting functions for feature group data.
#'
#' @param retMin Plot retention time in minutes (instead of seconds).
#' @param pch,type,lty Common plotting parameters passed to \emph{e.g.} \code{\link[graphics]{plot}}. For \code{plot}:
#'   if \code{pch=NULL} then values are automatically assigned.
#' @param col Colour(s) used. If \code{col=NULL} then colours are automatically generated.
#' @param colourBy Sets the automatic colour selection: \code{"none"} for a single colour or
#'   \code{"rGroups"}/\code{"fGroups"} for a distinct colour per replicate/feature group.
#' @param showLegend Plot a legend if \code{TRUE}.
#' @param which A character vector with replicate groups used for comparison. Set to \code{NULL} to ignore.
#'
#'   For \code{plotVenn}: alternatively a named \code{list} containing elements of \code{character} vectors with
#'   replicate groups to compare. For instance, \code{which=list(infl = c("influent-A", "influent-B"), effl =
#'   c("effluent-A", "effluent-B"))}, will compare the features in replicate groups \samp{"influent-A/B"} against those
#'   in \samp{"effluent-A/B"}. The names of the list are used for labelling in the plot, and will be made automatically
#'   if not specified.
#' @param \dots passed to \code{\link[base]{plot}} (\code{plot} and \code{plotChroms}), \code{\link[graphics]{lines}}
#'   (\code{plotInt}), \pkg{\link{VennDiagram}} plotting functions (\code{plotVenn}), \code{\link{chordDiagram}}
#'   (\code{plotChord}) or \code{\link[UpSetR]{upset}} (\code{plotUpSet}).
#' @param sets \setsWF For \code{plotInt}: if \code{TRUE} then feature intensities are plot per set (order follows the
#'   \link[=analysis-information]{analysis information}).
#'
#'   For \code{plotVenn}: If \code{TRUE} then the \code{which} argument changes its meaning and is used to specify the
#'   names of the sets to be compared.
#'
#' @inheritParams featureGroups-class
#'
#' @section Sets workflows: \setsWFChangedMethods{
#'
#'   \item \code{plotVenn} and \code{plotInt} allow to handle data per set. See the \code{sets} argument description.
#'
#'   }
#'
#' @seealso \code{\link{featureGroups-class}}, \code{\link{groupFeatures}}
#'
#' @name feature-plotting
NULL

#' @details \code{plot} Generates an \emph{m/z} \emph{vs} retention time
#'   plot for all featue groups. Optionally highlights unique/overlapping
#'   presence amongst replicate groups.
#' @param onlyUnique If \code{TRUE} and \code{colourBy="rGroups"} then only
#'   feature groups that are unique to a replicate group are plotted.
#' 
#' @rdname feature-plotting
#' @export
setMethod("plot", c(x = "featureGroups", y = "missing"), function(x, colourBy = c("none", "rGroups", "fGroups"),
                                                                  onlyUnique = FALSE, retMin = FALSE,
                                                                  showLegend = TRUE, col = NULL,
                                                                  pch = NULL, ...)
{
    rGroups <- replicateGroups(x)
    if (is.null(which))
        which <- rGroups
    
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertFlag, . ~ onlyUnique + retMin + showLegend, fixed = list(add = ac))
    colourBy <- checkmate::matchArg(colourBy, c("none", "rGroups", "fGroups"), add = ac)
    checkmate::reportAssertions(ac)
    
    if (length(x) == 0)
        noDataPlot()
    else
    {
        if (colourBy == "fGroups" || colourBy == "none")
        {
            if (is.null(col))
                col <- if (colourBy == "fGroups") colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(x)) else "black"
            if (is.null(pch))
                pch <- 16
            
            if (colourBy == "fGroups" && showLegend)
            {
                labels <- names(x)
                labCol <- rep(col, length.out = length(labels))
                labPch <- rep(pch, length.out = length(labels))
                names(labCol) <- labels; names(labPch) <- labels
            }
            else
                showLegend <- FALSE
        }
        else if (colourBy == "rGroups")
        {
            labels <- c(replicateGroups(x), "overlap")
            
            if (is.null(col))
                labCol <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(labels))
            else
                labCol <- rep(col, length.out = length(labels))
            names(labCol) <- labels
            
            if (is.null(pch))
            {
                # prefer closed symbols (15-25)
                ll <- length(labels)
                if (ll < (25-15))
                    labPch <- seq_len(ll) + 14
                else if (ll <= 25)
                    labPch <- seq_len(ll)
                else
                    labPch <- rep(16, ll) # just stick with one
            }
            else
                labPch <- rep(pch, length.out = length(labels))
            names(labPch) <- labels
            
            # get averaged intensities for each rGroup and omit initial name/rt/mz columns
            gTable <- as.data.table(x, average = TRUE)[, replicateGroups(x), with = FALSE]
            
            for (r in seq_len(nrow(gTable)))
            {
                present <- which(gTable[r] != 0)
                set(gTable, r, "label", if (length(present) > 1) "overlap" else labels[present])
            }
            
            if (onlyUnique)
            {
                keep <- gTable$label != "overlap"
                gTable <- gTable[keep]; x <- x[, keep]
            }
            
            # remove unassigned (e.g. rGroups without unique feature groups)
            labels <- labels[labels %in% gTable$label]
            
            col <- labCol[gTable$label]; pch <- labPch[gTable$label]
        }
        
        oldp <- par(no.readonly = TRUE)
        if (showLegend)
        {
            makeLegend <- function(x, y, ...)
            {
                return(legend(x, y, labels, col = labCol[labels], pch = labPch[labels],
                              text.col = labCol[labels],
                              xpd = NA, ncol = 1, cex = 0.75, bty = "n", ...))
            }
            
            plot.new()
            leg <- makeLegend(0, 0, plot = FALSE)
            lw <- (grconvertX(leg$rect$w, to = "ndc") - grconvertX(0, to = "ndc"))
            par(omd = c(0, 1 - lw, 0, 1), new = TRUE)
        }
        
        plot(if (retMin) x@groupInfo$rts / 60 else x@groupInfo$rts, x@groupInfo$mzs,
             xlab = if (retMin) "Retention time (min)" else "Retention time (sec.)",
             ylab = "m/z", col = col, pch = pch, ...)
        
        if (showLegend)
            makeLegend(par("usr")[2], par("usr")[4])
        
        par(oldp)
    }
})

#' @details \code{plotInt} Generates a line plot for the (averaged) intensity
#'   of feature groups within all analyses
#' @param xnames Plot analysis (or replicate group if \code{average=TRUE}) names on the x axis.
#' 
#' @rdname feature-plotting
#' @export
setMethod("plotInt", "featureGroups", function(obj, average = FALSE, normalized = FALSE, xnames = TRUE,
                                               showLegend = FALSE, pch = 20, type = "b", lty = 3, col = NULL, ...)
{
    aapply(checkmate::assertFlag, . ~ average + normalized + xnames + showLegend)
    doPlotFeatInts(obj, average, normalized, xnames, showLegend, pch, type, lty, col, ..., doSets = FALSE)    
})

#' @rdname feature-plotting
#' @export
setMethod("plotInt", "featureGroupsSet", function(obj, average = FALSE, normalized = FALSE, xnames = !sets,
                                                  showLegend = sets, pch = 20, type = "b", lty = 3, col = NULL, ...,
                                                  sets = FALSE)
{
    aapply(checkmate::assertFlag, . ~ average + normalized + xnames + showLegend + sets)
    doPlotFeatInts(obj, average, normalized, xnames, showLegend, pch, type, lty, col, ..., doSets = sets)    
})

setMethod("plotIntHash", "featureGroups", function(obj, average = FALSE, ...) makeHash(allArgs()))

#' @details \code{plotChord} Generates a chord diagram which can be used to
#'   visualize shared presence of feature groups between analyses or replicate
#'   groups. In addition, analyses/replicates sharing similar properties
#'   (\emph{e.g.} location, age, type) may be grouped to enhance visualization
#'   between these 'outer groups'.
#'
#' @param outerGroups Character vector of names to be used as outer groups. The
#'   values in the specified vector should be named by analysis names
#'   (\code{average} set to \code{FALSE}) or replicate group names
#'   (\code{average} set to \code{TRUE}), for instance: \code{c(analysis1 =
#'   "group1", analysis2 = "group1", analysis3 = "group2")}. Set to \code{NULL}
#'   to disable outer groups.
#' @param addIntraOuterGroupLinks If \code{TRUE} then links will be added within
#'   outer groups.
#'
#' @template plotChord-args
#'
#' @references \addCitations{circlize}{1}
#'
#' @rdname feature-plotting
#' @export
setMethod("plotChord", "featureGroups", function(obj, addSelfLinks = FALSE, addRetMzPlots = TRUE, average = FALSE,
                                                 outerGroups = NULL, addIntraOuterGroupLinks = FALSE, ...)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertFlag, . ~ addSelfLinks + addRetMzPlots + average + addIntraOuterGroupLinks,
           fixed = list(add = ac))
    checkmate::assertCharacter(outerGroups, min.chars = 1, min.len = 2, names = "unique", null.ok = TRUE)
    checkmate::reportAssertions(ac)
    
    if (length(obj) == 0)
        stop("Can't plot empty feature groups object")
    
    hasOuter <- !is.null(outerGroups)
    
    obj <- removeEmptyAnalyses(obj)
    
    anaInfo <- analysisInfo(obj)
    gInfo <- groupInfo(obj)
    
    if (average)
        snames <- unique(anaInfo$group)
    else
        snames <- anaInfo$analysis
    
    if (length(snames) < 2)
        stop(sprintf("Nothing to compare: need multiple (non-empty) %s!", if (average) "replicate groups" else "analyses"))
    
    if (hasOuter && !all(snames %in% names(outerGroups)))
        stop(sprintf("The following %s have no outerGroups assigned: %s", if (average) "replicate groups" else "analyses",
                     paste0(setdiff(snames, names(outerGroups)), collapse = ", ")))
    
    nsamp <- length(snames)
    
    chordTable <- rbindlist(lapply(seq_along(snames),
                                   function(sni) data.table(from = snames[sni], to = snames[seq(sni, length(snames))])))
    
    getGTable <- function(snlist = c())
    {
        if (average)
        {
            if (length(snlist) > 0)
                fgf <- replicateGroupFilter(obj, snlist, verbose = FALSE)
            else
                fgf <- obj
            if (length(fgf) == 0)
                return(data.table())
            return(averageGroups(fgf))
        }
        else
        {
            if (length(snlist) > 0)
                return(groupTable(obj[snlist]))
            return(groupTable(obj))
        }
    }
    
    getLinkScore <- function(sn1, sn2)
    {
        if (sn1 == sn2)
            return(0)
        
        gTable <- getGTable(c(sn1, sn2))
        if (nrow(gTable) == 0)
            return(0)
        
        # Count all feature groups that are present in both samples/groups
        return(sum(gTable[, sapply(.SD, function(rows) all(rows > 0))]))
    }
    
    chordTable[, value := as.integer(Vectorize(getLinkScore)(from, to))]
    
    if (addSelfLinks)
    {
        gt <- getGTable()
        uniqueLinkCount <- sapply(seq_along(snames),
                                  function(sni) sum(sapply(gt, function(ints) ints[sni] > 0 && all(ints[-sni] == 0))))
        chordTable[from == to, value := uniqueLinkCount[.GRP], by = from]
    }
    
    if (hasOuter)
    {
        chordTable[, groupFrom := outerGroups[from]]
        chordTable[, groupTo := outerGroups[to]]
        if (!addIntraOuterGroupLinks)
            chordTable[groupFrom == groupTo & from != to, value := 0] # clear links within same groups (except self links)
        setorder(chordTable, groupFrom)
        
        remainingSN <- unique(unlist(chordTable[value != 0, .(from, to)])) # assigned samples, others will be removed
        og <- outerGroups[remainingSN] # outer groups assigned to each remaining sample
        gaps <- rep(1, length(og)) # initialize gaps
        gaps[cumsum(sapply(unique(og), function(x) length(og[og == x])))] <- 8 # make gap bigger after each outer group
        circlize::circos.par(gap.after = gaps)
    }
    
    if (all(chordTable$value == 0))
        stop("Did not found any overlap! Nothing to plot.")
    
    tracks <- NULL
    if (hasOuter)
        tracks <- list(list(track.height = 0.1, track.margin = c(if (addRetMzPlots) 0.05 else 0.06, 0)))
    if (addRetMzPlots)
        tracks <- c(tracks, list(list(track.height = 0.1, track.margin = c(0.08, 0))))
    
    maxv <- max(if (hasOuter) chordTable[groupFrom != groupTo, value] else chordTable$value)
    colFunc <- circlize::colorRamp2(maxv * seq(0, 1, 0.25),
                                    c("blue4", "deepskyblue1", "green", "orange", "red"),
                                    transparency = 0.5)
    
    if (hasOuter && addIntraOuterGroupLinks)
    {
        colFuncWithin <- circlize::colorRamp2(range(chordTable[groupFrom == groupTo, value]),
                                              c("grey80", "grey60"), transparency = 0.7)
        linkColors <- chordTable[, ifelse(groupFrom == groupTo, colFuncWithin(value), colFunc(value))]
    }
    else
        linkColors <- chordTable[, colFunc(value)]
    
    cdf <- circlize::chordDiagram(chordTable[, 1:3], annotationTrack = c("grid", "axis"),
                                  preAllocateTracks = tracks,
                                  grid.col = getBrewerPal(nsamp, "Dark2"),
                                  col = linkColors,
                                  annotationTrackHeight = c(0.1, 0.09),
                                  ...)
    
    circlize::circos.track(track.index = length(tracks) + 1, panel.fun = function(x, y)
    {
        sector.index <- circlize::get.cell.meta.data("sector.index")
        xlim <- circlize::get.cell.meta.data("xlim")
        ylim <- circlize::get.cell.meta.data("ylim")
        circlize::circos.text(mean(xlim), mean(ylim), sector.index, col = "white", cex = 1, niceFacing = TRUE)
    }, bg.border = NA)
    
    if (addRetMzPlots)
    {
        retMz <- rbindlist(sapply(unique(cdf$rn), function(sn)
        {
            ftgrps <- colnames(getGTable(sn))
            return(gInfo[ftgrps, ])
        }, simplify = FALSE), idcol = "sname")
        retMz$rts <- retMz$rts / max(retMz$rts) # normalize
        
        circlize::circos.track(fa = retMz$sname, x = retMz$rts, y = retMz$mzs, ylim = c(0, max(retMz$mzs)), track.index = length(tracks),
                               panel.fun = function(x, y)
                               {
                                   x <- x / (max(x) / circlize::get.cell.meta.data("xrange"))
                                   circlize::circos.points(x, y, cex = 0.5, col = "blue", pch = 16)
                               })
    }
    
    if (hasOuter)
    {
        finalChordTable <- chordTable[from %in% cdf$rn]
        finalOuterGroups <- unique(finalChordTable$groupFrom)
        
        ogcol <- getBrewerPal(length(finalOuterGroups), "Paired")
        for (ogi in seq_along(finalOuterGroups))
            circlize::highlight.sector(unique(finalChordTable[groupFrom == finalOuterGroups[ogi], from]), track.index = 1, col = ogcol[ogi],
                                       text = finalOuterGroups[ogi], cex = 1, text.col = "white", niceFacing = TRUE)
    }
    
    circlize::circos.clear()
})

#' @details \code{plotChroms} Plots extracted ion chromatograms (EICs) of feature groups.
#'
#' @param rtWindow Retention time (in seconds) that will be subtracted/added to respectively the minimum and maximum
#'   retention time of the plotted feature groups. Thus, setting this value to a positive value will 'zoom out' on the
#'   retention time axis.
#' @param mzExpWindow In case the \emph{m/z} window to plot an EIC for a particular analysis is not known (\emph{i.e.}
#'   no feature was detected of the feature group to be plot and \code{onlyPresent=FALSE}) then the EIC \emph{m/z} range
#'   is estimated from the range for the complete feature group and expanded by the offset defined by
#'   \code{mzExpWindow}.
#' @param topMost Only plot EICs from features within this number of top most intense analyses. If \code{NULL} then all
#'   analyses are used for plotted.
#' @param topMostByRGroup If set to \code{TRUE} and \code{topMost} is set: only plot EICs for the top most features in
#'   each replicate group. For instance, when \code{topMost=1} and \code{topMostByRGroup=TRUE}, then EICs will be
#'   plotted for the most intense feature of each replicate group.
#' @param EICs Internal parameter for now and should be kept at \code{NULL} (default).
#' @param showPeakArea Set to \code{TRUE} to display integrated chromatographic peak ranges by filling (shading) their
#'   areas.
#' @param showFGroupRect Set to \code{TRUE} to mark the full retention/intensity range of all features within a feature
#'   group by drawing a rectangle around it.
#' @param title Character string used for title of the plot. If \code{NULL} a title will be automatically generated.
#' @param onlyPresent If \code{TRUE} then EICs will only be generated for analyses in which a particular feature group
#'   was detected. Disabling this option might be useful to see if any features were 'missed'.
#' @param annotate If set to \code{"ret"} and/or \code{"mz"} then retention and/or \emph{m/z} values will be drawn for
#'   each plotted feature group.
#' @param showProgress if set to \code{TRUE} then a text progressbar will be displayed when all EICs are being plot. Set
#'   to \code{"none"} to disable any annotation.
#'
#' @template plot-lim
#'
#' @rdname feature-plotting
#' @export
setMethod("plotChroms", "featureGroups", function(obj, rtWindow = 30, mzExpWindow = 0.001, retMin = FALSE, topMost = NULL,
                                                  topMostByRGroup = FALSE, EICs = NULL, showPeakArea = FALSE,
                                                  showFGroupRect = TRUE, title = NULL,
                                                  colourBy = c("none", "rGroups", "fGroups"),
                                                  showLegend = TRUE, onlyPresent = TRUE,
                                                  annotate = c("none", "ret", "mz"), showProgress = FALSE,
                                                  xlim = NULL, ylim = NULL, ...)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertNumber, . ~ rtWindow + mzExpWindow, lower = 0, finite = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertFlag, . ~ retMin + topMostByRGroup + showPeakArea + showFGroupRect + showLegend +
               onlyPresent + showProgress,
           fixed = list(add = ac))
    checkmate::assertCount(topMost, positive = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertString(title, null.ok = TRUE, add = ac)
    
    colourBy <- checkmate::matchArg(colourBy, c("none", "rGroups", "fGroups"), add = ac)
    annotate <- checkmate::matchArg(annotate, c("none", "ret", "mz"), several.ok = TRUE, add = ac)
    
    assertXYLim(xlim, ylim, add = ac)
    
    checkmate::reportAssertions(ac)
    
    if (length(obj) == 0)
    {
        noDataPlot()
        return(invisible(NULL))
    }
    
    if (showLegend && colourBy == "none")
        showLegend <- FALSE
    
    fTable <- featureTable(obj)
    gTable <- groupTable(obj)
    gInfo <- groupInfo(obj)
    gCount <- nrow(gInfo)
    gNames <- names(obj)
    anaInfo <- analysisInfo(obj)
    ftind <- groupFeatIndex(obj)
    
    rGroups <- unique(anaInfo$group)
    
    if (is.null(EICs))
        EICs <- getEICsForFGroups(obj, rtWindow, mzExpWindow, topMost, topMostByRGroup, onlyPresent)
    EICFGroups <- unique(unlist(sapply(EICs, names)))
    
    if (colourBy == "rGroups")
    {
        EICColors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(rGroups))
        names(EICColors) <- rGroups
    }
    else if (colourBy == "fGroups")
    {
        EICColors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(gCount)
        names(EICColors) <- gNames
    }
    else
        EICColors <- "blue"
    
    fillColors <- adjustcolor(EICColors, alpha.f = 0.35)
    names(fillColors) <- names(EICColors)
    
    # get overall retention/intensity limits
    plotLimits <- list(rtRange = c(0, 0), maxInt = 0)
    
    for (grpi in seq_len(gCount))
    {
        rtrs <- unlist(sapply(seq_len(nrow(anaInfo)), function(anai)
        {
            ana <- anaInfo$analysis[anai]
            fti <- ftind[[grpi]][anai]
            if (fti == 0)
                return(NA)
            return(unlist(fTable[[ana]][fti, .(retmin, retmax)]))
        }))
        
        if (any(!is.na(rtrs)))
        {
            rtmin <- min(rtrs, na.rm = TRUE)
            rtmax <- max(rtrs, na.rm = TRUE)
            if (plotLimits$rtRange[1] == 0 || plotLimits$rtRange[1] > rtmin)
                plotLimits$rtRange[1] <- rtmin
            if (plotLimits$rtRange[2] < rtmax)
                plotLimits$rtRange[2] <- rtmax
        }
        
        plotLimits$maxInt <- max(plotLimits$maxInt, max(gTable[[grpi]]))
    }
    
    plotLimits$rtRange <- plotLimits$rtRange + c(-rtWindow, rtWindow)
    if (retMin)
        plotLimits$rtRange <- plotLimits$rtRange / 60
    
    if (is.null(title))
    {
        # NOTE: plotChroms() for sets override default
        if (gCount == 1)
            title <- sprintf("Group '%s'\nrt: %.1f - m/z: %.4f", names(gTable)[1],
                             if (retMin) gInfo[1, "rts"] / 60 else gInfo[1, "rts"],
                             gInfo[1, "mzs"])
        else
            title <- sprintf("%d feature groups", gCount)
    }
    
    anaIndsToPlot <- sapply(gNames, function(grp)
    {
        anasWGroup <- names(EICs)[sapply(EICs, function(e) !is.null(e[[grp]]))]
        ret <- match(anasWGroup, anaInfo$analysis)
        if (onlyPresent)
            ret <- ret[gTable[[grp]][ret] != 0]
        return(ret)
    }, simplify = FALSE)
    
    rGroupsInPlot <- unique(anaInfo$group[unlist(anaIndsToPlot)])
    fGroupsInPlot <- gNames[gNames %in% EICFGroups]
    
    oldp <- par(no.readonly = TRUE)
    if (showLegend)
    {
        makeLegend <- function(x, y, ...)
        {
            texts <- if (colourBy == "rGroups") rGroupsInPlot else fGroupsInPlot
            return(legend(x, y, texts, col = EICColors[texts],
                          text.col = EICColors[texts], lty = 1,
                          xpd = NA, ncol = 1, cex = 0.75, bty = "n", ...))
        }
        
        plot.new()
        leg <- makeLegend(0, 0, plot = FALSE)
        lw <- (grconvertX(leg$rect$w, to = "ndc") - grconvertX(0, to = "ndc"))
        lw <- min(lw, 0.5) # don't make it too wide
        par(omd = c(0, 1 - lw, 0, 1), new = TRUE)
    }
    
    if (is.null(xlim))
        xlim <- plotLimits$rtRange
    if (is.null(ylim))
        ylim <- c(0, plotLimits$maxInt * 1.1)
    
    plot(0, type = "n", main = title, xlab = sprintf("Retention time (%s)", if (retMin) "min." else "sec."), ylab = "Intensity",
         xlim = xlim, ylim = ylim, ...)
    
    if (showProgress)
        prog <- openProgBar(0, gCount)
    
    for (grp in fGroupsInPlot)
    {
        grpi <- match(grp, fGroupsInPlot)
        
        fts <- rbindlist(sapply(anaInfo$analysis, function(ana)
        {
            fti <- ftind[[grpi]][match(ana, anaInfo$analysis)]
            if (fti == 0)
                return(data.table(ret = NA, retmin = NA, retmax = NA, mz = NA))
            return(fTable[[ana]][fti])
        }, simplify = F), fill = TRUE)
        
        rtRange <- c(min(fts[anaIndsToPlot[[grp]], retmin], na.rm = TRUE), max(fts[anaIndsToPlot[[grp]], retmax], na.rm = TRUE))
        # rtEICRange <- rtRange + c(-rtWindow, rtWindow)
        
        if (retMin)
            rtRange <- rtRange / 60
        
        for (ana in names(EICs))
        {
            anai <- match(ana, anaInfo$analysis)
            if (is.na(anai))
                next
            
            if (onlyPresent && gTable[[grp]][anai] == 0)
                next
            
            EIC <- EICs[[ana]][[grp]]
            if (is.null(EIC))
                next
            
            if (colourBy == "rGroups")
                colInd <- anaInfo$group[anai]
            else if (colourBy == "fGroups")
                colInd <- grp
            else
                colInd <- 1
            
            points(if (retMin) EIC$time / 60 else EIC$time, EIC$intensity, type = "l", col = EICColors[colInd])
            
            if (showPeakArea && !is.na(fts[["mz"]][anai]))
            {
                EICFill <- EIC[numGTE(EIC$time, fts[anai, retmin]) & numLTE(EIC$time, fts[anai, retmax]), ]
                if (retMin)
                    EICFill$time <- EICFill$time / 60
                EICFill <- EICFill[EICFill$time %inrange% xlim, ]
                # filling doesn't work if outside y plot range
                EICFill$intensity[EICFill$intensity < ylim[1]] <- ylim[1]
                EICFill$intensity[EICFill$intensity > ylim[2]] <- ylim[2]
                polygon(c(EICFill$time, rev(EICFill$time)), c(EICFill$intensity, rep(0, length(EICFill$intensity))),
                        col = fillColors[colInd], border = NA)
            }
        }
        
        if (showFGroupRect || !"none" %in% annotate)
        {
            intRange <- c(0, max(fts[anaIndsToPlot[[grp]], intensity], na.rm = TRUE))
            
            if (showFGroupRect)
                rect(rtRange[1], intRange[1], rtRange[2], intRange[2], border = "red", lty = "dotted")
            
            if (!"none" %in% annotate)
            {
                antxt <- character()
                rt <- mean(fts[anaIndsToPlot[[grp]], ret], na.rm = TRUE)
                if (retMin)
                    rt <- rt / 60
                
                if ("ret" %in% annotate)
                    antxt <- sprintf("%.1f", rt)
                if ("mz" %in% annotate)
                    antxt <- paste(antxt, sprintf("%.4f", gInfo[grp, "mzs"]), sep = "\n")
                
                if (nzchar(antxt))
                    text(rt, intRange[2] + ylim[2] * 0.02, antxt)
            }
        }
        
        if (showProgress)
            setTxtProgressBar(prog, grpi)
    }
    
    if (showLegend)
        makeLegend(par("usr")[2], par("usr")[4])
    
    if (showProgress)
    {
        setTxtProgressBar(prog, gCount)
        close(prog)
    }
    
    par(oldp)
})

setMethod("plotChromsHash", "featureGroups", function(obj, rtWindow = 30, mzExpWindow = 0.001, retMin = FALSE,
                                                      topMost = NULL, topMostByRGroup = FALSE, EICs = NULL,
                                                      showPeakArea = FALSE, showFGroupRect = TRUE,
                                                      title = NULL, colourBy = c("none", "rGroups", "fGroups"),
                                                      showLegend = TRUE, onlyPresent = TRUE,
                                                      annotate = c("none", "ret", "mz"), showProgress = FALSE,
                                                      xlim = NULL, ylim = NULL, ...)
{
    colourBy <- checkmate::matchArg(colourBy, c("none", "rGroups", "fGroups"))
    annotate <- checkmate::matchArg(annotate, c("none", "ret", "mz"), several.ok = TRUE)
    if ("none" %in% annotate)
        annotate <- "none"
    args <- allArgs(FALSE)
    makeHash(args[setdiff(names(args), "obj")], featureTable(obj), groupInfo(obj), analysisInfo(obj))
})

#' @details \code{plotVenn} plots a Venn diagram (using \pkg{\link{VennDiagram}}) outlining unique and shared feature
#'   groups between up to five replicate groups.
#' @template plotvenn-ret
#' 
#' @rdname feature-plotting
#' @export
setMethod("plotVenn", "featureGroups", function(obj, which = NULL, ...)
{
    rGroups <- replicateGroups(obj)
    if (is.null(which))
        which <- rGroups
    
    checkmate::assert(checkmate::checkSubset(which, rGroups, empty.ok = FALSE),
                      checkmate::checkList(which, "character", any.missing = FALSE),
                      .var.name = "which")
    
    fGroupsList <- lapply(which, function(w) replicateGroupFilter(obj, w, verbose = FALSE))
    
    if (is.list(which))
    {
        if (!checkmate::testNamed(which))
            names(which) <- sapply(which, paste0, collapse = "+", USE.NAMES = FALSE)
    }
    else
        names(which) <- which
    
    makeVennPlot(lapply(fGroupsList, names), names(which), lengths(fGroupsList),
                 function(obj1, obj2) intersect(obj1, obj2),
                 length, ...)
})

#' @rdname feature-plotting
#' @export
setMethod("plotVenn", "featureGroupsSet", function(obj, which = NULL, ..., sets = FALSE)
{
    checkmate::assertFlag(sets)
    if (sets)
    {
        mySets <- get("sets", pos = 2)(obj)
        if (is.null(which))
            which <- mySets
        else
            checkmate::assertSubset(which, mySets)
        ai <- analysisInfo(obj)
        which = sapply(which, function(s) ai[ai$set == s, "group"], simplify = FALSE)
    }
    callNextMethod(obj, which = which, ...)
})

#' @details \code{plotUpSet} plots an UpSet diagram (using the \code{\link[UpSetR]{upset}} function) outlining unique
#'   and shared feature groups between given replicate groups.
#'   
#' @template plotUpSet
#'
#' @rdname feature-plotting
#' @export
setMethod("plotUpSet", "featureGroups", function(obj, which = NULL, nsets = length(which),
                                                 nintersects = NA, ...)
{
    rGroups <- replicateGroups(obj)
    if (is.null(which))
        which <- rGroups
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertSubset(which, rGroups, empty.ok = FALSE, add = ac)
    checkmate::assertCount(nsets, positive = TRUE)
    checkmate::assertCount(nintersects, positive = TRUE, na.ok = TRUE)
    checkmate::reportAssertions(ac)
    
    if (length(obj) == 0)
        stop("Can't plot empty feature groups object")
    
    obj <- replicateGroupFilter(obj, which, verbose = FALSE)
    
    gt <- as.data.table(obj, average = TRUE)
    gt <- gt[, which, with = FALSE] # isolate relevant columns
    gt[, (which) := lapply(.SD, function(x) as.integer(x > 0))]
    
    if (sum(sapply(gt[, which, with = FALSE], function(x) any(x>0))) < 2)
        stop("Need at least two replicate groups with non-zero intensities")
    
    UpSetR::upset(gt, nsets = nsets, nintersects = nintersects, ...)
})

#' @details \code{plotVolcano} Plots Fold change data in a 'Volcano plot'.
#' @param FCParams A parameter list to calculate Fold change data. See \code{getFCParams} for more details.
#' @param averageFunc,normalized Used for intensity data treatment, see the documentation for the
#'   \code{\link[=as.data.table,featureGroups-method]{as.data.table method}}.
#' @rdname feature-plotting
#' @export
setMethod("plotVolcano", "featureGroups", function(obj, FCParams, showLegend = TRUE, averageFunc = mean,
                                                   normalized = FALSE, col = NULL, pch = 19, ...)
{
    ac <- checkmate::makeAssertCollection()
    assertFCParams(FCParams, obj, null.ok = FALSE, add = ac)
    checkmate::assertFlag(showLegend, add = ac)
    checkmate::assertFunction(averageFunc, add = ac)
    checkmate::assertFlag(normalized, add = ac)
    checkmate::reportAssertions(ac)
    
    ac <- checkmate::makeAssertCollection()
    checkmate::reportAssertions(ac)
    
    if (length(obj) == 0)
        stop("Can't plot empty feature groups object")
    
    if (is.null(col))
        col <- getBrewerPal(5, "Paired")
    names(col) <- c("increase", "decrease", "FC", "significant", "insignificant")
    
    gt <- as.data.table(obj, FCParams = FCParams, averageFunc = averageFunc, normalized = normalized)
    gt[, colour := col[classification]]
    
    oldp <- par(no.readonly = TRUE)
    if (showLegend)
    {
        makeLegend <- function(x, y, ...)
        {
            return(legend(x, y, names(col), col = col, pch = pch, text.col = col, xpd = NA, ncol = 1,
                          cex = 0.75, bty = "n", ...))
        }
        
        plot.new()
        leg <- makeLegend(0, 0, plot = FALSE)
        lw <- (grconvertX(leg$rect$w, to = "ndc") - grconvertX(0, to = "ndc"))
        par(omd = c(0, 1 - lw, 0, 1), new = TRUE)
    }
    
    plot(gt$FC_log, gt$PV_log, xlab = "log2 fold change", ylab = "-log10 p-value", 
         col = gt$colour, pch = pch, ...)
    abline(v = c(-FCParams$thresholdFC, FCParams$thresholdFC), col = "red", lty = 2, lwd = 1, h = -log10(FCParams$thresholdPV))
    
    if (showLegend)
        makeLegend(par("usr")[2], par("usr")[4])
    
    par(oldp)
    
    invisible(NULL)
})

#' @details \code{plotGraph} generates an interactive network plot which is used to explore internal standard (IS)
#'   assignments to each feature group. This requires the availability of IS assignments, see the documentation for
#'   \code{\link{normInts}} for details. The graph is rendered with \pkg{\link{visNetwork}}.
#'
#' @param onlyPresent Only plot feature groups of internal standards that are still present in the \code{featureGroups}
#'   input object (which may be otherwise be removed by \emph{e.g.} subsetting or
#'   \code{\link[=filter,featureGroups-method]{filter}}).
#' 
#' @template plotGraph
#'
#' @rdname feature-plotting
#' @export
setMethod("plotGraph", "featureGroups", function(obj, onlyPresent = TRUE)
{
    checkmate::assertFlag(onlyPresent)
    
    if (length(obj) == 0 || nrow(internalStandards(obj)) == 0)
        stop("No feature groups to plot or no internal standards assigned.", call. = FALSE)
    
    ISTDs <- internalStandards(obj)
    ISTDAssign <- internalStandardAssignments(obj)
    gNames <- names(obj)
    
    if (onlyPresent)
        ISTDAssign <- pruneList(lapply(ISTDAssign, function(ia) ia[ia %chin% gNames]), checkEmptyElements = TRUE)
    
    nodes <- data.table(id = union(names(ISTDAssign), unlist(ISTDAssign)))
    nodes[, group := fifelse(id %chin% names(ISTDAssign), "fGroup", "ISTD")]
    nodes[group == "ISTD", ISTD := paste0(ISTDs[group == id]$name, collapse = ","), by = "id"]
    nodes[group == "fGroup", ISTD := paste0(ISTDs[group %in% ISTDAssign[[id]]]$name, collapse = ","), by = id]
    nodes[, label := id]
    
    gInfo <- groupInfo(obj)
    sInfo <- if (isScreening(obj)) screenInfo(obj) else NULL
    nodes[id %chin% gNames, title := mapply(id, group, FUN = function(grp, type)
    {
        istds <- if (type == "ISTD")
            getStrListWithMax(ISTDs[group == grp]$name, 6, "/")
        else
            getStrListWithMax(ISTDs[group %chin% ISTDAssign[[grp]]]$name, 3, "/")
        ret <- sprintf("<b>%s</b><br>RT: %.2f<br>m/z: %.4f<br>ISTD: %s", grp, gInfo[grp, "rts"], gInfo[grp, "mzs"],
                       istds)
        if (!is.null(sInfo) && grp %chin% sInfo$group)
            ret <- paste0(ret, "<br>", "Suspect(s): ", getStrListWithMax(sInfo[group == grp]$name, 3, "/"))
        return(ret)
    })]
    nodes[is.na(title), title := sprintf("<b>%s</b> (removed)", id)]
    
    edges <- rbindlist(Map(names(ISTDAssign), ISTDAssign, f = function(grp, ia) data.table(from = ia, to = grp)))
    
    # based on default defined in visInteraction() --> decreased font-size
    titleStyle <- paste("position: fixed; visibility:hidden; padding: 5px; white-space: nowrap; font-family: verdana;",
                        "font-size:10px; font-color:#000000; background-color: #f5f4ed; -moz-border-radius: 3px;",
                        "-webkit-border-radius: 3px; border-radius: 3px; border: 1px solid #808074;",
                        "box-shadow: 3px 3px 10px rgba(0, 0, 0, 0.2);")
    
    edges[, value := abs(gInfo[from, "rts"] - gInfo[to, "rts"])]
    edges[, value := max(0.1, 1 - (value / max(value)))]
    
    visNetwork::visNetwork(nodes, edges,
                           submain = paste0("Explore connections by dragging/zooming/selecting.<br>",
                                            "Smaller retention time difference have wider edges.")) %>%
        visNetwork::visOptions(selectedBy = list(variable = "ISTD", multiple = TRUE),
                               highlightNearest = list(enabled = TRUE, hover = TRUE, algorithm = "hierarchical"),
                               nodesIdSelection = list(enabled = TRUE, main = "Select by feat group")) %>%
        visNetwork::visIgraphLayout(layout = "layout_with_lgl") %>%
        visNetwork::visEdges(arrows = "from", scaling = list(min = 0.5, max = 2)) %>%
        visNetwork::visInteraction(tooltipStyle = titleStyle, hideEdgesOnDrag = TRUE, hideEdgesOnZoom = TRUE) %>%
        visNetwork::visLegend()
})

#' @rdname feature-plotting
#' @export
setMethod("plotGraph", "featureGroupsSet", function(obj, onlyPresent = TRUE, set) plotGraph(unset(obj, set), onlyPresent = onlyPresent))
