#' @include main.R
NULL

#' Plotting of grouped features
#'
#' Various plotting functions for feature group data.
#'
#' @param obj,x \code{featureGroups} object to be used for plotting.
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
#' @param \dots passed to \code{\link[base]{plot}} (\code{plot} and \code{plotChroms}), \pkg{\link{VennDiagram}}
#'   plotting functions (\code{plotVenn}), \code{\link{chordDiagram}} (\code{plotChord}) or \code{\link[UpSetR]{upset}}
#'   (\code{plotUpSet}).
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
#'   \item \code{plotGraph} only plots data per set, and requires the \code{set} argument to be set.
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

setMethod("plotHash", "featureGroups", function(x, ...)
{
    return(makeHash(groupTable(x), ...))
})

#' @details \code{plotInt} Generates a line plot for the (averaged) intensity of feature groups within all analyses
#' @param xnames Plot analysis (or replicate group if \code{average=TRUE}) names on the x axis.
#' @param plotArgs,linesArgs A \code{list} with further arguments passed to \code{\link[base]{plot}} and
#'    \code{\link[graphics]{lines}}, respectively.
#'
#' @rdname feature-plotting
#' @export
setMethod("plotInt", "featureGroups", function(obj, average = FALSE, normalized = FALSE, xnames = TRUE,
                                               showLegend = FALSE, pch = 20, type = "b", lty = 3, col = NULL,
                                               plotArgs = NULL, linesArgs = NULL)
{
    aapply(checkmate::assertFlag, . ~ average + normalized + xnames + showLegend)
    aapply(checkmate::assertList, . ~ plotArgs + linesArgs, null.ok = TRUE)
    doPlotFeatInts(obj, average, normalized, xnames, showLegend, pch, type, lty, col, plotArgs, linesArgs,
                   doSets = FALSE)    
})

#' @rdname feature-plotting
#' @export
setMethod("plotInt", "featureGroupsSet", function(obj, average = FALSE, normalized = FALSE, xnames = !sets,
                                                  showLegend = sets, pch = 20, type = "b", lty = 3, col = NULL,
                                                  plotArgs = NULL, linesArgs = NULL, sets = FALSE)
{
    aapply(checkmate::assertFlag, . ~ average + normalized + xnames + showLegend + sets)
    aapply(checkmate::assertList, . ~ plotArgs + linesArgs, null.ok = TRUE)
    doPlotFeatInts(obj, average, normalized, xnames, showLegend, pch, type, lty, col, plotArgs, linesArgs,
                   doSets = sets)    
})

setMethod("plotIntHash", "featureGroups", function(obj, average = FALSE, ...) makeHash(groupTable(obj), average, ...))

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
    
    groupTab <- as.data.table(obj, average = average)
    
    getLinkScore <- function(sn1, sn2)
    {
        if (sn1 == sn2)
            return(0)
        return(sum(groupTab[[sn1]] > 0 & groupTab[[sn2]] > 0))
    }
    
    chordTable[, value := mapply(getLinkScore, from, to)]
    
    if (addSelfLinks)
    {
        gt <- if (average) averageGroups(obj) else groupTable(obj)
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
            ftgrps <- groupTab[get(sn) > 0]$group
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

setMethod("plotChordHash", "featureGroups", function(obj, ...)
{
    return(makeHash(groupTable(obj), ...))
})

#' @details \code{plotChroms} Plots extracted ion chromatograms (EICs) of feature groups.
#'
#' @param analysis,groupName \code{character} vector with the analyses/group names to be considered for plotting.
#'   Compared to subsetting the \code{featureGroups} object (\code{obj}) upfront this is slightly faster and (if
#'   \code{onlyPresent=FALSE}) allows plotting chromatograms for feature groups where none of the specified analyses
#'   contain the feature (which is impossible otherwise since subsetting leads to removal of 'empty' feature groups).
#' @param EICs Internal parameter for now and should be kept at \code{NULL} (default).
#' @param showPeakArea Set to \code{TRUE} to display integrated chromatographic peak ranges by filling (shading) their
#'   areas.
#' @param showFGroupRect Set to \code{TRUE} to mark the full retention/intensity range of all features within a feature
#'   group by drawing a rectangle around it.
#' @param title Character string used for title of the plot. If \code{NULL} a title will be automatically generated.
#' @param annotate If set to \code{"ret"} and/or \code{"mz"} then retention and/or \emph{m/z} values will be drawn for
#'   each plotted feature group.
#' @param intMax Method used to determine the maximum intensity plot limit. Should be \code{"eic"} (from EIC data) or
#'   \code{"feature"} (from feature data). Ignored if the \code{ylim} parameter is specified.
#' @param showProgress if set to \code{TRUE} then a text progressbar will be displayed when all EICs are being plot. Set
#'   to \code{"none"} to disable any annotation.
#'
#' @template EICParams-arg
#' @template plot-lim
#'
#' @rdname feature-plotting
#' @export
setMethod("plotChroms", "featureGroups", function(obj, analysis = analyses(obj), groupName = names(obj),
                                                  retMin = FALSE, showPeakArea = FALSE, showFGroupRect = TRUE,
                                                  title = NULL, colourBy = c("none", "rGroups", "fGroups"),
                                                  showLegend = TRUE, annotate = c("none", "ret", "mz"),
                                                  intMax = "eic", EICParams = getDefEICParams(), showProgress = FALSE,
                                                  xlim = NULL, ylim = NULL, EICs = NULL, ...)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertSubset, . ~ analysis + groupName, list(analyses(obj), names(obj)), empty.ok = TRUE,
           fixed = list(add = ac))
    aapply(checkmate::assertFlag, . ~ retMin + showPeakArea + showFGroupRect + showLegend + showProgress,
           fixed = list(add = ac))
    checkmate::assertString(title, null.ok = TRUE, add = ac)
    colourBy <- checkmate::matchArg(colourBy, c("none", "rGroups", "fGroups"), add = ac)
    annotate <- checkmate::matchArg(annotate, c("none", "ret", "mz"), several.ok = TRUE, add = ac)
    checkmate::assertChoice(intMax, c("eic", "feature"), add = ac)
    assertEICParams(EICParams, add = ac)
    assertXYLim(xlim, ylim, add = ac)
    checkmate::reportAssertions(ac)

    if (intMax == "feature" && !EICParams$onlyPresent)
        stop("intMax must be 'eic' when EICParams$onlyPresent == FALSE", call. = FALSE)
    
    if (showLegend && colourBy == "none")
        showLegend <- FALSE
    
    if (is.null(EICs))
        EICs <- getEICsForFGroups(obj, analysis, groupName, EICParams)
    else
    {
        # sync as much as possible with given EICParams
        EICs <- filterEICs(EICs, obj, analysis = analysis, groupName = groupName, topMost = EICParams$topMost,
                           topMostByRGroup = EICParams$topMostByRGroup, onlyPresent = EICParams$onlyPresent)
    }
    
    if (length(obj) == 0 || length(EICs) == 0)
    {
        noDataPlot()
        return(invisible(NULL))
    }
    
    gInfo <- groupInfo(obj)
    gCount <- length(groupName)
    anaInfo <- analysisInfo(obj)
    takeAnalysis <- analysis
    anaInfo <- anaInfo[analysis %chin% takeAnalysis & analysis %chin% names(EICs)]
    featTab <- as.data.table(getFeatures(obj))[group %chin% groupName]
    rGroups <- unique(anaInfo$group)
    
    if (colourBy == "rGroups")
    {
        EICColors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(rGroups))
        names(EICColors) <- rGroups
    }
    else if (colourBy == "fGroups")
    {
        EICColors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(gCount)
        names(EICColors) <- groupName
    }
    else
        EICColors <- "blue"
    
    fillColors <- adjustcolor(EICColors, alpha.f = 0.35)
    names(fillColors) <- names(EICColors)
    
    if (is.null(xlim))
    {
        xlim <- c(min(featTab$retmin) - EICParams$rtWindow, max(featTab$retmax) + EICParams$rtWindow)
        if (retMin)
            xlim <- xlim / 60
    }
    if (is.null(ylim))
    {
        if (intMax == "eic")
        {
            RTRangeGrp <- split(featTab[, .(retmin = min(retmin), retmax = max(retmax)), by = "group"], by = "group", keep.by = FALSE)
            plotIntMax <- max(unlist(lapply(EICs, function(aeic) Map(names(aeic), aeic, f = function(grp, eic)
            {
                rtr <- unlist(RTRangeGrp[[grp]])
                if (!is.null(rtr))
                {
                    if (!is.null(xlim))
                    {
                        xl <- if (retMin) xlim * 60 else xlim
                        rtr <- c(max(rtr[1], xl[1]), min(rtr[2], xl[2]))
                    }
                    eici <- eic[eic$time %between% rtr, "intensity"]
                    if (length(eici) > 0)
                        return(max(eici))
                }
                return(0)
            }))))
        }
        else # "feature"
            plotIntMax <- max(featTab$intensity)
        ylim <- c(0, plotIntMax * 1.1)
    }

    if (is.null(title))
    {
        # NOTE: plotChroms() for sets override default
        if (gCount == 1)
            title <- sprintf("Group '%s'\nrt: %.1f - m/z: %.4f", groupName[1],
                             if (retMin) gInfo[groupName[1], "rts"] / 60 else gInfo[groupName[1], "rts"],
                             gInfo[groupName[1], "mzs"])
        else
            title <- sprintf("%d feature groups", gCount)
    }
    
    if (showLegend)
    {
        makeLegend <- function(x, y, ...)
        {
            texts <- if (colourBy == "rGroups") rGroups else groupName
            return(legend(x, y, texts, col = EICColors[texts],
                          text.col = EICColors[texts], lty = 1,
                          xpd = NA, ncol = 1, cex = 0.75, bty = "n", ...))
        }
        
        plot.new()
        leg <- makeLegend(0, 0, plot = FALSE)
        lw <- (grconvertX(leg$rect$w, to = "ndc") - grconvertX(0, to = "ndc"))
        lw <- min(lw, 0.5) # don't make it too wide
        withr::local_par(list(omd = c(0, 1 - lw, 0, 1), new = TRUE))
    }
    
    plot(0, type = "n", main = title, xlab = sprintf("Retention time (%s)", if (retMin) "min." else "sec."),
         ylab = "Intensity", xlim = xlim, ylim = ylim, ...)
    
    effectiveXlim <- par("usr")[c(1, 2)]
    effectiveYlim <- par("usr")[c(3, 4)]
    
    if (showProgress)
        prog <- openProgBar(0, gCount)
    
    for (grp in groupName)
    {
        featTabGrp <- featTab[group == grp]
        
        for (ana in analysis)
        {
            EIC <- EICs[[ana]][[grp]]
            if (is.null(EIC))
                next
            
            featRow <- featTabGrp[analysis == ana]
            if (EICParams$onlyPresent && nrow(featRow) == 0)
                next
            
            if (colourBy == "rGroups")
                colInd <- anaInfo[match(ana, analysis)]$group
            else if (colourBy == "fGroups")
                colInd <- grp
            else
                colInd <- 1
            
            points(if (retMin) EIC$time / 60 else EIC$time, EIC$intensity, type = "l", col = EICColors[colInd])
            
            if (showPeakArea && nrow(featRow) != 0)
            {
                EICFill <- setDT(EIC[numGTE(EIC$time, featRow$retmin) & numLTE(EIC$time, featRow$retmax), ])
                if (retMin)
                    EICFill[, time := time / 60]
                EICFill <- EICFill[time %inrange% effectiveXlim]
                # filling doesn't work if outside y plot range
                EICFill[intensity < effectiveYlim[1], intensity := effectiveYlim[1]]
                EICFill[intensity > effectiveYlim[2], intensity := effectiveYlim[2]]
                polygon(c(EICFill$time, rev(EICFill$time)), c(EICFill$intensity, rep(0, length(EICFill$intensity))),
                        col = fillColors[colInd], border = NA)
            }
        }
        
        if (showFGroupRect || !"none" %in% annotate)
        {
            rtRange <- c(min(featTabGrp$retmin), max(featTabGrp$retmax))
            if (retMin)
                rtRange <- rtRange / 60
            maxInt <- max(featTabGrp$intensity)
            
            if (showFGroupRect)
                rect(rtRange[1], 0, rtRange[2], maxInt, border = "red", lty = "dotted")
            
            if (!"none" %in% annotate)
            {
                antxt <- character()
                rt <- mean(featTabGrp$ret)
                if (retMin)
                    rt <- rt / 60
                
                if ("ret" %in% annotate)
                    antxt <- sprintf("%.1f", rt)
                if ("mz" %in% annotate)
                    antxt <- paste(antxt, sprintf("%.4f", gInfo[grp, "mzs"]), sep = "\n")
                
                if (nzchar(antxt))
                    text(rt, maxInt + ylim[2] * 0.02, antxt)
            }
        }
        
        if (showProgress)
            setTxtProgressBar(prog, match(grp, groupName))
    }
    
    if (showLegend)
        makeLegend(par("usr")[2], par("usr")[4])
    
    if (showProgress)
    {
        setTxtProgressBar(prog, gCount)
        close(prog)
    }
})

setMethod("plotChromsHash", "featureGroups", function(obj, analysis = analyses(obj), groupName = names(obj),
                                                      retMin = FALSE, showPeakArea = FALSE, showFGroupRect = TRUE,
                                                      title = NULL, colourBy = c("none", "rGroups", "fGroups"),
                                                      showLegend = TRUE, annotate = c("none", "ret", "mz"),
                                                      intMax = "eic", EICParams = getDefEICParams(),
                                                      showProgress = FALSE, xlim = NULL, ylim = NULL, EICs = NULL, ...)
{
    colourBy <- checkmate::matchArg(colourBy, c("none", "rGroups", "fGroups"))
    annotate <- checkmate::matchArg(annotate, c("none", "ret", "mz"), several.ok = TRUE)
    if ("none" %in% annotate)
        annotate <- "none"
    args <- allArgs(FALSE)
    if (!is.null(EICs))
    {
        # omit data we don't need: speeds up hashing quite a bit
        # NOTE: only apply analysis/group filters, as the rest will slow down things considerably. Hence, this could
        # result in cache misses.
        EICs <- filterEICs(EICs, obj, analysis = analysis, groupName = groupName, topMost = NULL,
                           topMostByRGroup = FALSE, onlyPresent = FALSE)
    }
    anas <- analysis
    makeHash(args[setdiff(names(args), c("obj", "EICs"))], EICs, featureTable(obj)[analysis], groupInfo(obj)[groupName, ],
             analysisInfo(obj)[analysis %chin% anas])
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
        which = sapply(which, function(s) ai[set == s]$group, simplify = FALSE)
    }
    callNextMethod(obj, which = which, ...)
})

setMethod("plotVennHash", "featureGroups", function(obj, ...)
{
    return(makeHash(groupTable(obj), ...))
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

setMethod("plotUpSetHash", "featureGroups", function(obj, ...)
{
    return(makeHash(groupTable(obj), ...))
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
setMethod("plotGraph", "featureGroups", function(obj, onlyPresent = TRUE, width = NULL, height = NULL)
{
    checkmate::assertFlag(onlyPresent)
    
    ISTDs <- internalStandards(obj)
    ISTDAssign <- internalStandardAssignments(obj)
    gNames <- names(obj)
    
    if (onlyPresent)
        ISTDAssign <- pruneList(lapply(ISTDAssign, function(ia) ia[ia %chin% gNames]), checkEmptyElements = TRUE)
    
    nodes <- data.table(id = character(), label = character(), group = character())
    edges <- data.table()
    if (length(ISTDAssign) > 0)
    {
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
    }
    
    # based on default defined in visInteraction() --> decreased font-size
    titleStyle <- paste("position: fixed; visibility:hidden; padding: 5px; white-space: nowrap; font-family: verdana;",
                        "font-size:10px; font-color:#000000; background-color: #f5f4ed; -moz-border-radius: 3px;",
                        "-webkit-border-radius: 3px; border-radius: 3px; border: 1px solid #808074;",
                        "box-shadow: 3px 3px 10px rgba(0, 0, 0, 0.2);")
    
    if (nrow(edges) > 0)
    {
        edges[, value := abs(gInfo[from, "rts"] - gInfo[to, "rts"])]
        edges[, value := max(0.1, 1 - (value / max(value)))]
    }
    else
        edges[, c("from", "to") := character()]

    gr <- visNetwork::visNetwork(nodes, edges, width = width, height = height,
                                 submain = paste0("Explore connections by dragging/zooming/selecting.<br>",
                                                  "Smaller retention time difference have wider edges."))
    if (nrow(edges) > 0)
    {
        gr <- gr %>%
            visNetwork::visOptions(selectedBy = list(variable = "ISTD", multiple = TRUE),
                                   highlightNearest = list(enabled = TRUE, hover = TRUE, algorithm = "hierarchical"),
                                   nodesIdSelection = list(enabled = TRUE, main = "Select by feat group")) %>%
            visNetwork::visIgraphLayout(layout = "layout_with_lgl") %>%
            visNetwork::visEdges(arrows = "from", scaling = list(min = 0.5, max = 2)) %>%
            visNetwork::visInteraction(tooltipStyle = titleStyle, hideEdgesOnDrag = TRUE, hideEdgesOnZoom = TRUE) %>%
            visNetwork::visLegend()
    }
    
    return(gr)
})

#' @rdname feature-plotting
#' @param set \setsWF The set for which data must be plotted.
#' @export
setMethod("plotGraph", "featureGroupsSet", function(obj, onlyPresent = TRUE, set, ...) plotGraph(unset(obj, set), onlyPresent = onlyPresent, ...))
