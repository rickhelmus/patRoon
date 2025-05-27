# SPDX-FileCopyrightText: 2016 - 2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
NULL

#' Plotting of grouped features
#'
#' Various plotting functions for feature group data.
#'
#' @param obj,x \code{featureGroups} object to be used for plotting.
#' @param retMin Plot retention time in minutes (instead of seconds).
#' @param pch,type,lty Common plotting parameters passed to \emph{e.g.} \code{\link[graphics]{plot}}.
#'
#'   For \code{plot}: if \code{pch=NULL} then values are automatically assigned.
#'
#'   for \code{plotInt}: the \code{type} argument defaults to single points (\code{"p"}) when \code{regression=TRUE},
#'   and lines+points (\code{"b"}) otherwise. For non-aggregated data, it probably makes more sense to set \emph{e.g.}
#'   \code{type="p"} to avoid lines being drawn to points with equal x-values.
#' @param col Colour(s) used. If \code{col=NULL} then colours are automatically generated.
#' @param groupBy Specifies how results are grouped in the plot. Should be a name of a column in the
#'   \link[=analysis-information]{analysis information} table which is used to make analysis groups (\emph{e.g.}
#'   \code{"replicate"}). Set to \code{NULL} for no grouping. The \code{groupBy} argument can also be set to
#'   \code{"fGroups"} to group by feature groups (except \code{plotChord}).
#'
#'   For \code{plotInt}: see the \verb{Data aggregation in intensity plots} section for more details.
#'
#'   For \code{plotChord}: the grouping is used to generate 'outer groups'.
#' @param showLegend Plot a legend if \code{TRUE}.
#' @param which A character vector with the selection to compare (\emph{e.g.} replicates, as set by the \code{aggregate}
#'   argument). Set to \code{NULL} to select everything.
#' @param averageFunc,normalized Used for intensity data treatment, see the documentation for the
#'   \code{\link[=as.data.table,featureGroups-method]{as.data.table method}}.
#' @param analysis,groupName \code{character} vector with the analyses/group names to be considered for plotting.
#'   Compared to subsetting the \code{featureGroups} object (\code{obj}) upfront this is slightly faster. Furthermore,
#'   if \code{onlyPresent=FALSE} in \code{EICParams} or \code{EIMParams}, this allows plotting chromatograms for feature
#'   groups where none of the specified analyses contain the feature (which is impossible otherwise since subsetting
#'   leads to removal of 'empty' feature groups).
#'
#'   For \code{plotChroms} and IMS workflows: if \code{IMS!="both"} then the \code{analysis} and \code{groupName}
#'   arguments are adjusted for the remaining data after IMS selection.
#' @param EICs,EIMs Internal parameter for now and should be kept at \code{NULL} (default).
#' @param showPeakArea Set to \code{TRUE} to display integrated peak ranges by filling (shading) their areas.
#' @param showFGroupRect Set to \code{TRUE} to mark the full integration/intensity range of all features within a
#'   feature group by drawing a rectangle around it.
#' @param title Character string used for title of the plot. If \code{NULL} a title will be automatically generated.
#' @param annotate Set to \code{"ret"}, \code{"mz"} and/or \code{"mob"} to annotate peaks with the retention time,
#'   \emph{m/z} and/or ion mobility (if available) feature group values, respectively.
#' @param showProgress if set to \code{TRUE} then a text progressbar will be displayed when all EICs are being plot. Set
#'   to \code{"none"} to disable any annotation.
#' @param \dots passed to \code{\link[base]{plot}} (\code{plot}, \code{plotChroms}, \code{plotTICs} and
#'   \code{plotBPCs}), \pkg{\link{VennDiagram}} plotting functions (\code{plotVenn}), \code{\link{chordDiagram}}
#'   (\code{plotChord}) or \code{\link[UpSetR]{upset}} (\code{plotUpSet}).
#'
#' @templateVar consider for plotting
#' @template IMS-arg
#'
#' @template plot-lim
#'
#' @inheritParams featureGroups-class
#'
#' @templateVar what \code{plotChroms}, \code{plotMobilograms}, \code{plotTICs} and \code{plotBPCs}
#' @template uses-msdata
#'
#' @section Sets workflows: \setsWFChangedMethods{
#'
#'   \item \code{plotGraph} only plots data per set, and requires the \code{set} argument to be set.
#'
#'   }
#'
#'   In sets workflows the \link[=analysis-information]{analysis information} contains an additional \code{"set"}
#'   column, which can be used for arguments that involve grouping of analyses. For instance, if \code{groupBy="set"}
#'   then plotting data is grouped per set.
#'
#' @author Rick Helmus <\email{r.helmus@@uva.nl}> and Ricardo Cunha <\email{cunha@@iuta.de}> (\code{plotTICs} and
#'   \code{plotBPCs} functions)
#'
#' @seealso \code{\link{featureGroups-class}}, \code{\link{groupFeatures}}
#'
#' @name feature-plotting
NULL

#' @details \code{plot} Generates an \emph{m/z} \emph{vs} retention time plot for all featue groups. Optionally
#'   highlights unique/overlapping presence amongst replicates.
#' @param onlyUnique If \code{TRUE} and \code{groupBy} is set to a column name of the
#'   \link[=analysis-information]{analysis information} then only feature groups that are unique to a replicate are
#'   plotted.
#'
#' @rdname feature-plotting
#' @export
setMethod("plot", c(x = "featureGroups", y = "missing"), function(x, groupBy = NULL, onlyUnique = FALSE,
                                                                  retMin = FALSE, showLegend = TRUE, IMS = "maybe",
                                                                  col = NULL, pch = NULL, ...)
{
    assertIMSArg(IMS)
    x <- prepIMSFGroupsForPlot(x, IMS)
    
    anaInfo <- analysisInfo(x)

    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertFlag, . ~ onlyUnique + retMin + showLegend, fixed = list(add = ac))
    assertAnaInfoBy(groupBy, anaInfo, TRUE, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (length(x) == 0)
        noDataPlot()
    else
    {
        if (is.null(groupBy))
        {
            if (is.null(col))
                col <- "black"
            if (is.null(pch))
                pch <- 16
            showLegend <- FALSE
        }
        else if (groupBy == "fGroups")
        {
            if (is.null(col))
                col <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(x))
            if (is.null(pch))
                pch <- 16
            
            if (showLegend)
            {
                labels <- names(x)
                labCol <- rep(col, length.out = length(labels))
                labPch <- rep(pch, length.out = length(labels))
                names(labCol) <- labels; names(labPch) <- labels
            }
        }
        else
        {
            allGroups <- unique(anaInfo[[groupBy]])
            
            labels <- c(allGroups, "overlap")
            
            if (is.null(col))
                labCol <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(labels))
            else
                labCol <- rep(col, length.out = length(labels))
            names(labCol) <- labels
            
            if (is.null(pch))
            {
                # prefer closed symbols (15-25)
                ll <- length(labels)
                if (ll < (25 - 15))
                    labPch <- seq_len(ll) + 14
                else if (ll <= 25)
                    labPch <- seq_len(ll)
                else
                    labPch <- rep(16, ll) # just stick with one
            }
            else
                labPch <- rep(pch, length.out = length(labels))
            names(labPch) <- labels
            
            # get averaged intensities and omit initial name/rt/mz columns
            gTable <- as.data.table(x, average = groupBy)[, getADTIntCols(allGroups), with = FALSE]
            
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
        
        plot(if (retMin) x@groupInfo$ret / 60 else x@groupInfo$ret, x@groupInfo$mz,
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

#' @details \code{plotInt} Generates a line plot with feature intensities.
#'
#' @param average Controls plot data averaging: see the \verb{Data aggregation in intensity plots} section.
#' @param xBy Controls x-value grouping in the plot: see the \verb{Data aggregation in intensity plots} section.
#' @param areas Set to \code{TRUE} to use feature areas instead of peak intensities for plotting.
#' @param xNames Plot names (\code{xNames=TRUE}) or a number sequence (\code{xNames=FALSE}) on the x axis.
#' @param regression If \code{TRUE} then a regression line is plotted for each group, using the x values set by
#'   \code{xBy} (see \verb{Data aggregation in intensity plots}). if \code{showLegend=TRUE} then the R-squared value is
#'   show in the legend for each group (unless there are multiple feature groups and \code{groupBy != "fGroups"}).
#' @param plotArgs,linesArgs A \code{list} with further arguments passed to \code{\link[base]{plot}} and
#'   \code{\link[graphics]{lines}}, respectively.
#'
#' @section Data aggregation in intensity plots: the \code{average}, \code{xBy} and \code{groupBy} arguments control how
#'   data is aggregated in intensity plots: \itemize{
#'
#'   \item \code{average}: controls the averaging of feature intensities prior to plotting.
#'
#'   \item \code{xBy}: can map the x value of individual points to analysis metadata. For example, exposure time or
#'   sample location. Non-numeric values are allowed (unless \code{regression=TRUE}).
#'
#'   \item \code{groupBy}: controls the grouping of points in the plot. Equal groups are plotted in sequence so they can
#'   be connected with lines and are coloured equally. Examples include experiment type or feature groups.
#'
#'   }
#'
#'   The following values are valid: \itemize{
#'
#'   \item \code{FALSE} (\code{average}) or \code{NULL} (\code{xBy} and \code{groupBy}): aggregation is disabled.
#'
#'   \item \code{TRUE} (only \code{average}): results are averaged for each replicate
#'
#'   \item a name of a column in the \link[=analysis-information]{analysis information}: results are aggregated for
#'   analyses with the same table column value.
#'
#'   \item \code{"fGroups"} (only \code{groupBy}): plots are grouped by feature groups.
#'
#'   }
#'
#' @rdname feature-plotting
#' @export
setMethod("plotInt", "featureGroups", function(obj, average = FALSE, averageFunc = mean, areas = FALSE,
                                               normalized = FALSE, xBy = NULL, xNames = TRUE, groupBy = "fGroups",
                                               regression = FALSE, showLegend = FALSE, IMS = "maybe", pch = 20,
                                               type = if (regression) "p" else "b", lty = 3,
                                               col = NULL, plotArgs = NULL, linesArgs = NULL)
{
    assertIMSArg(IMS)
    obj <- prepIMSFGroupsForPlot(obj, IMS)
    
    anaInfo <- analysisInfo(obj)
    
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertFlag, . ~ areas + normalized + xNames + regression + showLegend,
           fixed = list(add = ac))
    checkmate::assertFunction(averageFunc, add = ac)
    assertAnaInfoBy(xBy, anaInfo, FALSE, null.ok = TRUE, add = ac)
    assertAnaInfoBy(groupBy, anaInfo, TRUE, null.ok = TRUE, add = ac)
    aapply(checkmate::assertList, . ~ plotArgs + linesArgs, null.ok = TRUE, fixed = list(add = ac))
    checkmate::reportAssertions(ac)

    averageBy <- assertAndPrepareAnaInfoBy(average, anaInfo, FALSE)
    
    if (is.null(xBy))
        xBy <- averageBy
    else
        checkAnaInfoAggrGrouping(anaInfo, "averaged", averageBy, xBy)
    
    if (is.null(groupBy))
    {
        if (showLegend)
            showLegend <- FALSE
    }
    else
        checkAnaInfoAggrGrouping(anaInfo, "averaged", averageBy, groupBy)
    
    if (length(obj) == 0)
    {
        noDataPlot()
        return(invisible(NULL))
    }

    if (normalized)
        obj <- maybeAutoNormalizeFGroups(obj)
    
    intTab <- if (isFALSE(average))
        copy(groupTable(obj, areas, normalized))
    else
        averageGroups(obj, areas, normalized, averageBy, averageFunc)
    
    intTab[, avgGroup := unique(anaInfo[[averageBy]])]
    intTab[, x := anaInfo[[xBy]][match(avgGroup, anaInfo[[averageBy]])]]
    intTab[, xnum := match(x, unique(anaInfo[[xBy]]))]
    if (!is.null(groupBy) && groupBy != "fGroups")
        intTab[, xgroup := anaInfo[[groupBy]][match(avgGroup, anaInfo[[averageBy]])]]
    
    if (is.null(groupBy))
    {
        if (is.null(col))
            col <- "black"
    }
    else
    {
        labs <- if (groupBy == "fGroups") names(obj) else unique(intTab$xgroup)
        if (is.null(col))
            col <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(labs))
        else
            col <- rep(col, length.out = length(labs))
        names(col) <- labs
    }

    xNum <- is.numeric(intTab$x)
    
    if (regression && !xNum)
    {
        stop("The data in the analysis information column specified by xBy should be numeric for regression calculations",
             call. = FALSE)
    }

    regList <- if (regression)
    {
        sapply(intTab[, names(obj), with = FALSE], function(ints)
        {
            if (is.null(groupBy) || groupBy == "fGroups")
                calcFeatureRegression(intTab$x, ints)
            else
            {
                sapply(unique(intTab$xgroup), function(xg)
                {
                    wh <- intTab[, .I[xgroup == xg]]
                    calcFeatureRegression(intTab$x[wh], ints[wh])
                }, simplify = FALSE)
            }
        }, simplify = FALSE)
    }
    else
        NULL

    oldp <- par(no.readonly = TRUE)
    if (showLegend)
    {
        makeLegend <- function(x, y, ...)
        {
            leg <- names(col)
            if (regression && (groupBy == "fGroups" || length(regList) == 1))
            {
                RSQs <- if (groupBy == "fGroups") sapply(regList, "[[", "RSQ") else sapply(regList[[1]], "[[", "RSQ")
                RSQs <- sprintf("%.2f", RSQs)
                leg <- mapply(leg, RSQs, FUN = function(l, r) as.expression(bquote(.(l) ~ "(" ~ R^2 ~ .(r) ~ ")")))
            }
            return(legend(x, y, leg, col = col, pch = pch, text.col = col, xpd = NA, ncol = 1,
                          cex = 0.75, bty = "n", ...))
        }
        
        plot.new()
        leg <- makeLegend(0, 0, plot = FALSE)
        lw <- (grconvertX(leg$rect$w, to = "ndc") - grconvertX(0, to = "ndc"))
        par(omd = c(0, 1 - lw, 0, 1), new = TRUE)
    }
    
    maxX <- if (xNum) max(intTab$x) else uniqueN(intTab$x)
    do.call(plot, c(list(x = NULL, xlim = c(0, maxX), ylim = c(0, max(intTab[, names(obj), with = FALSE])),
                         type = "n", xlab = if (xNum) xBy else "", ylab = "Intensity", xaxt = "n"), plotArgs))
    
    if (xNum)
        axis(1)
    else if (xNames)
        axis(1, seq_len(uniqueN(intTab$x)), unique(intTab$x), las = 2)
    else
        axis(1, seq_len(maxX), seq_len(maxX))

    linesArgs <- c(list(type = type, pch = pch, lty = lty), linesArgs)
    usr <- par("usr")
    
    makeLine <- function(grp, col, xgrp = NULL)
    {
        y <- if (!is.null(xgrp)) intTab[xgroup == xgrp][[grp]] else intTab[[grp]]
        irows <- if (is.null(xgrp)) seq_len(nrow(intTab)) else intTab[, .I[xgroup == xgrp]]
        x <- if (xNum) intTab$x[irows] else intTab$xnum[irows]
        do.call(lines, c(list(x = x, y = y, col = col), linesArgs))
        if (regression)
        {
            rl <- if (!is.null(xgrp)) regList[[grp]][[xgrp]] else regList[[grp]]
            lm <- if (!is.null(rl)) rl[["lm"]]
            slope <- if (!is.null(rl)) rl$slope else NA
            if (!is.null(lm) && !is.na(slope))
            {
                # from https://stackoverflow.com/a/10046370
                clip(min(x), max(x), min(y), max(y))
                abline(lm, col = col)
                do.call("clip", as.list(usr))  # reset to plot region (from ?clip examples)
            }
        }
    }
    
    for (grp in names(obj))
    {
        if (is.null(groupBy) || groupBy == "fGroups")
            makeLine(grp, col[if (is.null(groupBy)) 1 else grp])
        else
        {
            for (xgrp in unique(intTab$xgroup))
                makeLine(grp, col[xgrp], xgrp = xgrp)
        }
    }
    
    if (showLegend)
        makeLegend(par("usr")[2], par("usr")[4])
    
    par(oldp)
})

# NOTE: anaInfo is included in cases where xBy/groupBy args are used
setMethod("plotIntHash", "featureGroups", function(obj, ...) makeHash(groupTable(obj), analysisInfo(obj), ...))

#' @details \code{plotChord} Generates a chord diagram which can be used to
#'   visualize shared presence of feature groups between analyses or replicate
#'   groups. In addition, analyses/replicates sharing similar properties
#'   (\emph{e.g.} location, age, type) may be grouped to enhance visualization
#'   between these 'outer groups'.
#'
#' @param addIntraOuterGroupLinks If \code{TRUE} then links will be added within
#'   outer groups.
#'
#' @template plotChord-args
#'
#' @references \addCitations{circlize}
#'
#' @rdname feature-plotting
#' @export
setMethod("plotChord", "featureGroups", function(obj, addSelfLinks = FALSE, addRetMzPlots = TRUE, aggregate = FALSE,
                                                 groupBy = NULL, addIntraOuterGroupLinks = FALSE, IMS = "maybe",  ...)
{
    assertIMSArg(IMS)
    obj <- prepIMSFGroupsForPlot(obj, IMS)
    
    anaInfo <- analysisInfo(obj)
    
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertFlag, . ~ addSelfLinks + addRetMzPlots + addIntraOuterGroupLinks,
           fixed = list(add = ac))
    aggregateBy <- assertAndPrepareAnaInfoBy(aggregate, anaInfo, FALSE, add = ac)
    assertAnaInfoBy(groupBy, anaInfo, FALSE, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (length(obj) == 0)
        stop("Can't plot empty feature groups object")
    
    obj <- removeEmptyAnalyses(obj)
    
    gInfo <- groupInfo(obj)
    
    snames <- unique(anaInfo[[aggregateBy]])
    
    if (length(snames) < 2)
    {
        stop(sprintf("There are no features to compare among the given %s.",
                     if (isFALSE(aggregate)) "analyses" else "groups to aggregate"))
    }
    
    nsamp <- length(snames)
    
    chordTable <- rbindlist(lapply(seq_along(snames),
                                   function(sni) data.table(from = snames[sni], to = snames[seq(sni, length(snames))])))
    
    groupTab <- as.data.table(obj, average = aggregate)
    
    getLinkScore <- function(sn1, sn2)
    {
        if (sn1 == sn2)
            return(0)
        return(sum(groupTab[[getADTIntCols(sn1)]] > 0 & groupTab[[getADTIntCols(sn2)]] > 0))
    }
    
    chordTable[, value := mapply(getLinkScore, from, to)]
    
    if (addSelfLinks)
    {
        gt <- if (!isFALSE(aggregate)) averageGroups(obj, by = aggregateBy) else groupTable(obj)
        uniqueLinkCount <- sapply(seq_along(snames),
                                  function(sni) sum(sapply(gt, function(ints) ints[sni] > 0 && all(ints[-sni] == 0))))
        chordTable[from == to, value := uniqueLinkCount[.GRP], by = from]
    }
    
    if (!is.null(groupBy))
    {
        checkAnaInfoAggrGrouping(anaInfo, "aggregated", aggregateBy, groupBy)
        
        getOG <- function(s) anaInfo[match(s, anaInfo[[aggregateBy]])][[groupBy]]
        ogLookup <- anaInfo[, .(get(groupBy), get(aggregateBy))]
        chordTable[, groupFrom := getOG(from)]
        chordTable[, groupTo := getOG(to)]
        if (!addIntraOuterGroupLinks)
            chordTable[groupFrom == groupTo & from != to, value := 0] # clear links within same groups (except self links)
        setorder(chordTable, groupFrom)
        
        remainingSN <- unique(unlist(chordTable[value != 0, .(from, to)])) # assigned samples, others will be removed
        og <- getOG(remainingSN) # outer groups assigned to each remaining sample
        gaps <- rep(1, length(og)) # initialize gaps
        gaps[cumsum(sapply(unique(og), function(x) sum(og == x)))] <- 8 # make gap bigger after each outer group
        circlize::circos.par(gap.after = gaps)
    }
    
    if (all(chordTable$value == 0))
        stop("Did not found any overlap! Nothing to plot.")
    
    tracks <- NULL
    if (!is.null(groupBy))
        tracks <- list(list(track.height = 0.1, track.margin = c(if (addRetMzPlots) 0.05 else 0.06, 0)))
    if (addRetMzPlots)
        tracks <- c(tracks, list(list(track.height = 0.1, track.margin = c(0.08, 0))))
    
    maxv <- max(if (!is.null(groupBy)) chordTable[groupFrom != groupTo, value] else chordTable$value)
    colFunc <- circlize::colorRamp2(maxv * seq(0, 1, 0.25),
                                    c("blue4", "deepskyblue1", "green", "orange", "red"),
                                    transparency = 0.5)
    
    if (!is.null(groupBy) && addIntraOuterGroupLinks)
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
            ftgrps <- groupTab[get(getADTIntCols(sn)) > 0]$group
            return(gInfo[match(ftgrps, group), c("ret", "mz"), with = FALSE])
        }, simplify = FALSE), idcol = "sname")
        retMz[, ret := ret / max(ret)] # normalize
        
        circlize::circos.track(fa = retMz$sname, x = retMz$ret, y = retMz$mz, ylim = c(0, max(retMz$mz)), track.index = length(tracks),
                               panel.fun = function(x, y)
                               {
                                   x <- x / (max(x) / circlize::get.cell.meta.data("xrange"))
                                   circlize::circos.points(x, y, cex = 0.5, col = "blue", pch = 16)
                               })
    }
    
    if (!is.null(groupBy))
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
#' @param intMax Method used to determine the maximum intensity plot limit. Should be \code{"eic"} (from EIC data) or
#'   \code{"feature"} (from feature data). Ignored if the \code{ylim} parameter is specified.
#'
#' @template EICParams-arg
#'
#' @rdname feature-plotting
#' @export
setMethod("plotChroms", "featureGroups", function(obj, analysis = analyses(obj), groupName = names(obj),
                                                  retMin = FALSE, showPeakArea = FALSE, showFGroupRect = TRUE,
                                                  title = NULL, groupBy = NULL, showLegend = TRUE,
                                                  annotate = c("none", "ret", "mz", "mob"), intMax = "eic",
                                                  EICParams = getDefEICParams(), showProgress = FALSE, IMS = "maybe",
                                                  xlim = NULL, ylim = NULL, EICs = NULL, ...)
{
    ac <- checkmate::makeAssertCollection()
    assertPlotEIXArgs(obj, analysis, groupName, showPeakArea, showFGroupRect, title, groupBy, showLegend, annotate,
                      showProgress, xlim, ylim, ac)
    checkmate::assertFlag(retMin, add = ac)
    checkmate::assertChoice(intMax, c("eic", "feature"), add = ac)
    assertEICParams(EICParams, add = ac)
    assertIMSArg(IMS, add = ac)
    checkmate::reportAssertions(ac)
    
    if (intMax == "feature" && !EICParams$onlyPresent)
        stop("intMax must be 'eic' when EICParams$onlyPresent == FALSE", call. = FALSE)
    
    if (intMax == "eic")
        intMax <- "eix" # for makeEIXPlot()
    
    obj <- prepIMSFGroupsForPlot(obj, IMS)
    if (IMS != "both")
    {
        analysis <- intersect(analysis, analyses(obj))
        groupName <- intersect(groupName, names(obj))
    }
    
    if (is.null(EICs))
        EICs <- getFeatureEIXs(obj, type = "EIC", analysis, groupName, EICParams, pad = TRUE)
    else
    {
        # sync as much as possible with given EICParams
        EICs <- filterEIXs(EICs, obj, analysis = analysis, groupName = groupName, topMost = EICParams$topMost,
                           topMostByReplicate = EICParams$topMostByReplicate, onlyPresent = EICParams$onlyPresent)
    }

    # prepare EICs for plotting
    
    EICs <- lapply(EICs, function(ea)
    {
        ax <- attr(ea, "allXValues")
        lapply(ea, function(eg)
        {
            if (retMin)
                eg$time <- eg$time/60
            return(eg)
        })
    })

    takeAnalysis <- analysis
    anaInfo <- analysisInfo(obj)
    anaInfo <- anaInfo[analysis %chin% takeAnalysis & analysis %chin% names(EICs)]
    gInfo <- groupInfo(obj)
    gCount <- length(groupName)
    featTab <- as.data.table(getFeatures(obj))[group %chin% groupName]
    
    setnames(featTab, c("ret", "retmin", "retmax"), c("x", "xmin", "xmax"))
    
    if (retMin)
    {
        EICParams$window <- EICParams$window / 60
        featTab[, c("x", "xmin", "xmax") := .(x/60, xmin/60, xmax/60)]
        gInfo <- copy(gInfo)
        gInfo[, ret := ret / 60]
    }

    makeEIXPlot(featTab, anaInfo, gInfo, showPeakArea, showFGroupRect, title, groupBy, showLegend, intMax,
                showProgress, annotate, xlim, ylim, EICs, EICParams$window, EICParams$onlyPresent,
                sprintf("Retention time (%s)", if (retMin) "min." else "sec."), ...)
})

setMethod("plotChromsHash", "featureGroups", function(obj, analysis = analyses(obj), groupName = names(obj),
                                                      retMin = FALSE, showPeakArea = FALSE, showFGroupRect = TRUE,
                                                      title = NULL, groupBy = NULL, showLegend = TRUE,
                                                      annotate = c("none", "ret", "mz", "mob"),
                                                      intMax = "eic", EICParams = getDefEICParams(),
                                                      showProgress = FALSE, IMS = "maybe", xlim = NULL, ylim = NULL,
                                                      EICs = NULL, ...)
{
    annotate <- checkmate::matchArg(annotate, c("none", "ret", "mz", "mob"), several.ok = TRUE)
    if ("none" %in% annotate)
        annotate <- "none"
    args <- allArgs(FALSE)
    if (!is.null(EICs))
    {
        # omit data we don't need: speeds up hashing quite a bit
        # NOTE: only apply analysis/group filters, as the rest will slow down things considerably. Hence, this could
        # result in cache misses.
        # NOTE: we also ignore changes in analysis and groupName by selection from IMS arg
        EICs <- filterEIXs(EICs, obj, analysis = analysis, groupName = groupName, topMost = NULL,
                           topMostByReplicate = FALSE, onlyPresent = FALSE)
    }
    anas <- analysis
    makeHash(args[setdiff(names(args), c("obj", "EICs"))], EICs, featureTable(obj)[analysis], groupInfo(obj)[group %chin% groupName],
             analysisInfo(obj)[analysis %chin% anas])
})

#' @details \code{plotMobilograms} Plots extracted ion mobilograms (EIMs) of feature groups.
#' @template EIMParams-arg
#' @rdname feature-plotting
#' @export
setMethod("plotMobilograms", "featureGroups", function(obj, analysis = analyses(obj), groupName = names(obj),
                                                       showPeakArea = FALSE, showFGroupRect = TRUE, title = NULL,
                                                       groupBy = NULL, showLegend = TRUE,
                                                       annotate = c("none", "ret", "mz", "mob"),
                                                       EIMParams = getDefEIMParams(), showProgress = FALSE, xlim = NULL,
                                                       ylim = NULL, EIMs = NULL, ...)
{
    ac <- checkmate::makeAssertCollection()
    assertPlotEIXArgs(obj, analysis, groupName, showPeakArea, showFGroupRect, title, groupBy, showLegend, annotate,
                      showProgress, xlim, ylim, ac)
    assertEIMParams(EIMParams, add = ac)
    checkmate::reportAssertions(ac)
    
    if (!hasMobilities(obj))
        stop("There are no mobilities assigned to features.", call. = FALSE)
    
    if (is.null(EIMs))
        EIMs <- getFeatureEIXs(obj, type = "EIM", analysis, groupName, EIMParams)
    else
    {
        # sync as much as possible with given EIMParams
        EIMs <- filterEIXs(EIMs, obj, analysis = analysis, groupName = groupName, topMost = EIMParams$topMost,
                           topMostByReplicate = EIMParams$topMostByReplicate, onlyPresent = EIMParams$onlyPresent)
    }
    
    takeAnalysis <- analysis
    anaInfo <- analysisInfo(obj)
    anaInfo <- anaInfo[analysis %chin% takeAnalysis & analysis %chin% names(EIMs)]
    gInfo <- groupInfo(obj)
    gCount <- length(groupName)
    featTab <- as.data.table(getFeatures(obj))[group %chin% groupName]
    
    setnames(featTab, c("mobility", "mobmin", "mobmax"), c("x", "xmin", "xmax"))

    makeEIXPlot(featTab, anaInfo, gInfo, showPeakArea, showFGroupRect, title, groupBy, showLegend, "eix",
                showProgress, annotate, xlim, ylim, EIMs, EIMParams$window, EIMParams$onlyPresent, "Mobility", ...)
})

setMethod("plotMobilogramHash", "featureGroups", function(obj, analysis = analyses(obj), groupName = names(obj),
                                                          showPeakArea = FALSE, showFGroupRect = TRUE, title = NULL,
                                                          groupBy = NULL, showLegend = TRUE,
                                                          annotate = c("none", "ret", "mz", "mob"),
                                                          EIMParams = getDefEIMParams(), showProgress = FALSE, xlim = NULL,
                                                          ylim = NULL, EIMs = NULL, ...)
{
    annotate <- checkmate::matchArg(annotate, c("none", "ret", "mz", "mob"), several.ok = TRUE)
    if ("none" %in% annotate)
        annotate <- "none"
    args <- allArgs(FALSE)
    if (!is.null(EIMs))
    {
        # omit data we don't need: speeds up hashing quite a bit
        # NOTE: only apply analysis/group filters, as the rest will slow down things considerably. Hence, this could
        # result in cache misses.
        EIMs <- filterEIXs(EIMs, obj, analysis = analysis, groupName = groupName, topMost = NULL,
                           topMostByReplicate = FALSE, onlyPresent = FALSE)
    }
    anas <- analysis
    makeHash(args[setdiff(names(args), c("obj", "EIMs"))], EIMs, featureTable(obj)[analysis], groupInfo(obj)[group %chin% groupName],
             analysisInfo(obj)[analysis %chin% anas])
})

#' @details \code{plotVenn} plots a Venn diagram (using \pkg{\link{VennDiagram}}) outlining unique and shared feature
#'   groups between up to five replicates.
#' @template plotvenn-ret
#' 
#' @rdname feature-plotting
#' @export
setMethod("plotVenn", "featureGroups", function(obj, which = NULL, aggregate = TRUE, IMS = "maybe", ...)
{
    assertIMSArg(IMS)
    obj <- prepIMSFGroupsForPlot(obj, IMS)
    
    anaInfo <- analysisInfo(obj)
    aggregate <- assertAndPrepareAnaInfoBy(aggregate, anaInfo, FALSE)
    groups <- unique(anaInfo[[aggregate]])
    
    checkmate::assert(checkmate::checkSubset(which, groups, empty.ok = FALSE),
                      checkmate::checkNull(which),
                      .var.name = "which")

    if (is.null(which))
        which <- groups

    fGroupsList <- lapply(which, function(w) obj[anaInfo[get(aggregate) == w]$analysis])
    
    makeVennPlot(lapply(fGroupsList, names), which, lengths(fGroupsList), intersect, length, ...)
})

setMethod("plotVennHash", "featureGroups", function(obj, ...)
{
    return(makeHash(groupTable(obj), ...))
})

#' @details \code{plotUpSet} plots an UpSet diagram (using the \code{\link[UpSetR]{upset}} function) outlining unique
#'   and shared feature groups between given replicates.
#'   
#' @template plotUpSet
#'
#' @rdname feature-plotting
#' @export
setMethod("plotUpSet", "featureGroups", function(obj, which = NULL, aggregate = TRUE, IMS = "maybe", nsets = NULL,
                                                 nintersects = NA, ...)
{
    assertIMSArg(IMS)
    obj <- prepIMSFGroupsForPlot(obj, IMS)
    
    anaInfo <- analysisInfo(obj)
    aggregate <- assertAndPrepareAnaInfoBy(aggregate, anaInfo, FALSE)
    groups <- unique(anaInfo[[aggregate]])
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assert(checkmate::checkSubset(which, groups, empty.ok = FALSE),
                      checkmate::checkNull(which),
                      .var.name = "which", add = ac)
    checkmate::assertCount(nsets, positive = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertCount(nintersects, positive = TRUE, na.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)
    
    if (is.null(which))
        which <- groups
    
    if (length(obj) == 0)
        stop("Can't plot empty feature groups object")

    obj <- obj[anaInfo[get(aggregate) %in% which]$analysis]

    gt <- as.data.table(obj, average = aggregate)
    gt <- gt[, getADTIntCols(which), with = FALSE] # isolate relevant columns
    setnames(gt, stripADTIntSuffix(names(gt)))
    gt[, (names(gt)) := lapply(.SD, function(x) as.integer(x > 0))]
    
    if (sum(sapply(gt, function(x) any(x > 0))) < 2)
        stop("Need at least two groups with non-zero intensities", call. = FALSE)
    
    if (is.null(nsets))
        nsets <- length(which)
    
    UpSetR::upset(gt, nsets = nsets, nintersects = nintersects, ...)
})

setMethod("plotUpSetHash", "featureGroups", function(obj, ...)
{
    return(makeHash(groupTable(obj), ...))
})

#' @details \code{plotVolcano} Plots Fold change data in a 'Volcano plot'.
#' @param FCParams A parameter list to calculate Fold change data. See \code{getFCParams} for more details.
#' @rdname feature-plotting
#' @export
setMethod("plotVolcano", "featureGroups", function(obj, FCParams, showLegend = TRUE, averageFunc = mean,
                                                   normalized = FALSE, IMS = "maybe", col = NULL, pch = 19, ...)
{
    ac <- checkmate::makeAssertCollection()
    assertFCParams(FCParams, obj, null.ok = FALSE, add = ac)
    checkmate::assertFlag(showLegend, add = ac)
    checkmate::assertFunction(averageFunc, add = ac)
    checkmate::assertFlag(normalized, add = ac)
    assertIMSArg(IMS, add = ac)
    checkmate::reportAssertions(ac)

    obj <- prepIMSFGroupsForPlot(obj, IMS)
    
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
            ret <- sprintf("<b>%s</b><br>RT: %.2f<br>m/z: %.4f<br>ISTD: %s", grp, gInfo[grp == group]$ret,
                           gInfo[grp == group]$mz, istds)
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
        edges[, value := abs(gInfo[match(from, group)]$ret - gInfo[match(to, group)]$ret)]
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

#' @details \code{plotTICs} Plots the total ion chromatogram/s (TICs) of the analyses.
#' @rdname feature-plotting
#' @export
setMethod("plotTICs", "featureGroups", function(obj, retentionRange = NULL, MSLevel = 1, retMin = FALSE, title = NULL,
                                                groupBy = NULL, showLegend = TRUE, xlim = NULL, ylim = NULL, ...)
{
    plotTICs(obj@features, retentionRange, MSLevel, retMin, title, groupBy, showLegend, xlim, ylim, ...)
})

#' @details \code{plotTICs} Plots the base peak chromatogram/s (BPCs) of the analyses.
#' @rdname feature-plotting
#' @export
setMethod("plotBPCs", "featureGroups", function(obj, retentionRange = NULL, MSLevel = 1, retMin = FALSE, title = NULL,
                                                groupBy = NULL, showLegend = TRUE, xlim = NULL, ylim = NULL, ...)
{
    plotBPCs(obj@features, retentionRange, MSLevel, retMin, title, groupBy, showLegend, xlim, ylim, ...)
})
