# based on http://www.cureffi.org/2013/09/23/a-quick-intro-to-chemical-informatics-in-r/
getRCDKStructurePlot <- function(molecule, width = 500, height = 500, trim = TRUE, transparent = TRUE)
{
    # drawing unconnected is not supported. Assume these are salts and we can just draw largest instead
    # UNDONE: is this still relevant?
    # if (!is.connected(molecule))
    #     molecule <- get.largest.component(molecule)
    
    img <- rcdk::view.image.2d(molecule, rcdk::get.depictor(width, height)) # get Java representation into an image matrix.
    img <- magick::image_read(img)
    
    if (!isEmptyMol(molecule))
    {
        if (trim)
            img <- magick::image_trim(img)
        if (transparent)
            img <- magick::image_transparent(img, "white")
    }
    
    return(img)
}

getBrewerPal <- function(n, name)
{
    maxn <- RColorBrewer::brewer.pal.info[name, "maxcolors"]
    if (n > maxn)
        return(colorRampPalette(RColorBrewer::brewer.pal(maxn, name))(n))
    
    if (n < 3)
        return(RColorBrewer::brewer.pal(3, name)[seq_len(n)])
    
    return(RColorBrewer::brewer.pal(n, name))
}

makeVennPlot <- function(plotObjects, categories, areas, intersectFunc,
                         intersectLenFunc, ...)
{
    nobj <- length(plotObjects)
    areas <- unname(areas) # do.call below won't work with names
    
    if (all(areas == 0))
        stop("Cannot plot Venn when all objects are empty")
    
    fill <- getBrewerPal(nobj, "Paired")
    vennArgs <- list(category = categories, lty = rep("blank", nobj), alpha = rep(0.5, nobj),
                     cex = 1.5, cat.cex = 1.5, fill = fill)
    vennArgs <- modifyList(vennArgs, list(...))
    
    getIntersectCounts <- function(inters) sapply(inters,
                                                  function(i) intersectLenFunc(Reduce(intersectFunc, plotObjects[i])))
    
    grid::grid.newpage() # need to clear plot region manually
    # plot.new()
    
    if (nobj == 1)
        gRet <- do.call(draw.single.venn, c(list(area = areas), vennArgs))
    else if (nobj == 2)
    {
        icounts <- getIntersectCounts(list(c(1, 2)))
        gRet <- do.call(VennDiagram::draw.pairwise.venn,
                        c(areas, icounts, list(rotation.degree = if (areas[1] < areas[2]) 180 else 0), vennArgs))
    }
    else if (nobj == 3)
    {
        icounts <- getIntersectCounts(list(c(1, 2), c(2, 3), c(1, 3), c(1, 2, 3)))
        gRet <- do.call(VennDiagram::draw.triple.venn, c(areas, icounts, vennArgs))
    }
    else if (nobj == 4)
    {
        icounts <- getIntersectCounts(list(c(1, 2), c(1, 3), c(1, 4), c(2, 3), c(2, 4), c(3, 4),
                                           c(1, 2, 3), c(1, 2, 4), c(1, 3, 4), c(2, 3, 4), c(1, 2, 3, 4)))
        gRet <- do.call(VennDiagram::draw.quad.venn, c(areas, icounts, vennArgs))
    }
    else if (nobj == 5)
    {
        icounts <- getIntersectCounts(list(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(2, 3), c(2, 4), c(2, 5), c(3, 4),
                                           c(3, 5), c(4, 5), c(1, 2, 3), c(1, 2, 4), c(1, 2, 5), c(1, 3, 4), c(1, 3, 5),
                                           c(1, 4, 5), c(2, 3, 4), c(2, 3, 5), c(2, 4, 5), c(3, 4, 5),
                                           c(1, 2, 3, 4), c(1, 2, 3, 5), c(1, 2, 4, 5), c(1, 3, 4, 5), c(2, 3, 4, 5),
                                           c(1, 2, 3, 4, 5)))
        gRet <- do.call(VennDiagram::draw.quintuple.venn, c(areas, icounts, vennArgs))
    }
    
    invisible(list(gList = gRet, areas = areas, intersectionCounts = icounts))
}

getMSPlotData <- function(spec)
{
    hasFragInfo <- !is.null(spec[["formula"]])
    plotData <- copy(spec)
    
    # default colour/line width
    plotData[, c("colour", "lwd", "legend") := .("grey", 1, "unassigned")]
    
    if (hasFragInfo)
    {
        isAnnotated <- plotData[!is.na(formula), which = TRUE]
        if (!is.null(spec[["mergedBy"]]))
        {
            plotData[isAnnotated, legend := sapply(mergedBy, wrapStr, width = 10)]
            
            mbsUnique <- unique(plotData$legend)
            # order from small to big based on number of commas
            mbsUnique <- mbsUnique[order(sapply(mbsUnique, countCharInStr, ch = ",", USE.NAMES = FALSE))]
            mbCombCols <- setNames(getBrewerPal(length(mbsUnique), "Paired"), mbsUnique)
            
            plotData[isAnnotated, c("colour", "lwd") := .(mbCombCols[match(legend, mbsUnique)], 2)]
        }
        else
            plotData[isAnnotated, c("colour", "lwd", "legend") := .("blue", 2, "assigned")] # nothing merged, just mark all annotated blue
    }
    
    # mark precursor
    plotData[precursor == TRUE, c("colour", "lwd", "legend") := .("red", 2, "precursor")]
    
    return(plotData)
}

makeScoresPlot <- function(scoreTable, mcn, useGGPlot2)
{
    scores <- setnames(transpose(scoreTable), "score")
    scores[, type := names(scoreTable)]
    scores <- scores[!is.na(score)]
    
    if (length(mcn) > 1)
    {
        scores[, merged := "consensus"]
        for (n in mcn)
        {
            withM <- which(grepl(paste0("-", n), scores[["type"]], fixed = TRUE))
            set(scores, withM, "merged", n)
            set(scores, withM, "type", gsub(paste0("-", n), "", scores[["type"]][withM]))
        }
    }
    
    if (!useGGPlot2)
    {
        oldp <- par(no.readonly = TRUE)
        maxStrW <- max(strwidth(unique(scores$type), units = 'in', cex = 0.9)) + 0.5
        omai <- par("mai")
        par(mai = c(maxStrW, 0.5, omai[3], 0))
        
        if (length(mcn) > 1)
            cols <- getBrewerPal(length(unique(scores$merged)), "Paired")
        else
            cols <- getBrewerPal(nrow(scores), "Paired")
        
        bpargs <- list(las = 2, col = cols, border = cols, cex.axis = 0.9, xpd = TRUE)
        
        if (length(mcn) > 1)
        {
            scSplit <- split(scores, by = "type", keep.by = FALSE)
            scSplit <- sapply(names(scSplit), function(mb) setnames(scSplit[[mb]], "score", mb), simplify = FALSE) # assign column names
            
            plotTab <- Reduce(function(left, right)
            {
                merge(left, right, by = "merged", all = TRUE)
            }, scSplit)
            
            plot.new()
            
            makeLegend <- function(x, y, ...) legend(x, y, unique(scores$merged), col = cols, lwd = 1, xpd = NA, ncol = 1,
                                                     cex = 0.75, bty = "n", ...)
            
            # auto legend positioning: https://stackoverflow.com/a/34624632/9264518
            leg <- makeLegend(0, 0, plot = FALSE)
            lw <- (grconvertX(leg$rect$w, to = "ndc") - grconvertX(0, to = "ndc"))
            par(omd = c(0, 1 - lw, 0, 1), new = TRUE)
            bpvals <- as.matrix(plotTab[, -"merged"])
            bp <- do.call(barplot, c(list(bpvals, beside = TRUE), bpargs))
            bpsc <- as.vector(bpvals)
            makeLegend(par("usr")[2], par("usr")[4])
        }
        else
        {
            bp <- do.call(barplot, c(list(scores$score, names.arg = scores$type), bpargs))
            bpsc <- scores$score
        }
        
        text(bp, bpsc, labels = round(bpsc, 2), pos = 3, cex = 0.8, xpd = TRUE)
        
        par(oldp)
    }
    else
    {
        scorePlot <- ggplot(scores, aes_string(x = "type", y = "score")) +
            cowplot::theme_cowplot(font_size = 12) +
            theme(axis.title.y = element_blank(), axis.title.x = element_blank(), # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                  legend.position = "top", legend.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) +
            guides(colour = guide_legend(nrow = 3, ncol = 2, byrow = TRUE))
        
        
        if (length(mcn) > 1)
            scorePlot <- scorePlot + geom_bar(stat = "identity", position = "dodge",
                                              aes_string(colour = "merged", fill = "merged"))
        else
        {
            scorePlot <- scorePlot +
                geom_bar(stat = "identity", aes_string(colour = "type", fill = "type")) +
                theme(legend.position = "none")
        }
        
        return(scorePlot)
    }
}

# spec may be annotated
makeMSPlot <- function(spec, xlim, ylim, ..., extraHeightInch = 0)
{
    plotData <- getMSPlotData(spec)
    if (!is.null(plotData[["formula"]]))
        plotData[!is.na(formula), formWidth := strwidth(formula, units = "inches")]
    
    if (is.null(xlim))
        xlim <- range(plotData$mz) * c(0.9, 1.1)
    else
        plotData <- plotData[numGTE(mz, xlim[1]) & numLTE(mz, xlim[2])] # remove any peaks outisde plotting range
    
    doLegend <- !is.null(plotData[["legend"]]) && any(!is.na(plotData[["legend"]]) & nzchar(plotData[["legend"]]))
    if (doLegend)
    {
        makeLegend <- function(x, y, ...)
        {
            legTab <- unique(plotData, by = "legend")[!is.na(legend)]
            # order from small to big based on number of commas
            legTab <- legTab[order(sapply(legTab$legend, countCharInStr, ch = ",", USE.NAMES = FALSE))]
            legend(x, y, legTab$legend, col = legTab$colour, xpd = NA, bty = "n",
                   text.col = legTab$colour, lty = 1, cex = 0.75, ...)
        }
        
        plot.new()
        leg <- makeLegend(0, 0, plot = FALSE)
        lw <- grconvertX(leg$rect$w, to = "ndc") - grconvertX(0, to = "ndc")
        oldp <- par(omd = c(0, 1 - lw, 0, 1), new = TRUE)
    }
    
    if (is.null(ylim))
    {
        formWidths <- if (!is.null(plotData[["formWidth"]])) plotData$formWidth else rep(0, nrow(plotData))
        formWidths[is.na(formWidths)] <- 0
        
        # see how much extra vertical space is needed by formula labels
        # get character widths (assuming that height of vertically plotted text is the same)
        pheight <- par("din")[2] # plot dev height in inches
        relFormHeights <- formWidths / pheight
        
        ym <- max(plotData$intensity) # 'regular' plot height in user coordinates, not considering any labels etc
        relIntHeights <- plotData$intensity / ym
        maxRelH <- max(relFormHeights + relIntHeights)
        
        ym <- ym * maxRelH * 1.05 # enlarge y limit and add some extra spacing
        
        if (extraHeightInch > 0)
            ym <- ym * (1 + (extraHeightInch / pheight))
        
        ylim <- c(0, ym)
    }
    
    plot(0, xlab = "m/z", ylab = "Intensity", xlim = xlim, ylim = ylim,
         type = "n", bty = "l", ...)
    
    segments(plotData[["mz"]], 0, plotData[["mz"]], plotData[["intensity"]],
             col = plotData[["colour"]], lwd = plotData[["lwd"]])
    
    if (!is.null(plotData[["formula"]]))
    {
        pdf <- plotData[!is.na(formula)]
        if (nrow(pdf) > 0)
            text(pdf[["mz"]], pdf[["intensity"]] + (ylim[2] * 0.02),
                 subscriptFormula(pdf[["formula"]]), srt = 90, adj = 0)
    }
    
    if (doLegend)
    {
        makeLegend(par("usr")[2], par("usr")[4])
        par(oldp)
    }
}

makeMSPlotGG <- function(spec, ...)
{
    plotData <- getMSPlotData(spec)
    
    # BUG: throws errors when parse=TRUE and all labels are empty
    if (!is.null(plotData$formula) && any(!is.na(plotData$formula)))
    {
        plotData[!is.na(formula), formula := subscriptFormula(formula, parse = FALSE)]
        ret <- ggplot(plotData, aes_string(x = "mz", y = 0, label = "formula")) +
            ggrepel::geom_text_repel(aes_string(y = "intensity", angle = 0), min.segment.length = 0.1, parse = TRUE,
                                     nudge_y = grid::convertUnit(grid::unit(5, "mm"), "npc", valueOnly = TRUE), size = 3.2,
                                     na.rm = TRUE)
    }
    else
        ret <- ggplot(plotData, aes_string(x = "mz", y = 0))
    
    ret <- ret + xlim(range(spec$mz) * c(0.9, 1.1)) +
        geom_segment(aes_string(xend = "mz", yend = "intensity", colour = "legend", size = "lwd")) +
        scale_size(range = c(0.5, 2), guide = FALSE) + xlab("m/z") + ylab("Intensity") +
        cowplot::theme_cowplot(font_size = 12) + theme(legend.position = "bottom", legend.title = element_blank())
    
    return(ret)
}

plotDendroWithClusters <- function(dendro, ct, pal, colourBranches, showLegend, ...)
{
    if (colourBranches)
    {
        ct <- ct[order.dendrogram(dendro)] # re-order for dendrogram
        nclust <- length(unique(ct[ct != 0])) # without unassigned (in case of dynamicTreeCut)
        cols <- getBrewerPal(nclust, pal)
        dendro <- dendextend::color_branches(dendro, clusters = ct, col = cols[unique(ct)]) # unique: fixup colour order
        lcols <- dendextend::get_leaves_branches_col(dendro)
        dendextend::labels_colors(dendro) <- lcols
    }
    
    if (colourBranches && showLegend)
    {
        mr <- par("mar")
        mr[4] <- max(5.5, mr[4])
        withr::with_par(list(mar = mr),
                        {
                            plot(dendro, ...)
                            legend("topright", legend = seq_len(nclust),
                                   bty = "n", cex = 1, fill = cols, inset = c(-0.22, 0), xpd = NA,
                                   ncol = 2, title = "cluster")
                        })
    }
    else
        plot(dendro, ...)
}

doPlotSilhouettes <- function(clust, distm, kSeq, pch = 16, type = "b", ...)
{
    if (length(distm) == 1)
        stop("Need >2 clustered compounds")
    
    meanws <- sapply(kSeq, function(k)
    {
        sil <- silhouette(cutree(clust, k = k), distm)
        summary(sil)$avg.width
    })
    
    plot(x = kSeq, y = meanws, pch = pch, type = type, xaxt = "none",
         xlab = "cluster k", ylab = "mean silhouette width", ...)
    axis(1, kSeq)
    abline(v = kSeq[which.max(meanws)], lty = 2)
}
