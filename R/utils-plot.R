# SPDX-FileCopyrightText: 2016 - 2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

prepIMSFGroupsForPlot <- function(fGroups, IMS)
{
    if (is.null(IMS) || length(fGroups) == 0 || !hasMobilities(fGroups))
        return(fGroups)
    
    ret <- selectIMSFilter(fGroups, IMS = IMS, verbose = FALSE)
    if (length(ret) == 0) # this could only happened due to the IMS subsetting
        warning("No feature groups available to plot, perhaps you want to change the IMS argument?", call. = FALSE)
    return(ret)
}

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

saveRCDKStructure <- function(molecule, format, out, width = 500L, height = 500L, transparent = TRUE)
{
    if (!"fillToFit" %in% formalArgs(rcdk::get.depictor))
    {
        stopifnot(format == "svg")
        # HACK: until https://github.com/CDK-R/cdkr/pull/142 is merged we fallback to some direct CDK hackery
        depg <- rJava::.jnew("org/openscience/cdk/depict/DepictionGenerator")
        depg <- rJava::.jcall(depg, "Lorg/openscience/cdk/depict/DepictionGenerator;", "withZoom", 1.3)
        depg <- rJava::.jcall(depg, "Lorg/openscience/cdk/depict/DepictionGenerator;", "withSize", as.numeric(width),
                              as.numeric(height))
        depg <- rJava::.jcall(depg, "Lorg/openscience/cdk/depict/DepictionGenerator;", "withFillToFit")
        depg <- rJava::.jcall(depg, "Lorg/openscience/cdk/depict/DepictionGenerator;", "withAtomColors")
        dep <- rJava::.jcall(depg, "Lorg/openscience/cdk/depict/Depiction;", "depict", molecule)
        cat(dep$toSvgStr(), file = out)    
    }
    else
    {
        dep <- rcdk::get.depictor(width, height, style = "cow", fillToFit = TRUE)
        mi <- rJava::.jnew("org/guha/rcdk/view/MoleculeImage", molecule, dep)
        writeBin(rJava::.jcall(mi, "[B", "getBytes", as.integer(width), as.integer(height), format), out)
    }
    
    if (transparent)
    {
        # HACK: make it transparent
        lines <- readLines(out)
        lines <- sub("fill='#FFFFFF'", "fill='none'", lines)
        writeLines(lines, out)
    }
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
    {
        gRet <- do.call(VennDiagram::draw.single.venn, c(list(area = areas), vennArgs))
        icounts <- numeric()
    }
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
    else
    {
        warning("Venn plot only supports up to 5 categories", call. = FALSE)
        return(invisible(NULL))
    }
    
    invisible(list(gList = gRet, areas = areas, intersectionCounts = icounts))
}

makeEIXPlot <- function(featPlotTab, anaInfo, gInfo, showPeakArea, showFGroupRect, title, groupBy, showLegend, intMax,
                        showProgress, annotate, xlim, ylim, EIXs, window, onlyPresent, xlab, ...)
{
    # NOTE: when called by plotChroms(), all RTs in featPlotTab, gInfo and EIXs may have been converted to minutes
    
    if (showLegend && is.null(groupBy))
        showLegend <- FALSE
    
    if (nrow(featPlotTab) == 0 || length(EIXs) == 0)
    {
        noDataPlot()
        return(invisible(NULL))
    }
    
    gNames <- unique(featPlotTab$group)
    gCount <- length(gNames)

    if (is.null(groupBy))
        EIXColors <- "blue"
    else if (groupBy == "fGroups")
    {
        EIXColors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(gCount)
        names(EIXColors) <- gNames
    }
    else
    {
        colgrps <- unique(anaInfo[[groupBy]])
        EIXColors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(colgrps))
        names(EIXColors) <- colgrps
    }
    
    fillColors <- adjustcolor(EIXColors, alpha.f = 0.35)
    names(fillColors) <- names(EIXColors)
    
    if (is.null(title))
    {
        if (gCount == 1)
        {
            title <- sprintf("Group '%s'\nrt: %.1f - m/z: %.4f", gNames[1], gInfo[group == gNames]$ret,
                             gInfo[group == gNames]$mz)
            if (!is.null(gInfo[["mobility"]]))
                title <- sprintf("%s - mob: %.2f", title, gInfo[group == gNames]$mobility)
        }
        else
            title <- sprintf("%d feature groups", gCount)
    }
    
    if (!"none" %in% annotate)
    {
        featPlotTab[, annotation := {
            antxt <- character()
            if ("ret" %in% annotate)
                antxt <- sprintf("%.1f", gInfo$ret[match(group[1], gInfo$group)])
            if ("mz" %in% annotate)
                antxt <- paste(antxt, sprintf("%.4f", gInfo$mz[match(group[1], gInfo$group)]), sep = "\n")
            if ("mob" %in% annotate && !is.null(gInfo[["mobility"]]))
                antxt <- paste(antxt, sprintf("%.2f", gInfo$mobility[match(group[1], gInfo$group)]), sep = "\n")
            antxt
        }, by = "group"]
    }
    
    if (is.null(xlim))
    {
        xlim <- if (anyNA(featPlotTab$xmin))
            range(unlist(lapply(EIXs, function(ea) lapply(ea, "[[", 1))))
        else
            c(min(featPlotTab$xmin) - window, max(featPlotTab$xmax) + window)
    }
    
    if (is.null(ylim))
    {
        if (intMax == "eix")
        {
            xRangeGrp <- split(featPlotTab[, .(xmin = min(xmin), xmax = max(xmax)), by = "group"], by = "group", keep.by = FALSE)
            plotIntMax <- max(unlist(lapply(EIXs, function(aeix) Map(names(aeix), aeix, f = function(grp, eix)
            {
                xr <- unlist(xRangeGrp[[grp]])
                if (!is.null(xr))
                {
                    eixi <- if (is.na(xr[1]))
                        eix$intensity
                    else
                    {
                        xr <- c(max(xr[1], xlim[1]), min(xr[2], xlim[2]))
                        eix[eix[[1]] %between% xr, "intensity"]
                    }
                    if (length(eixi) > 0)
                        return(max(eixi))
                }
                return(0)
            }))))
        }
        else # "feature"
            plotIntMax <- max(featPlotTab$intensity)
        ylim <- c(0, plotIntMax * 1.1)
    }
    
    if (showLegend)
    {
        makeLegend <- function(x, y, ...)
        {
            return(legend(x, y, names(EIXColors), col = EIXColors, text.col = EIXColors, lty = 1, xpd = NA, ncol = 1,
                          cex = 0.75, bty = "n", ...))
        }
        
        plot.new()
        leg <- makeLegend(0, 0, plot = FALSE)
        lw <- (grconvertX(leg$rect$w, to = "ndc") - grconvertX(0, to = "ndc"))
        lw <- min(lw, 0.5) # don't make it too wide
        withr::local_par(list(omd = c(0, 1 - lw, 0, 1), new = TRUE))
    }
    
    plot(0, type = "n", main = title, xlab = xlab, ylab = "Intensity", xlim = xlim, ylim = ylim, ...)
    
    effectiveXlim <- par("usr")[c(1, 2)]
    effectiveYlim <- par("usr")[c(3, 4)]
    
    if (showProgress)
        prog <- openProgBar(0, gCount)
    
    for (grp in gNames)
    {
        featTabGrp <- featPlotTab[group == grp]
        for (ana in names(EIXs))
        {
            if (is.null(EIXs[[ana]]) || is.null(EIXs[[ana]][[grp]]))
                next
            EIX <- EIXs[[ana]][[grp]]
            featRow <- featTabGrp[analysis == ana]
            if (onlyPresent && nrow(featRow) == 0)
                next
            
            colInd <- if (is.null(groupBy))
                1
            else if (groupBy == "fGroups")
                grp
            else
                anaInfo[match(ana, analysis)][[groupBy]]
            
            points(EIX[[1]], EIX$intensity, type = "l", col = EIXColors[colInd])
            
            if (showPeakArea && nrow(featRow) != 0 && !is.na(featRow$xmin))
            {
                EIXFill <- setDT(EIX[numGTE(EIX[[1]], featRow$xmin) & numLTE(EIX[[1]], featRow$xmax), ])
                EIXFill <- EIXFill[EIXFill[[1]] %inrange% effectiveXlim]
                # filling doesn't work if outside y plot range
                EIXFill[intensity < effectiveYlim[1], intensity := effectiveYlim[1]]
                EIXFill[intensity > effectiveYlim[2], intensity := effectiveYlim[2]]
                polygon(c(EIXFill[[1]], rev(EIXFill[[1]])), c(EIXFill$intensity, rep(0, length(EIXFill$intensity))),
                        col = fillColors[colInd], border = NA)
            }
        }

        if ((showFGroupRect || !is.null(featTabGrp[["annotation"]])) && !anyNA(featTabGrp$xmin))
        {
            xRange <- c(min(featTabGrp$xmin), max(featTabGrp$xmax))
            maxInt <- if (intMax == "eix")
            {
                eixag <- EIXs[[ana]][[grp]]
                ints <- eixag[eixag[[1]] %between% xRange, "intensity"]
                if (length(ints) > 0) max(ints) else 0
            }
            else
                max(featTabGrp$intensity)
            
            if (showFGroupRect)
                rect(xRange[1], 0, xRange[2], maxInt, border = "red", lty = "dotted")
            
            if (!is.null(featTabGrp[["annotation"]]) && nzchar(featTabGrp$annotation))
                text(mean(featTabGrp$x), maxInt + ylim[2] * 0.02, featTabGrp$annotation)
        }
        
        if (showProgress)
            setTxtProgressBar(prog, match(grp, gNames))
    }
    
    if (showLegend)
        makeLegend(par("usr")[2], par("usr")[4])
    
    if (showProgress)
    {
        setTxtProgressBar(prog, gCount)
        close(prog)
    }
    invisible(NULL)
}

getMSPlotData <- function(spec, marklwd, markWhich = NULL)
{
    hasFragInfo <- !is.null(spec[["ion_formula"]])
    plotData <- copy(spec)
    
    # default colour/line width
    plotData[, c("colour", "lwd", "legend") := .("grey", 1, "unassigned")]
    
    if (is.null(markWhich))
        markWhich <- if (hasFragInfo) which(plotData$annotated) else seq_len(nrow(plotData))
    else
        markWhich <- plotData[mergedBy %in% markWhich, which = TRUE]
    
    if (!is.null(spec[["mergedBy"]]))
    {
        plotData[nzchar(mergedBy) & !is.na(mergedBy), legend := sapply(mergedBy, wrapStr, width = 10)]
        plotData[markWhich, lwd := marklwd]
        
        mbsUnique <- unique(plotData$legend)
        # order from small to big based on number of commas
        mbsUnique <- mbsUnique[order(sapply(mbsUnique, countCharInStr, ch = ",", USE.NAMES = FALSE))]
        mbCombCols <- setNames(getBrewerPal(length(mbsUnique), "Paired"), mbsUnique)
        
        plotData[, colour := mbCombCols[match(legend, mbsUnique)]]
    }
    else if (hasFragInfo)
        plotData[markWhich, c("colour", "lwd", "legend") := .("blue", marklwd, "assigned")] # nothing merged, just mark all annotated blue
    
    # mark precursor
    plotData[precursor == TRUE, c("colour", "lwd", "legend") := .("red", marklwd, "precursor")]
    
    return(plotData)
}

makeScoresPlot <- function(scoreTable, mcn)
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
    
    isMerged <- length(mcn) > 1 && uniqueN(scores$merged) > 1
    
    oldp <- par(no.readonly = TRUE)
    maxStrW <- max(strwidth(unique(scores$type), units = 'in', cex = 0.9)) + 0.5
    omai <- par("mai")
    par(mai = c(maxStrW, 0.5, omai[3], 0))
    
    if (isMerged)
        cols <- getBrewerPal(length(unique(scores$merged)), "Paired")
    else
        cols <- getBrewerPal(nrow(scores), "Paired")
    
    bpargs <- list(las = 2, col = cols, border = cols, cex.axis = 0.9, xpd = TRUE)
    
    if (isMerged)
    {
        scSplit <- split(scores, by = "type", keep.by = FALSE)
        scSplit <- sapply(names(scSplit), function(mb) setnames(scSplit[[mb]], "score", mb), simplify = FALSE) # assign column names
        
        plotTab <- Reduce(function(left, right)
        {
            merge(left, right, by = "merged", all = TRUE, sort = FALSE)
        }, scSplit)
        
        plot.new()
        
        makeLegend <- function(x, y, ...) legend(x, y, unique(plotTab$merged), col = cols, lwd = 1, xpd = NA, ncol = 1,
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

# spec may be annotated
makeMSPlot <- function(plotData, mincex, xlim, ylim, ylab = "Intensity", ..., mol = NULL, maxMolSize, molRes)
{
    if (is.null(xlim))
        xlim <- range(plotData$mz) * c(0.9, 1.1)
    else
    {
        plotData <- plotData[numGTE(mz, xlim[1]) & numLTE(mz, xlim[2])] # remove any peaks outside plotting range
        if (nrow(plotData) == 0)
            return(noDataPlot())
    }
    
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
        expand <- if (!is.null(mol))
            1.75
        else if (!is.null(plotData[["ion_formula"]]) && any(!is.na(plotData$ion_formula)))
            1.5
        else
            1.2
        
        if (min(plotData$intensity) < 0) # mirror plot?
        {
            # extend both vertical directions 
            if (max(plotData$intensity) > 0)
                ylim <- range(plotData$intensity) * expand
            else # only bottom plot
                ylim <- c(min(plotData$intensity), abs(min(plotData$intensity)))
        }
        else
        {
            ylim <- if (nrow(plotData) > 1)
                range(plotData$intensity)
            else
                c(0, max(plotData$intensity))
            ylim <- ylim * c(1, expand)
        }
    }

    ticks <- pretty(c(-plotData$intensity, plotData$intensity))
    
    plot(0, xlab = "m/z", ylab = ylab, xlim = xlim, ylim = ylim, yaxt = "n", type = "n", bty = "l", ...)
    
    segments(plotData[["mz"]], 0, plotData[["mz"]], plotData[["intensity"]], col = plotData[["colour"]],
             lwd = plotData[["lwd"]])

    axis(2, at = ticks, labels = abs(ticks))
    
    if (min(plotData$intensity) < 0) # add horizontal line for mirror plots
        abline(h = 0, col = "darkgrey")
    
    annPlotData <- NULL
    if (!is.null(plotData[["ion_formula"]]))
    {
        annPlotData <- plotData[!is.na(ion_formula)]
        if (nrow(annPlotData) > 0)
        {
            formsSub <- subscriptFormula(annPlotData$ion_formula)
            
            calcFormHeights <- function(cex)
            {
                # get widths in inches which are device/direction independent
                formWidths <- strwidth(formsSub, cex = cex, units = "inch")
                # labels are rotated 90 deg: convert widths to heights
                return(yinch(formWidths))
            }
            
            curcex <- par("cex")
            formHeights <- calcFormHeights(curcex)
            totHeights <- abs(annPlotData$intensity) + formHeights
            
            # use formula text height from highest total height 
            maxFormHeight <- formHeights[which.max(totHeights)]
            
            # how high the text can be depends on ylim (+some spacing)
            maxHeight <- min(abs(ylim[1] - min(annPlotData$intensity)),
                             ylim[2] - max(annPlotData$intensity)) * 0.98
            
            # scale cex if necessary
            cex <- max(min(curcex, curcex * (maxHeight / maxFormHeight)), mincex)
            
            if (cex > 0)
            {
                formHeights <- calcFormHeights(cex) # update with new cex
                yoffsets <- ylim[2] * 0.02
                yoffsets <- ifelse(annPlotData$intensity < 0, -(yoffsets + formHeights), yoffsets)
                text(annPlotData$mz, annPlotData$intensity + yoffsets, formsSub, srt = 90, adj = 0, cex = cex)
                annPlotData[, formHeight := formHeights]
            }
            else
                annPlotData[, formHeight := 0]
        }
    }
    
    if (doLegend)
    {
        makeLegend(par("usr")[2], par("usr")[4])
        par(oldp)
    }
    
    # draw structure
    if (!is.null(mol))
    {
        img <- getRCDKStructurePlot(mol[[1]], molRes[1], molRes[2])
        
        dpi <- (par("cra")/par("cin"))[1]
        
        xl <- par("usr")[2]
        yl <- par("usr")[4]
        imgInfo <- magick::image_info(img)
        
        imgPlotW <- xinch(imgInfo$width / dpi)
        imgPlotH <- yinch(imgInfo$height / dpi)
        
        maxW <- maxMolSize[1] * xl
        if (imgPlotW > maxW)
        {
            hresize <- imgPlotW / maxW
            imgPlotH <- imgPlotH / hresize
            imgPlotW <- maxW
        }
        
        # offset a little
        startx <- par("usr")[1] + 0.01 * xl
        yl <- yl * 0.99
        
        # calculate available height
        pd <- plotData[numGTE(mz, startx) & numLTE(mz, startx + imgPlotW)]
        availHeight <- ylim[2] - if(nrow(pd) > 0) max(pd$intensity) else 0
        
        if (!is.null(annPlotData))
        {
            # also take annotations into account
            pd <- annPlotData[numGTE(mz, startx) & numLTE(mz, startx + imgPlotW)]
            if (nrow(pd) > 0)
                availHeight <- min(availHeight, ylim[2] - max(pd$intensity + pd$formHeight))
        }

        maxH <- min(maxMolSize[2] * yl, availHeight)
        if (maxH > 0)
        {
            if (imgPlotH > maxH)
            {
                wresize <- imgPlotH / maxH
                imgPlotW <- imgPlotW / wresize
                imgPlotH <- maxH
            }
            rasterImage(img, startx, yl - imgPlotH, startx + imgPlotW, yl)
        }
    }
}

getMSPlotDataOverlay <- function(specs, mirror, normalize, marklwd, markWhich)
{
    specs <- lapply(specs, copy)
    if (normalize)
        specs <- lapply(specs, function(sp) sp[, intensity := normalize(intensity, minMax = FALSE)])
    
    if (mirror && length(specs) == 2)
        specs[[2]][, intensity := -intensity]
    
    combinedSpec <- rbindlist(specs, fill = TRUE) # columns may be different due to fragInfos
    setorderv(combinedSpec, "intensity")
    
    plotData <- getMSPlotData(combinedSpec, marklwd, markWhich)
    return(plotData)
}

makeMSPlotOverlay <- function(plotData, title, mincex, xlim, ylim, ..., mol = NULL, maxMolSize = NULL,
                              molRes = NULL)
{
    makeMSPlot(plotData, mincex, xlim, ylim, ylab = "Normalized intensity",
               main = title, ..., mol = mol, maxMolSize = maxMolSize, molRes = molRes)
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

textPlot <- function(txt)
{
    withr::with_par(list(mar = c(0, 2, 0, 0)), {
        plot(1:10, 1:10, ann = FALSE, xaxt = "n", yaxt = "n", xlab = "", ylab = "", type = "n", adj = 1, bty = "n") # empty dummy plot
        text(1, 5, txt, adj = 0, cex = 0.8)
    })
}

noDataPlot <- function() textPlot("no data to plot")

withSVGLite <- withr::with_(function(file, ...) svglite::svglite(file, ...), function(old) dev.off())

doPlotHeaders <- function(obj, what = "tic", retentionRange, MSLevel, retMin = FALSE, title = NULL, 
                          groupBy = NULL, showLegend = TRUE, xlim = NULL, ylim = NULL, ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertString(title, null.ok = TRUE, add = ac)
    assertAnaInfoBy(groupBy, obj, FALSE, null.ok = TRUE, add = ac)
    assertXYLim(xlim, ylim, add = ac)
    checkmate::reportAssertions(ac)
    
    if (showLegend && is.null(groupBy))
        showLegend <- FALSE
    
    if ("tic" %in% what)
        data <- getTICs(obj, retentionRange, MSLevel)
    else if ("bpc" %in% what)
        data <- getBPCs(obj, retentionRange, MSLevel)
    else
        data <- data.table()
    
    if (nrow(data) == 0)
    {
        noDataPlot()
        return(invisible(NULL))
    }
    
    obj <- obj[obj$analysis %chin% unique(data$analysis), ]
    anas <- obj$analysis
    anaCount <- length(anas)
    reps <- unique(obj$replicate)
    
    if (is.null(groupBy))
        plotColors <- "blue"
    else
    {
        colgrps <- unique(obj[[groupBy]])
        plotColors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(colgrps))
        names(plotColors) <- colgrps
    }
    
    fillColors <- adjustcolor(plotColors, alpha.f = 0.35)
    names(fillColors) <- names(plotColors)
    
    if (is.null(xlim))
    {
        xlim <- c(min(data$ret), max(data$ret))
    }
    if (is.null(ylim))
    {
        plotIntMax <- max(data$intensity)
        ylim <- c(0, plotIntMax * 1.1)
    }
    
    if (is.null(title))
    {
        if (anaCount == 1)
            title <- sprintf("Analysis '%s'", anas[1])
        else
            title <- sprintf("%d analyses", anaCount)
    }
    
    if (showLegend)
    {
        makeLegend <- function(x, y, ...)
        {
            return(legend(x, y, names(plotColors), col = plotColors, text.col = plotColors, lty = 1, xpd = NA, ncol = 1,
                          cex = 0.75, bty = "n", ...))
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
    
    for (ana in anas)
    {
        anadt <- data[analysis %in% ana, ]
        
        if (nrow(anadt) == 0)
            next

        if (is.null(groupBy))
            colInd <- 1
        else
            colInd <- obj[match(ana, analysis)][[groupBy]]
        
        points(if (retMin) anadt$ret / 60 else anadt$ret, anadt$intensity, type = "l", col = plotColors[colInd])
    }
    
    if (showLegend)
        makeLegend(par("usr")[2], par("usr")[4])
    
}
