#' @include main.R
#' @include report-html.R
NULL

getHTMLReportPlotPath <- function(outPath)
{
    plotPath <- file.path(outPath, "report_files", "plots")
    mkdirp(plotPath)
    return(plotPath)
}

getHTMLReportFinalPlotPath <- function(out, selfContained)
{
    if (selfContained)
        return(knitr::image_uri(out))
    return(file.path("report_files", "plots", basename(out)))
}

makeHTMLReportPlot <- function(out, outPath, selfContained, code, ...)
{
    if (FALSE)
    {
        # UNDONE: while embedding the SVG directly would be nice, this seems to give major headaches with scaling,
        # especially with Firefox... For now just base64 it :(
        
        if (selfContained)
        {
            svgstr <- svglite::svgstring(standalone = FALSE, fix_text_size = FALSE, ...)
            on.exit(dev.off(), add = TRUE)
            force(code)
            ret <- as.character(svgstr())
            
            # replace fixed width/height properties to allow proper scaling (see https://stackoverflow.com/a/45144857).
            # NOTE: use sub so that only header (first occurrence) is modified. Furthermore, note that the svglite css class
            # is changed in report.Rmd.
            # ret <- sub("width=\\'[[:graph:]]+\\'", "width='100%'", ret)
            # ret <- sub("height=\\'[[:graph:]]+\\'", "height='auto'", ret)
            return(ret)
        }
    }
    
    withSVGLite(file.path(getHTMLReportPlotPath(outPath), out), standalone = TRUE, code = code, ...)
    return(getHTMLReportFinalPlotPath(out, selfContained))
    
    # UNDONE: object tag makes text selectable but messes up layout...
    # return(paste0("<object data='", out, "' type='image/svg+xml' width=500 height=300></object>"))
}

genHTMLReportPlotsChromsLarge <- function(fGroups, outPath, EICs, selfContained)
{
    cat("Generate large chromatograms...\n")
    # UNDONE: parallel option
    doApply("sapply", TRUE, names(fGroups), function(grp)
    {
        doProgress()
        makeHTMLReportPlot(paste0("chrom_large-", grp, ".svg"), outPath, selfContained, {
            # UNDONE: params
            mar <- par("mar")
            par(mar = c(4.1, mar[2], 0.2, 0.2))
            plotChroms(fGroups, groupName = grp, rtWindow = 30, mzExpWindow = 0.005, retMin = TRUE, topMost = 1,
                       topMostByRGroup = TRUE, EICs = EICs, colourBy = "rGroups", title = "",
                       onlyPresent = TRUE, bty = "l")
        }, width = 6, height = 4, bg = "transparent", pointsize = 16)
    }, simplify = FALSE)
}

genHTMLReportPlotsChromsSmall <- function(fGroups, outPath, EICs, selfContained)
{
    cat("Generate small chromatograms...\n")
    # UNDONE: parallel option
    doApply("sapply", TRUE, names(fGroups), function(grp)
    {
        doProgress()
        makeHTMLReportPlot(paste0("chrom_small-", grp, ".svg"), outPath, selfContained, {
            # UNDONE: params
            par(mai = c(0, 0, 0, 0), lwd = 10)
            plotChroms(fGroups, groupName = grp, rtWindow = 30, mzExpWindow = 0.005, retMin = TRUE, topMost = 1,
                       topMostByRGroup = TRUE, EICs = EICs, showFGroupRect = FALSE, showPeakArea = TRUE,
                       title = "", onlyPresent = TRUE, bty = "n")
        }, width = 12, height = 4, bg = "transparent", pointsize = 16)
    }, simplify = FALSE)
}

genHTMLReportPlotsChromsFeatures <- function(fGroups, outPath, EICs, selfContained)
{
    anas <- analyses(fGroups)
    cat("Generate individual feature chromatograms...\n")
    # UNDONE: parallel option
    doApply("sapply", TRUE, names(fGroups), function(grp)
    {
        doProgress()
        # whana <- fGroups[[grp]] > 0
        whana <- rep(TRUE, length(anas))
        sapply(anas[whana], function(ana)
        {
            makeHTMLReportPlot(paste0("chrom-", grp, "-", ana, ".svg"), outPath, selfContained, {
                # UNDONE: params
                mar <- par("mar")
                par(mar = c(4.1, mar[2], 0.2, 0.2))
                plotChroms(fGroups, analysis = ana, groupName = grp, rtWindow = 30, mzExpWindow = 0.005, retMin = TRUE,
                           EICs = EICs, showFGroupRect = FALSE, showPeakArea = TRUE, title = "", onlyPresent = FALSE,
                           bty = "l")
            }, width = 6, height = 4, bg = "transparent", pointsize = 20, scaling = 1)
        })
    }, simplify = FALSE)
}

genHTMLReportPlotsStructs <- function(fGroups, compounds, outPath, selfContained)
{
    scrStructInfo <- if (isScreening(fGroups)) screenInfo(fGroups)[, c("SMILES", "InChIKey"), with = FALSE] else NULL
    compStructInfo <- if (!is.null(compounds) && length(compounds) != 0) as.data.table(compounds)[, c("SMILES", "InChIKey"), with = FALSE] else NULL
    
    structInfo <- rbindlist(list(scrStructInfo, compStructInfo))
    if (nrow(structInfo) > 0)
    {
        structInfo <- unique(structInfo, by = "InChIKey")
        cat("Generate structures...\n")
        # UNDONE: parallel option
        return(doApply("Map", TRUE, structInfo$InChIKey, structInfo$SMILES, f = function(ik, smi)
        {
            pf <- file.path(getHTMLReportPlotPath(outPath), paste0("struct-", ik, ".svg"))
            saveRCDKStructure(getMoleculesFromSMILES(smi)[[1]], "svg", pf, 100, 100)
            return(getHTMLReportFinalPlotPath(pf, selfContained))
            doProgress()
        }))
    }
    return(list())
}

genHTMLReportPlotsFormulas <- function(formulas, MSPeakLists, outPath, selfContained)
{
    cat("Generate formula annotation plots...\n")
    # UNDONE: parallel option
    return(doApply("Map", TRUE, groupNames(formulas), annotations(formulas), f = function(grp, ann)
    {
        ret <- list()
        
        ret$spectra <- sapply(seq_len(nrow(ann)), function(index)
        {
            if (is.null(MSPeakLists[[grp]][["MSMS"]]))
                return("")
            makeHTMLReportPlot(sprintf("form-spec-%s-%d.svg", grp, index), outPath, selfContained, {
                plotSpectrum(formulas, index, grp, MSPeakLists = MSPeakLists)
            }, width = 7, height = 4)
        })
        
        ret$scores <- sapply(seq_len(nrow(ann)), function(index)
        {
            makeHTMLReportPlot(sprintf("form-scores-%s-%d.svg", grp, index), outPath, selfContained, {
                plotScores(formulas, index, grp) # UNDONE: params
            }, width = 6, height = 5)
        })
        
        doProgress()
        
        return(ret)
    }))
}

genHTMLReportPlotsCompounds <- function(compounds, MSPeakLists, formulas, outPath, selfContained)
{
    cat("Generate compound annotation plots...\n")
    # UNDONE: parallel option
    return(doApply("Map", TRUE, groupNames(compounds), annotations(compounds), f = function(grp, ann)
    {
        ret <- list()
        
        ret$spectra <- sapply(seq_len(nrow(ann)), function(index)
        {
            makeHTMLReportPlot(sprintf("comp-spec-%s-%d.svg", grp, index), outPath, selfContained, {
                mar <- par("mar")
                par(mar = c(mar[1], mar[2], 0.2, 0.2))
                plotSpectrum(compounds, index, grp, MSPeakLists, formulas, FALSE, title = "")
            }, width = 7, height = 4, pointsize = 16)
        })
        
        ret$scores <- sapply(seq_len(nrow(ann)), function(index)
        {
            makeHTMLReportPlot(sprintf("comp-scores-%s-%d.svg", grp, index), outPath, selfContained, {
                mar <- par("mar")
                par(mar = c(mar[1], mar[2], 0.2, 0.2))
                plotScores(compounds, index, grp) # UNDONE: params
            }, width = 7, height = 5, pointsize = 16)
        })
        
        doProgress()
        
        return(ret)
    }))
}

genHTMLReportPlotsCompsCluster <- function(compsCluster, outPath, selfContained)
{
    cat("Generate compound cluster plots...\n")
    # UNDONE: parallel option
    return(doApply("Map", TRUE, groupNames(compsCluster), cutClusters(compsCluster), f = function(grp, ct)
    {
        ret <- list()
        
        ret$dendro <- makeHTMLReportPlot(sprintf("comp-clust-dendro-%s.svg", grp), outPath, selfContained, {
            plot(compsCluster, groupName = grp)
        }, width = 12, height = 4, pointsize = 16)
        
        ret$mcs <- sapply(sort(unique(ct)), function(cli)
        {
            makeHTMLReportPlot(sprintf("comp-clust-mcs-%s-%d.svg", grp, cli), outPath, selfContained, {
                plotStructure(compsCluster, groupName = grp, cluster = cli, 100, 100)
            }, width = 5, height = 4)
        })
        
        doProgress()

        return(ret)
    }))
}

genHTMLReportPlotsComponents <- function(fGroups, components, outPath, EICs, selfContained)
{
    cInfo <- componentInfo(components)
    isIntCl <- inherits(components, "componentsIntClust")
    
    cat("Generate component plots...\n")
    ret <- list()
    
    if (isIntCl || inherits(components, "componentsSpecClust"))
    {
        ret$dendro <- makeHTMLReportPlot("compon-dendro.svg", outPath, selfContained, {
            plot(components)
        })
    }
    
    # UNDONE: parallel option
    ret$components <- pruneList(doApply("Map", TRUE, names(components), componentTable(components), f = function(cn, ct)
    {
        if (!any(ct$group %chin% names(fGroups)))
            return(NULL)
        
        pl <- list()
        
        pl$chrom <- makeHTMLReportPlot(paste0("compon-chrom-", cn, ".svg"), outPath, selfContained, {
            mar <- par("mar")
            par(mar = c(mar[1], mar[2], 0.2, 0.2))
            plotChroms(components, cn, fGroups, retMin = TRUE, title = "", EICs = EICs) # UNDONE: params
        }, width = 6, height = 4, bg = "transparent", pointsize = 16)
        
        pl$spec <- makeHTMLReportPlot(paste0("compon-spec-", cn, ".svg"), outPath, selfContained, {
            mar <- par("mar")
            par(mar = c(mar[1], mar[2], 0.2, 0.2))
            cir <- cInfo[name == cn] 
            plotSpectrum(components, cn)
        }, width = 7, height = 4, pointsize = 16)
        
        if (isIntCl)
        {
            pl$profileRel <- makeHTMLReportPlot(paste0("compon-int_rel-", cn, ".svg"), outPath, selfContained, {
                plotInt(components, index = cn, main = "normalized")
            }, width = 7, height = 6, pointsize = 16)
            
            pl$profileAbs <- makeHTMLReportPlot(paste0("compon-int_abs-", cn, ".svg"), outPath, selfContained, {
                plotInt(fGroups[, components[[cn]]$group], average = clusterProperties(components)$average,
                        main = "absolute")
            }, width = 7, height = 6, pointsize = 16)
        }
        
        doProgress()
        
        return(pl)
    }))
    
    return(ret)
}

generateHTMLReportPlots <- function(fGroups, MSPeakLists, formulas, compounds, compsCluster, components, TPs, outPath, EICs,
                                    selfContained)
{
    ret <- list()
    
    cat("Genarate summary plots...")
    ret$overview$chroms <- makeHTMLReportPlot("chroms.svg", outPath, selfContained, {
        par(mai = c(0.9, 0.8, 0.6, 0.1))
        # UNDONE: params
        plotChroms(fGroups, rtWindow = 30, mzExpWindow = 0.005, retMin = TRUE, topMost = 1, topMostByRGroup = FALSE,
                   EICs = EICs, showPeakArea = TRUE, showFGroupRect = FALSE, colourBy = "fGroups", showLegend = FALSE,
                   onlyPresent = TRUE)
    }, width = 10, height = 4)
    
    ret$overview$retMZ <- makeHTMLReportPlot("retmz.svg", outPath, selfContained, {
        par(mai = c(0.9, 0.8, 0.1, 0.1))
        # UNDONE: params
        plot(fGroups, colourBy = "fGroups", showLegend = FALSE, retMin = TRUE)
    }, width = 10, height = 4)
    
    rGroupLenNonEmpty <- length(replicateGroups(removeEmptyAnalyses(fGroups)))
    rGroupLen <- length(replicateGroups(fGroups))
    anyOverlap <- rGroupLen > 1 &&
        length(unique(fGroups, which = replicateGroups(fGroups), outer = TRUE)) < length(fGroups)
    
    ret$overview$chord <- ret$overview$venn  <- ret$overview$UpSet <- NULL
    if (anyOverlap && rGroupLenNonEmpty > 1)
    {
        if (rGroupLenNonEmpty > 2)
        {
            ret$overview$chord <- makeHTMLReportPlot("chord.svg", outPath, selfContained, {
                # UNDONE: params(?)
                plotChord(fGroups, average = TRUE)
            }, width = 7, height = 7)
        }
        if (rGroupLen < 6)
        {
            ret$overview$venn <- makeHTMLReportPlot("venn.svg", outPath, selfContained, {
                # UNDONE: params(?)
                plotVenn(fGroups)
            }, width = 7, height = 7)
        }
        
        # UpSet
        ret$overview$UpSet <- makeHTMLReportPlot("upset.svg", outPath, selfContained, {
            # UNDONE: params(?)
            print(plotUpSet(fGroups))
        }, width = 7, height = 7)
    }
    cat(" Done!\n")

    ret$chromsLarge <- genHTMLReportPlotsChromsLarge(fGroups, outPath, EICs, selfContained)
    ret$chromsSmall <- genHTMLReportPlotsChromsSmall(fGroups, outPath, EICs, selfContained)
    ret$chromsFeatures <- genHTMLReportPlotsChromsFeatures(fGroups, outPath, EICs, selfContained)
    
    ret$structs <- genHTMLReportPlotsStructs(fGroups, compounds, outPath, selfContained)
    if (!is.null(formulas))
        ret$formulas <- genHTMLReportPlotsFormulas(formulas, MSPeakLists, outPath, selfContained)
    if (!is.null(compounds))
        ret$compounds <- genHTMLReportPlotsCompounds(compounds, MSPeakLists, formulas, outPath, selfContained)
    if (!is.null(compsCluster))
        ret$compsCluster <- genHTMLReportPlotsCompsCluster(compsCluster, outPath, selfContained)
    if (!is.null(components) && !inherits(components, "componentsTPs"))
        ret$components <- genHTMLReportPlotsComponents(fGroups, components, outPath, EICs, selfContained)
    
    return(ret)
}

reportHTMLUtils$methods(
    plotImg = function(p)
    {
        if (properties$selfContained)
            return(paste0("<img src=", knitr::image_uri(p), "></img>"))
        return(paste0("<img src='", p, "'></img>"))
    },
    genISTDGraph = function(set = NULL)
    {
        args <- list(objects$fGroups, width = "100%", height = "90%") # NOTE: bit less height to avoid vertical scrollbar
        if (!is.null(set))
            args <- c(args, set = set)
        do.call(plotGraph, args)
    },
    genTPGraphs = function()
    {
        cInfo <- componentInfo(objects$components)
        pars <- parents(objects$TPs)
        hwidgets <- lapply(seq_len(nrow(cInfo)), function(i)
        {
            DOMID <- paste0('TPGraph_', cInfo$name[i])
            TPInd <- match(cInfo$parent_name[i], pars$name, nomatch = NA)
            if (is.na(TPInd))
                return(htmltools::div(id = DOMID))
        
            # NOTE: bit less height to avoid scrollbar in card    
            gr <- plotGraph(objects$TPs, which = TPInd, components = objects$components, structuresMax = 10,
                            width = "100%", height = "97%") # UNDONE: params
            gr$elementId <- DOMID
            return(htmlwidgets::onRender(gr, htmlwidgets::JS("function(el, x)
{
    el.style.display = 'none'; // hide by default
    // auto fit TP graphs on window resize
    let ro = new ResizeObserver(() => { document.getElementById('graph' + el.id).chart.fit(); });
    ro.observe(el);
}")))
        })
        return(do.call(htmltools::tagList, hwidgets))
    },
    genCompClustsImgs = function()
    {
        elements <- unlist(Map(names(plots$compsCluster), plots$compsCluster, f = function(grp, pl)
        {
            lapply(seq_along(pl$mcs), function(i)
            {
                img <- htmltools::img(src = pl$mcs[[i]], class = paste0("mcs mcs-", grp),
                                      style = "max-height: 250px; display: none;")
            })
        }), recursive = FALSE)
        elements <- do.call(htmltools::tagList, elements)
        return(htmltools::div(style = list(display = "flex"), elements))
    },
    genIntClustHeatMap = function()
    {
        plotHeatMap(objects$components, interactive = TRUE) # UNDONE: make interactive configfurable
    },
    genComponNTGraph = function(s)
    {
        args <- list(objects$components, onlyLinked = TRUE, width = "100%", height = "100%")
        if (!is.null(s))
            args <- c(args, list(set = s))
        do.call(plotGraph, args)
    }
)
