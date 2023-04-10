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

makeHTMLReportPlotOld <- function(out, outPath, selfContained, code, ...)
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

makeHTMLReportPlot <- function(outPrefix, outPath, selfContained, func, args, parParams = list(), cache = TRUE,
                               doPrint = FALSE, ...)
{
    out <- if (cache)
    {
        hash <- do.call(paste0(func, "Hash"), args)
        out <- paste0(outPrefix, "-", hash, ".svg")
    }
    else
        out <- paste0(outPrefix, ".svg")
    ppath <- file.path(getHTMLReportPlotPath(outPath), out)
    
    if (!cache || !file.exists(ppath))
    {
        withSVGLite(ppath, standalone = TRUE, code = {
            if (length(parParams))
                do.call(par, parParams)
            if (doPrint)
                print(do.call(func, args))
            else
                do.call(func, args)
        }, ...)
    }
    
    return(getHTMLReportFinalPlotPath(out, selfContained))
}

genHTMLReportPlotsChromsLarge <- function(fGroups, outPath, EICs, EICParams, selfContained)
{
    cat("Generate large chromatograms...\n")
    # UNDONE: parallel option
    doApply("sapply", TRUE, names(fGroups), function(grp)
    {
        doProgress()
        makeHTMLReportPlot("chrom_large-", outPath, selfContained, "plotChroms",
                           list(fGroups, groupName = grp, retMin = TRUE, EICs = EICs, EICParams = EICParams,
                                colourBy = "rGroups", title = "", bty = "l"),
                           parParams = list(mar = c(4.1, 4.1, 0.2, 0.2)),
                           width = 6, height = 4, bg = "transparent", pointsize = 16)
    }, simplify = FALSE)
}

genHTMLReportPlotsChromsSmall <- function(fGroups, outPath, EICs, EICParams, selfContained)
{
    cat("Generate small chromatograms...\n")
    # UNDONE: parallel option
    doApply("sapply", TRUE, names(fGroups), function(grp)
    {
        doProgress()
        makeHTMLReportPlot("chrom_small", outPath, selfContained, "plotChroms",
                           list(fGroups, groupName = grp, retMin = TRUE, EICs = EICs,
                                EICParams = modifyList(EICParams, list(topMost = 1, topMostByRGroup = FALSE,
                                                                       onlyPresent = TRUE)),
                                showFGroupRect = FALSE, showPeakArea = TRUE, title = "", bty = "n"),
                           parParams = list(mai = c(0, 0, 0, 0), lwd = 10), width = 12, height = 4, bg = "transparent",
                           pointsize = 16)
    }, simplify = FALSE)
}

genHTMLReportPlotsChromsFeatures <- function(fGroups, outPath, EICs, EICParams, selfContained)
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
            makeHTMLReportPlot("chrom_feat", outPath, selfContained, "plotChroms",
                               list(fGroups, analysis = ana, groupName = grp, retMin = TRUE, EICs = EICs,
                                    EICParams = modifyList(EICParams, list(topMost = NULL, onlyPresent = FALSE),
                                                           keep.null = TRUE), showFGroupRect = FALSE,
                                    showPeakArea = TRUE, title = "", bty = "l"),
                               parParams = list(mar = c(4.1, 4.1, 0.2, 0.2)), width = 6, height = 4, bg = "transparent",
                               pointsize = 20, scaling = 1)
        })
    }, simplify = FALSE)
}

genHTMLReportPlotsIntPlots <- function(fGroups, outPath, selfContained)
{
    cat("Generate intensity plots...\n")
    # UNDONE: parallel option
    
    mainArgs <- list(average = TRUE, normalize = TRUE, plotArgs = list(bty = "l"))
    if (isFGSet(fGroups))
        mainArgs <- c(mainArgs, list(sets = TRUE))
    else
        mainArgs <- c(mainArgs, list(col = "black"))
    
    doApply("sapply", TRUE, names(fGroups), function(grp)
    {
        doProgress()
        makeHTMLReportPlot("int_plot", outPath, selfContained, "plotInt",  c(list(fGroups[, grp]), mainArgs),
                           parParams = list(mar = c(4.1, 4.1, 1, 0.1)), width = 8, height = 4, bg = "transparent",
                           pointsize = 16)
    }, simplify = FALSE)
}

genHTMLReportPlotsStructs <- function(fGroups, compounds, outPath, selfContained)
{
    scrStructInfo <- if (isScreening(fGroups)) screenInfo(fGroups)[, c("SMILES", "InChIKey"), with = FALSE] else NULL
    compStructInfo <- NULL
    if (!is.null(compounds) && length(compounds) != 0)
    {
        compStructInfo <- as.data.table(compounds)[, c("group", "SMILES", "InChIKey"), with = FALSE]
        compStructInfo[, index := seq_len(.N), by = "group"]
        compStructInfo <- compStructInfo[index <= 25][, -c("group", "index")] # UNDONE
    }
    
    structInfo <- rbindlist(list(scrStructInfo, compStructInfo))
    if (nrow(structInfo) > 0)
    {
        structInfo <- unique(structInfo[!is.na(SMILES)], by = "InChIKey")
        cat("Generate structures...\n")
        # UNDONE: parallel option
        return(doApply("Map", TRUE, structInfo$InChIKey, structInfo$SMILES, f = function(ik, smi)
        {
            # NOTE: we use the InChIKey here instead of makeHash()
            pf <- file.path(getHTMLReportPlotPath(outPath), paste0("struct-", ik, ".svg"))
            if (!file.exists(pf))
                saveRCDKStructure(getMoleculesFromSMILES(smi)[[1]], "svg", pf, 100, 100)
            doProgress()
            return(getHTMLReportFinalPlotPath(pf, selfContained))
        }))
    }
    return(list())
}

genHTMLReportPlotsMSPeakLists <- function(MSPeakLists, outPath, selfContained)
{
    cat("Generate MS spectra...\n")
    
    if (length(MSPeakLists) == 0)
        return(list())
    
    # UNDONE: parallel option
    return(doApply("sapply", TRUE, groupNames(MSPeakLists), function(grp)
    {
        ret <- list()
        
        args <- list(MSPeakLists, groupName = grp, title = "")
        pp <- list(mar = c(4.1, 4.1, 0.2, 0.2))
        
        ret$MS <- makeHTMLReportPlot("spec-MS", outPath, selfContained, "plotSpectrum",
                                     c(list(MSLevel = 1), args), parParams = pp,
                                     width = 7, height = 4)
        
        ret$MSMS <- if (!is.null(MSPeakLists[[grp]][["MSMS"]]))
        {
            makeHTMLReportPlot("spec-MSMS", outPath, selfContained, "plotSpectrum",
                               c(list(MSLevel = 2), args), parParams = pp,
                               width = 7, height = 4)
        }
        else
            ""
        
        doProgress()
        
        return(ret)
    }, simplify = FALSE))
}

genHTMLReportPlotsFormulas <- function(formulas, MSPeakLists, outPath, selfContained)
{
    cat("Generate formula annotation plots...\n")
    
    if (length(formulas) == 0)
        return(list())
    
    # UNDONE: parallel option
    return(doApply("Map", TRUE, groupNames(formulas), annotations(formulas), f = function(grp, ann)
    {
        ret <- list()
        
        if (nrow(ann) > 25)
            ann <- ann[seq_len(25)] # UNDONE
        
        ret$spectra <- sapply(seq_len(nrow(ann)), function(index)
        {
            if (is.null(MSPeakLists[[grp]][["MSMS"]]))
                return("")
            makeHTMLReportPlot("form-spec", outPath, selfContained, "plotSpectrum",
                               list(formulas, index, grp, MSPeakLists = MSPeakLists), width = 7, height = 4)
        })
        
        ret$scores <- sapply(seq_len(nrow(ann)), function(index)
        {
            # UNDONE: params
            makeHTMLReportPlot("form-scores", outPath, selfContained, "plotScores", list(formulas, index, grp), 
                               width = 6, height = 5)
        })
        
        doProgress()
        
        return(ret)
    }))
}

genHTMLReportPlotsCompounds <- function(compounds, MSPeakLists, formulas, outPath, selfContained)
{
    cat("Generate compound annotation plots...\n")
    
    if (length(compounds) == 0)
        return(list())
    
    # UNDONE: parallel option
    return(doApply("Map", TRUE, groupNames(compounds), annotations(compounds), f = function(grp, ann)
    {
        ret <- list()

        if (nrow(ann) > 25)
            ann <- ann[seq_len(25)] # UNDONE
        
        ret$spectra <- sapply(seq_len(nrow(ann)), function(index)
        {
            makeHTMLReportPlot("comp-spec", outPath, selfContained, "plotSpectrum",
                               list(compounds, index, grp, MSPeakLists, formulas, FALSE, title = ""),
                               parParams = list(mar = c(4.1, 4.1, 0.2, 0.2)), width = 7, height = 4, pointsize = 16)
        })
        
        ret$scores <- sapply(seq_len(nrow(ann)), function(index)
        {
            # UNDONE: params
            makeHTMLReportPlot("comp-scores", outPath, selfContained, "plotScores", list(compounds, index, grp),
                               parParams = list(mar = c(4.1, 4.1, 0.4, 0.2)), width = 7, height = 4, pointsize = 16)
        })
        
        doProgress()
        
        return(ret)
    }))
}

genHTMLReportPlotsCompsCluster <- function(compsCluster, outPath, selfContained)
{
    cat("Generate compound cluster plots...\n")
    
    if (length(compsCluster) == 0)
        return(list())
    
    # UNDONE: parallel option
    return(doApply("Map", TRUE, groupNames(compsCluster), cutClusters(compsCluster), f = function(grp, ct)
    {
        ret <- list()
        
        ret$dendro <- makeHTMLReportPlot("comp-clust-dendro", outPath, selfContained, "plot",
                                         list(compsCluster, groupName = grp), width = 12, height = 4, pointsize = 16)
        
        ret$mcs <- sapply(sort(unique(ct)), function(cli)
        {
            makeHTMLReportPlot("comp-clust-mcs", outPath, selfContained, "plotStructure",
                               list(compsCluster, groupName = grp, cluster = cli, 100, 100), width = 5, height = 4)
        })
        
        doProgress()

        return(ret)
    }))
}

genHTMLReportPlotsComponents <- function(fGroups, components, outPath, EICs, EICParams, selfContained)
{
    cat("Generate component plots...\n")
    
    if (length(components) == 0)
        return(list())

    cInfo <- componentInfo(components)
    isIntCl <- inherits(components, "componentsIntClust")
    ret <- list()
    
    if (isIntCl || inherits(components, "componentsSpecClust"))
        ret$dendro <- makeHTMLReportPlot("compon-dendro.svg", outPath, selfContained, "plot", list(components))
    
    # UNDONE: parallel option
    ret$components <- pruneList(doApply("Map", TRUE, names(components), componentTable(components), f = function(cn, ct)
    {
        if (!any(ct$group %chin% names(fGroups)))
            return(NULL)
        
        pl <- list()
        
        pl$chrom <- makeHTMLReportPlot("compon-chrom", outPath, selfContained, "plotChroms",
                                       list(components, cn, fGroups, retMin = TRUE, title = "", EICs = EICs,
                                            EICParams = EICParams),
                                       parParams = list(mar = c(4.1, 4.1, 0.2, 0.2)), width = 6, height = 4,
                                       bg = "transparent", pointsize = 16)

        pl$spec <- makeHTMLReportPlot("compon-spec", outPath, selfContained, "plotSpectrum", list(components, cn),
                                      parParams = list(mar = c(4.1, 4.1, 0.2, 0.2)), width = 7, height = 4,
                                      pointsize = 16)
        
        if (isIntCl)
        {
            pl$profileRel <- makeHTMLReportPlot("compon-int_rel", outPath, selfContained, "plotInt",
                                                list(components, index = cn,
                                                     plotArgs = list(main = "normalized", bty = "l")),
                                                width = 7, height = 6, pointsize = 16)
            
            pl$profileAbs <- makeHTMLReportPlot("compon-int_abs", outPath, selfContained, "plotInt",
                                                list(fGroups[, components[[cn]]$group],
                                                     average = clusterProperties(components)$average,
                                                     plotArgs = list(main = "absolute", bty = "l")),
                                                width = 7, height = 6, pointsize = 16)
        }
        
        doProgress()
        
        return(pl)
    }))
    
    return(ret)
}

genHTMLReportPlotsTPs <- function(fGroups, components, MSPeakLists, formulas, compounds, outPath, EICs, selfContained)
{
    if (is.null(MSPeakLists) || length(components) == 0 || length(MSPeakLists) == 0)
        return(list())

    scr <- if (isScreening(fGroups)) screenInfo(fGroups) else NULL
    
    cat("Generate TP similarity plots...\n")
    # UNDONE: parallel option
    return(doApply("Map", TRUE, names(components), componentTable(components), split(componentInfo(components), seq_len(length(components))), f = function(cmpName, cmpTab, cmpInfoRow)
    {
        ret <- mapply(split(cmpTab, seq_len(nrow(cmpTab))), seq_len(nrow(cmpTab)), FUN = function(ctRow, ctInd)
        {
            if (!all(c(ctRow$group, cmpInfoRow$parent_group) %chin% names(fGroups)))
                return("") # fGroups was probably subset
            
            # try to plot a mirror spectrum: use compounds if possible, otherwise try formulas or finally peak lists
            plSpecArgs <- list()
            
            if (isScreening(fGroups))
            {
                scrParRow <- scr[name == cmpInfoRow$parent_name & group == cmpInfoRow$parent_group]
                scrTPRow <- scr[name == ctRow$TP_name & group == ctRow$group]
                if (!is.null(compounds) && !is.null(scr[["compRank"]]) &&
                    all(c(cmpInfoRow$parent_group, ctRow$group) %chin% groupNames(compounds)) &&
                    nrow(scrTPRow) == 1 && !is.na(scrParRow$compRank) && !is.na(scrTPRow$compRank))
                {
                    plSpecArgs <- list(obj = compounds, formulas = formulas,
                                       index = c(scrParRow$compRank, scrTPRow$compRank), MSPeakLists = MSPeakLists,
                                       plotStruct = FALSE)
                }
                else if (!is.null(formulas) && !is.null(scr[["formRank"]]) &&
                         all(c(cmpInfoRow$parent_group, ctRow$group) %chin% groupNames(formulas)) &&
                         nrow(scrTPRow) == 1 && !is.na(scrParRow$formRank) && !is.na(scrTPRow$formRank) &&
                         !is.null(MSPeakLists[[cmpInfoRow$parent_group]][["MSMS"]]) &&
                         !is.null(MSPeakLists[[ctRow$group]][["MSMS"]]))
                {
                    plSpecArgs <- list(obj = formulas,
                                       index = c(scrParRow$formRank, scrTPRow$formRank),
                                       MSPeakLists = MSPeakLists)
                }
            }
            
            if (length(plSpecArgs) == 0 && !is.null(MSPeakLists[[cmpInfoRow$parent_group]][["MSMS"]]) &&
                !is.null(MSPeakLists[[ctRow$group]][["MSMS"]]))
            {
                # no formulas/compounds, try peak lists
                plSpecArgs <- list(obj = MSPeakLists, MSLevel = 2)
            }
            
            if (length(plSpecArgs) == 0)
                return("")
            
            # UNDONE
            specSimParams <- getDefSpecSimParams()
            return(makeHTMLReportPlot("spec_sim", outPath, selfContained, "plotSpectrum",
                                      c(plSpecArgs, list(groupName = c(cmpInfoRow$parent_group, ctRow$group),
                                                         specSimParams = specSimParams, title = "")),
                                      parParams = list(mar = c(mar[1], mar[2], 0.2, 0.2)), width = 10, height = 5,
                                      pointsize = 16))
        })
        
        doProgress()
        
        return(ret)
    }))
}

generateHTMLReportPlots <- function(fGroups, MSPeakLists, formulas, compounds, compsCluster, components, TPs, outPath,
                                    EICs, EICParams, selfContained)
{
    ret <- list()
    
    cat("Genarate summary plots...")
    
    # UNDONE: params
    
    ret$overview$chroms <- makeHTMLReportPlot("chroms", outPath, selfContained, "plotChroms",
                                              list(fGroups, retMin = TRUE, EICs = EICs,
                                                   EICParams = modifyList(EICParams, list(topMost = 1,
                                                                                          topMostByRGroup = FALSE,
                                                                                          onlyPresent = TRUE)),
                                                   showPeakArea = TRUE, showFGroupRect = FALSE, colourBy = "fGroups",
                                                   showLegend = FALSE),
                                              parParams = list(mai = c(0.9, 0.8, 0.6, 0.1)), width = 10, height = 4)
    
    ret$overview$retMZ <- makeHTMLReportPlot("retmz", outPath, selfContained, "plot",
                                             list(fGroups, colourBy = "fGroups", showLegend = FALSE, retMin = TRUE),
                                             parParams = list(mai = c(0.9, 0.8, 0.1, 0.1)), width = 10, height = 4)
    
    rGroupLenNonEmpty <- length(replicateGroups(removeEmptyAnalyses(fGroups)))
    rGroupLen <- length(replicateGroups(fGroups))
    anyOverlap <- rGroupLen > 1 &&
        length(unique(fGroups, which = replicateGroups(fGroups), outer = TRUE)) < length(fGroups)
    
    ret$overview$chord <- ret$overview$venn  <- ret$overview$UpSet <- NULL
    if (anyOverlap && rGroupLenNonEmpty > 1)
    {
        if (rGroupLenNonEmpty > 2)
        {
            ret$overview$chord <- makeHTMLReportPlot("chord", outPath, selfContained, "plotChord",
                                                     list(fGroups, average = TRUE), width = 7, height = 7)
        }
        if (rGroupLen < 6)
        {
            ret$overview$venn <- makeHTMLReportPlot("venn", outPath, selfContained, "plotVenn", list(fGroups),
                                                    width = 7, height = 7)
        }
        
        # UpSet
        ret$overview$UpSet <- makeHTMLReportPlot("upset", outPath, selfContained, "plotUpSet", list(fGroups),
                                                 doPrint = TRUE, width = 7, height = 7)
    }
    cat(" Done!\n")

    ret$chromsLarge <- genHTMLReportPlotsChromsLarge(fGroups, outPath, EICs, EICParams, selfContained)
    ret$chromsSmall <- genHTMLReportPlotsChromsSmall(fGroups, outPath, EICs, EICParams, selfContained)
    ret$chromsFeatures <- genHTMLReportPlotsChromsFeatures(fGroups, outPath, EICs, EICParams, selfContained)
    
    ret$intPlots <- genHTMLReportPlotsIntPlots(fGroups, outPath, selfContained)

    gNames <- names(fGroups)
    if (!is.null(compounds))
        compounds <- compounds[gNames]
    
    ret$structs <- genHTMLReportPlotsStructs(fGroups, compounds, outPath, selfContained)

    if (!is.null(MSPeakLists))
        ret$MSPeakLists <- genHTMLReportPlotsMSPeakLists(MSPeakLists[, gNames], outPath, selfContained)    
    if (!is.null(formulas))
        ret$formulas <- genHTMLReportPlotsFormulas(formulas[gNames], MSPeakLists, outPath, selfContained)
    if (!is.null(compounds))
        ret$compounds <- genHTMLReportPlotsCompounds(compounds, MSPeakLists, formulas, outPath, selfContained)
    if (!is.null(compsCluster))
        ret$compsCluster <- genHTMLReportPlotsCompsCluster(compsCluster[gNames], outPath, selfContained)
    if (!is.null(components))
    {
        if (!inherits(components, "componentsTPs"))
            ret$components <- genHTMLReportPlotsComponents(fGroups, components, outPath, EICs, EICParams, selfContained)
        else
            ret$TPs <- genHTMLReportPlotsTPs(fGroups, components, MSPeakLists, formulas, compounds, outPath, EICs,
                                             selfContained)
    }
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
        plotHeatMap(objects$components, interactive = TRUE) # UNDONE: make interactive configurable?
    },
    genComponNTGraph = function(s)
    {
        args <- list(objects$components, onlyLinked = TRUE, width = "100%", height = "100%")
        if (!is.null(s))
            args <- c(args, list(set = s))
        do.call(plotGraph, args)
    }
)
