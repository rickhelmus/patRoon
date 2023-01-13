#' @include main.R
#' @include report-html.R
NULL

makeHTMLReportPlot <- function(out, outPath, selfContained, code, ...)
{
    if (FALSE)
    {
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
    else if (selfContained)
    {
        # UNDONE: while embedding the SVG directly would be nice, this seems to give major headaches with scaling,
        # especially with Firefox... For now just base64 it :(
        withSVGLite(out, standalone = TRUE, code = code, ...)
        return(knitr::image_uri(out))
    }
    
    destPath <- file.path(outPath, "report_files", "plots")
    mkdirp(destPath)
    withSVGLite(file.path(destPath, out), standalone = TRUE, code = code, ...)
    return(file.path("report_files", "plots", out))
    
    # UNDONE: object tag makes text selectable but messes up layout...
    # return(paste0("<object data='", out, "' type='image/svg+xml' width=500 height=300></object>"))
}

genHTMLReportPlotsChroms <- function(fGroups, outPath, EICs, selfContained)
{
    sapply(names(fGroups), function(grp)
    {
        makeHTMLReportPlot(paste0("chrom-", grp, ".svg"), outPath, selfContained, {
            # UNDONE: params
            mar <- par("mar")
            par(mar = c(mar[1], mar[2], 0.2, 0.2))
            plotChroms(fGroups[, grp], 30, 0.005, TRUE, 1, TRUE, EICs, colourBy = "rGroups", title = "",
                       onlyPresent = TRUE, bty = "l")
        }, width = 6, height = 4, bg = "transparent", pointsize = 16)
    }, simplify = FALSE)
}

genHTMLReportPlotsStructs <- function(fGroups, compounds, outPath, selfContained)
{
    scrStructInfo <- if (isScreening(fGroups)) screenInfo(fGroups)[, c("SMILES", "InChIKey"), with = FALSE] else NULL
    compStructInfo <- if (!is.null(compounds)) as.data.table(compounds)[, c("SMILES", "InChIKey"), with = FALSE] else NULL
    
    structInfo <- rbindlist(list(scrStructInfo, compStructInfo))
    if (nrow(structInfo) > 0)
    {
        structInfo <- unique(structInfo, by = "InChIKey")
        return(Map(structInfo$InChIKey, structInfo$SMILES, f = function(ik, smi)
        {
            makeHTMLReportPlot(paste0("struct-", ik, ".svg"), outPath, selfContained, {
                mol <- getMoleculesFromSMILES(smi, emptyIfFails = TRUE)[[1]]
                withr::with_par(list(mar = rep(0, 4)), plot(getRCDKStructurePlot(mol, 150, 150)))
            }, width = 5, height = 4)
        }))
    }
    return(list())
}

genHTMLReportPlotsFormulas <- function(formulas, MSPeakLists, outPath, selfContained)
{
    return(Map(groupNames(formulas), annotations(formulas), f = function(grp, ann)
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
        
        return(ret)
    }))
}

genHTMLReportPlotsCompounds <- function(compounds, MSPeakLists, formulas, outPath, selfContained)
{
    return(Map(groupNames(compounds), annotations(compounds), f = function(grp, ann)
    {
        ret <- list()
        
        ret$spectra <- sapply(seq_len(nrow(ann)), function(index)
        {
            makeHTMLReportPlot(sprintf("comp-spec-%s-%d.svg", grp, index), outPath, selfContained, {
                # par(cex = 1.5) # , cex.lab = 1.5, cex.axis = 1.5
                plotSpectrum(compounds, index, grp, MSPeakLists, formulas, FALSE)
            }, width = 7, height = 4, pointsize = 14)
        })
        
        ret$scores <- sapply(seq_len(nrow(ann)), function(index)
        {
            makeHTMLReportPlot(sprintf("comp-scores-%s-%d.svg", grp, index), outPath, selfContained, {
                plotScores(compounds, index, grp) # UNDONE: params
            }, width = 6, height = 5, pointsize = 14)
        })
        
        return(ret)
    }))
}

generateHTMLReportPlots <- function(fGroups, MSPeakLists, formulas, compounds, components, TPs, outPath, EICs,
                                    selfContained)
{
    ret <- list()
    
    ret$overview$chroms <- makeHTMLReportPlot("chroms.svg", outPath, selfContained, {
        par(mai = c(0.9, 0.8, 0.6, 0.1))
        # UNDONE: params
        plotChroms(fGroups, 30, 0.005, TRUE, 1, FALSE, EICs, TRUE, FALSE, colourBy = "fGroups", showLegend = FALSE,
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
    
    ret$chroms <- genHTMLReportPlotsChroms(fGroups, outPath, EICs, selfContained)
    ret$structs <- genHTMLReportPlotsStructs(fGroups, compounds, outPath, selfContained)
    if (!is.null(formulas))
        ret$formulas <- genHTMLReportPlotsFormulas(formulas, MSPeakLists, outPath, selfContained)
    if (!is.null(compounds))
        ret$compounds <- genHTMLReportPlotsCompounds(compounds, MSPeakLists, formulas, outPath, selfContained)
    
    return(ret)
}

reportHTMLGenerator$methods(
    plotImg = function(p)
    {
        if (properties$selfContained)
            return(paste0("<img src=", knitr::image_uri(p), "></img>"))
        return(paste0("<img src='", p, "'></img>"))
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
    }
)
