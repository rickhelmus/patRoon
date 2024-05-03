#' @include main.R
#' @include report-html.R
NULL

genHTMLReportPlotsComponents <- function(fGroups, components, settings, outPath, EICs, EICParams, parallel)
{
    cat("Generate component plots...\n")
    
    if (length(components) == 0)
        return(list())
    
    isIntCl <- inherits(components, "componentsIntClust")
    ret <- list()
    
    if (isIntCl || inherits(components, "componentsSpecClust"))
        ret$dendro <- makeHTMLReportPlot("compon-dendro.svg", outPath, "plot", list(components))
    
    ret$components <- pruneList(doApply("Map", parallel, names(components), componentTable(components), f = function(cn, ct)
    {
        if (!any(ct$group %chin% names(fGroups)))
            return(NULL)
        
        pl <- list()
        
        pl$chrom <- makeHTMLReportPlot("compon-chrom", outPath, "plotChroms",
                                       list(components, cn, fGroups, retMin = settings$features$retMin, title = "",
                                            intMax = settings$features$chromatograms$intMax, EICs = EICs,
                                            EICParams = EICParams),
                                       parParams = list(mar = c(4.1, 4.1, 0.2, 0.2)), width = 6, height = 4,
                                       bg = "transparent", pointsize = 16)
        
        pl$spec <- makeHTMLReportPlot("compon-spec", outPath, "plotSpectrum", list(components, cn),
                                      parParams = list(mar = c(4.1, 4.1, 0.2, 0.2)), width = 7, height = 4,
                                      pointsize = 16)
        
        if (isIntCl)
        {
            pl$profileRel <- makeHTMLReportPlot("compon-int_rel", outPath, "plotInt",
                                                list(components, index = cn,
                                                     plotArgs = list(main = "normalized", bty = "l")),
                                                width = 7, height = 6, pointsize = 16)
            
            pl$profileAbs <- makeHTMLReportPlot("compon-int_abs", outPath, "plotInt",
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


reportHTMLUtils$methods(
    
    getComponObj = function()
    {
        # remove any data from missing fGroups (i.e. those removed after creating the components)
        return(objects$components[, names(objects$fGroups)])
    },
    
    genMainTableComponents = function()
    {
        ctab <- as.data.table(getComponObj())
        setnames(ctab, "name", "cmpName")
        # HACK: if all intensities are one than these are dummy values (eg with intclust components)
        if (all(ctab$intensity == 1))
            ctab[, intensity := NULL]
        
        ftab <- getFGReactTab(objects, settings)
        
        tab <- merge(ctab, ftab, by = "group", sort = FALSE)
        
        makeMainResultsReactableNew(tab, "Components", settings$features$retMin, plots, initView = "Components",
                                    initTabFunc = "initTabComponents")
    },
    
    genComponentInfoTable = function()
    {
        tab <- componentInfo(getComponObj())
        tab <- removeDTColumnsIfPresent(tab, c("links", "size"))
        
        for (col in intersect(c("neutral_mass", "mz_increment"), names(tab)))
            set(tab, j = col, value = round(tab[[col]], 5))
        for (col in intersect(c("cmp_ret", "cmp_retsd", "cmp_ppm", "ret_increment", "ret_min", "ret_max", "ret_range"),
                              names(tab)))
        {
            if (col != "cmp_ppm" && settings$features$retMin)
                set(tab, j = col, value = tab[[col]] / 60)
            set(tab, j = col, value = round(tab[[col]], 2))
        }
        
        ptab <- makePropTab(tab, NULL, "name")
        makePropReactable(ptab, "componentInfoTab", "name", minPropWidth = 120, minValWidth = 150)
    },
    
    genComponentsSidebar = function()
    {
        bd <- bslib::card_body(
            padding = 0,
            pruneUI(bslib::accordion,
                    maybeInclUI(settings$features$chromatograms$large, bslib::accordion_panel(
                        "Chromatogram",
                        bslib::card_body_fill(htmltools::img(id = "chrom_view-component"))
                    )),
                    bslib::accordion_panel(
                        "Spectrum",
                        bslib::card_body_fill(htmltools::img(id = "spectrum_view-component"))
                    ),
                    maybeInclUI(hasComponentsIntClust(), bslib::accordion_panel(
                        "Profile",
                        bslib::card_body_fill(
                            htmltools::img(id = "profileRel_view-component"),
                            htmltools::img(id = "profileAbs_view-component")
                        )
                    )),
                    maybeInclUI(hasComponentInfo(), bslib::accordion_panel(
                        "Component info",
                        bsCardBodyNoFill(
                            genComponentInfoTable()
                        )
                    ))
            )
        )
        cmpNames <- names(getComponObj())
        makeSideBar("Components", "Component", "Compon-select", "updateCompon", cmpNames, cmpNames,
                    "advanceCompon", bd)
    },
    
    genDetailsComponentsUI = function()
    {
        list(
            genComponentsSidebar(),
            makeFGTableCard(genMainTableComponents(), "detailsMainTableOnly", "Components")
        )
    },
    
    genIntClustHeatMap = function()
    {
        plotHeatMap(objects$components, interactive = TRUE) # UNDONE: make interactive configurable?
    },
    
    genComponNTGraph = function(s)
    {
        if (!is.null(s) && length(objects$components[, sets = s]) == 0)
            return(htmltools::div("No components for this set.")) # NOTE: plotGraph() will throw an error if no results
        
        args <- list(objects$components, onlyLinked = TRUE, width = "100%", height = "100%")
        if (!is.null(s))
            args <- c(args, list(set = s))
        do.call(plotGraph, args)
    }
)

