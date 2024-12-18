# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

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
    genFGTableComponents = function()
    {
        tab <- getFGTable(objects$fGroups, ",", settings$features$retMin, settings$features$aggregateConcs,
                          settings$features$aggregateTox)
        groupDefs <- getFGGroupDefs(tab, "component", replicateGroups(objects$fGroups))
        colDefs <- getFeatGroupColDefs(tab)
        
        ctab <- as.data.table(objects$components)
        setnames(ctab, "name", "component")
        ctab <- removeDTColumnsIfPresent(ctab, c(
            # already present from features table
            "ret", "mz",
            
            # general component properties
            "mz_increment", "rt_increment", "ret_min", "ret_max", "ret_range",
            replicateGroups(objects$fGroups), # nontarget: presence in replicate groups
            "cmp_ret", "cmp_mz", "cmp_retsd", "neutral_mass", "cmp_ppm", "analysis", "size",
            
            # for graphs
            "links"
        ))
        if (all(ctab$intensity == 1))
            ctab[, intensity := NULL]
        tab <- merge(tab, ctab, by = "group", sort = FALSE)
        
        if (!is.null(tab[["set"]]))
        {
            colDefs$set <- reactable::colDef(filterInput = function(values, name) reactSelectFilter("detailsTabComponents",
                                                                                                    values, name))
        }
        
        cmpGrpCols <- setdiff(names(ctab), c("component", "group"))
        if (length(cmpGrpCols) > 0)
            groupDefs <- c(groupDefs, list(reactable::colGroup("component", cmpGrpCols, headerStyle = getFGColSepStyle())))
        
        makeFGReactable(tab, "detailsTabComponents", colDefs = colDefs, groupDefs = groupDefs, visible = FALSE,
                        plots = plots, settings = settings, objects = objects, groupBy = "component")
    },
    
    genComponentInfoTable = function()
    {
        tab <- componentInfo(objects$components)
        tab <- removeDTColumnsIfPresent(tab, c("links", "size"))
        
        for (col in intersect(c("neutral_mass", "mz_increment"), names(tab)))
        {
            # NOTE: neutral_mass may be a character, eg for CAMERA where multiple masses are collapsed
            if (is.numeric(tab[[col]]))
                set(tab, j = col, value = round(tab[[col]], 5))
        }
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

