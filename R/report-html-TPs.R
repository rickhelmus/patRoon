# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
#' @include report-html.R
NULL

genHTMLReportPlotsTPs <- function(fGroups, components, MSPeakLists, formulas, compounds, settings, specSimParams,
                                  outPath, EICs, parallel)
{
    if (is.null(MSPeakLists) || length(components) == 0 || length(MSPeakLists) == 0)
        return(list())
    
    scr <- if (isScreening(fGroups)) screenInfo(fGroups) else NULL
    
    cat("Generate TP similarity plots...\n")
    return(doApply("Map", parallel, names(components), componentTable(components), split(componentInfo(components), seq_len(length(components))), f = function(cmpName, cmpTab, cmpInfoRow)
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
                
                if (nrow(scrParRow) > 0 && nrow(scrTPRow) > 0) # check if hit wasn't filtered away
                {
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
            }
            
            if (length(plSpecArgs) == 0 && !is.null(MSPeakLists[[cmpInfoRow$parent_group]][["MSMS"]]) &&
                !is.null(MSPeakLists[[ctRow$group]][["MSMS"]]))
            {
                # no formulas/compounds, try peak lists
                plSpecArgs <- list(obj = MSPeakLists, MSLevel = 2)
            }
            
            if (length(plSpecArgs) == 0)
                return("")
            
            return(makeHTMLReportPlot("spec_sim", outPath, "plotSpectrum",
                                      c(plSpecArgs, list(groupName = c(cmpInfoRow$parent_group, ctRow$group),
                                                         specSimParams = specSimParams, title = "")),
                                      parParams = list(mar = c(4.1, 4.1, 0.2, 0.2)), width = 10, height = 5,
                                      pointsize = 16))
        })
        
        doProgress()
        
        return(ret)
    }))
}

reportHTMLUtils$methods(
    genFGTableTPs = function()
    {
        fromTPs <- objects$components@fromTPs
        
        tabTPsFeat <- getFGTable(objects$fGroups, if (fromTPs) NULL else ",", settings$features$retMin,
                                 settings$features$aggregateConcs, settings$features$aggregateTox)
        
        tabCompon <- as.data.table(objects$components)
        tabCompon <- tabCompon[parent_group %chin% names(objects$fGroups)]
        tabCompon <- subsetDTColumnsIfPresent(tabCompon, c("name", "parent_name", "parent_group", "group", "TP_retDir",
                                                           "TP_name", "retDir", "retDiff", "mzDiff", "formulaDiff",
                                                           "specSimilarity", "mergedBy"))
        tabCompon[, cmpIndex := seq_len(.N), by = "name"]
        
        if (fromTPs)
        {
            tabTPs <- merge(tabCompon, tabTPsFeat, by.x = c("group", "TP_name"), by.y = c("group", "susp_name"),
                            sort = FALSE)
            setnames(tabTPs, c("name", "TP_name", "parent_name"),
                     c("component", "susp_name", "parent_susp_name"), skip_absent = TRUE)
        }
        else
        {
            tabTPs <- merge(tabCompon, tabTPsFeat, by = "group", sort = FALSE)
            setnames(tabTPs, "name", "component")
            tabTPs <- removeDTColumnsIfPresent(tabTPs, "susp_name")
        }
        
        parAggr <- function(parentCol) htmlwidgets::JS(sprintf("function(values, rows, groupRows)
{
    const v = rows[0]['%s'];
    return (v == null) ? '' : '<i>' + v + '</i>';
}", parentCol))
        
        parAggred <- function(level) htmlwidgets::JS(sprintf("function(cellInfo, state)
{
    return (cellInfo.level === %d) ? cellInfo.value : '';
}", level))
        
        # NOTE: below values may be in components but then from suspect list
        tabTPs[, c("parent_ret", "parent_mz") := groupInfo(objects$fGroups)[parent_group, ]]
        
        for (col in intersect(c("parent_ret", "retDiff"), names(tabTPs)))
        {
            if (settings$features$retMin)
                tabTPs[, (col) := get(col) / 60]
            tabTPs[, (col) := round(get(col), 2)]
        }
        if (!is.null(tabTPs[["specSimilarity"]]))
            tabTPs[, specSimilarity := round(specSimilarity, 2)]
        for (col in c("parent_mz", "mzDiff"))
            tabTPs[, (col) := round(get(col), 5)]
        
        # add parent intensities & screening info
        rgs <- replicateGroups(objects$fGroups)
        tabTPsPar <- unique(subsetDTColumnsIfPresent(tabTPsFeat,
                                                     c("group", rgs,
                                                       paste0("susp_", c("estIDLevel", "d_rt", "d_mz", "sets",
                                                                         "InChIKey")))),
                            by = "group")
        setnames(tabTPsPar, paste0("parent_", names(tabTPsPar)))
        tabTPs <- merge(tabTPs, tabTPsPar, by = "parent_group", sort = FALSE, all.x = TRUE)
        
        groupBy <- if (fromTPs) c("component", "susp_name") else "component"
        groupDefs <- getFGGroupDefs(tabTPs, groupBy, rgs)
        # squeeze in TP column
        groupDefs <- c(groupDefs[1:2],
                       list(reactable::colGroup("TP", columns = intersect(c("TP_name", "retDiff", "mzDiff",
                                                                            "formulaDiff", "retDir",
                                                                            "specSimilarity", "mergedBy"),
                                                                          names(tabTPs)),
                                                headerStyle = getFGColSepStyle())),
                       groupDefs[seq(3, length(groupDefs))])
        
        colDefs <- getFeatGroupColDefs(tabTPs)
        
        # set parent 'aggregates': actual value of parent feature group
        for (col in grep("^parent_", names(tabTPs), value = TRUE))
        {
            colTP <- sub("^parent_", "", col)
            colDefs[[colTP]]$aggregate <- parAggr(col)
            colDefs[[colTP]]$aggregated <- parAggred(0)
            colDefs[[colTP]]$html <- TRUE
            colDefs[[col]] <- reactable::colDef(show = FALSE)
        }
        # similar for TP retDir
        colDefs$retDir <- reactable::colDef(aggregate = parAggr("TP_retDir"), aggregated = parAggred(1), html = TRUE)
        colDefs$TP_retDir <- reactable::colDef(show = FALSE)
        
        # these are grouped
        if (!is.null(tabTPs[["TP_name"]]))
            colDefs$TP_name <- reactable::colDef("name")
        colDefs$retDiff <- reactable::colDef(name = "\U0394 ret")
        colDefs$mzDiff <- reactable::colDef(name = "\U0394 mz")
        if (!is.null(tabTPs[["formulaDiff"]]))
            colDefs$formulaDiff <- reactable::colDef(name = "\U0394 formula",
                                                     cell = function(value) htmltools::span(dangerouslySetInnerHTML = list("__html" = subscriptFormulaHTML(value, charges = FALSE))))
        
        # InChIKeys are only there for internal usage
        if (!is.null(tabTPs[["parent_susp_InChIKey"]]))
            colDefs$parent_susp_InChIKey <- reactable::colDef(show = FALSE)
        # same for cmpIndex
        colDefs$cmpIndex <- reactable::colDef(show = FALSE)
        
        makeFGReactable(tabTPs, "detailsTabTPs", FALSE, plots, settings = settings, objects = objects,
                        groupBy = groupBy, colDefs = colDefs, groupDefs = groupDefs)
    },

    genTPSimTable = function()
    {
        tab <- as.data.table(objects$components)
        tab[, cmpID := paste0(name[1], "-", seq_len(.N)), by = "name"] # make unique IDs
        tab[, name := NULL]
        cols <- data.table::copy(names(tab))
        
        for (col in cols)
        {
            if (grepl("^specSimilarity", col))
                set(tab, j = col, value = round(tab[[col]], 2))
            else if (col %chin% c("fragmentMatches", "neutralLossMatches"))
                set(tab, j = col, value = round(tab[[col]], 0))
            else if (col != "cmpID")
                set(tab, j = col, value = NULL)
        }
        
        ptab <- makePropTab(tab, if (isFGSet(objects$fGroups)) sets(objects$fGroups) else NULL, "cmpID")
        makePropReactable(ptab, "similarityTab", "cmpID", minPropWidth = 150, minValWidth = 85)
    },
    
    genTPGraphs = function()
    {
        cInfo <- componentInfo(objects$components)
        pars <- parents(objects$TPs)
        
        # NOTE: bit less height to avoid scrollbar in card
        mainArgs <- list(objects$TPs, components = objects$components, width = "100%", height = "97%")
        if (inherits(objects$TPs, "transformationProductsStructure"))
            mainArgs <- c(mainArgs, list(structuresMax = settings$TPs$graphStructuresMax))
        
        hwidgets <- lapply(seq_len(nrow(cInfo)), function(i)
        {
            DOMID <- paste0('TPGraph_', cInfo$name[i])
            TPInd <- match(cInfo$parent_name[i], pars$name, nomatch = NA)
            if (is.na(TPInd))
                return(htmltools::div(id = DOMID))
            gr <- do.call(plotGraph, c(mainArgs, list(which = TPInd)))
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

