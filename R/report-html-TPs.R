#' @include main.R
#' @include report-html.R
NULL

genHTMLReportPlotsTPs <- function(fGroups, components, MSPeakLists, formulas, compounds, settings, specSimParams,
                                  outPath, EICs, parallel)
{
    if (is.null(MSPeakLists) || length(components) == 0 || length(MSPeakLists) == 0)
        return(list())
    
    scr <- if (isScreening(fGroups)) screenInfo(fGroups) else NULL
    
    cmpTabSplit <- split(as.data.table(components, candidates = TRUE), by = "name")
    cInfoSplit <- split(componentInfo(components), seq_len(length(components)))

    cat("Generate TP similarity plots...\n")
    return(doApply("Map", parallel, names(components), cmpTabSplit, cInfoSplit, f = function(cmpName, cmpTab, cmpInfoRow)
    {
        ret <- mapply(split(cmpTab, seq_len(nrow(cmpTab))), seq_len(nrow(cmpTab)), FUN = function(ctRow, ctInd)
        {
            if (!all(c(ctRow$group, cmpInfoRow$parent_group) %chin% names(fGroups)))
                return("") # fGroups was probably subset
            
            # try to plot a mirror spectrum: use compounds if possible, otherwise try formulas or finally peak lists
            plSpecArgs <- list()
            if (!is.null(compounds) && !is.null(cmpInfoRow[["parent_InChIKey"]]) && !is.null(ctRow[["InChIKey"]]) &&
                all(c(cmpInfoRow$parent_group, ctRow$group) %chin% groupNames(compounds)))
            {
                pwh <- match(getIKBlock1(cmpInfoRow$parent_InChIKey), compounds[[cmpInfoRow$parent_group]]$UID)
                tpwh <- match(getIKBlock1(ctRow$InChIKey), compounds[[ctRow$group]]$UID)
                if (!is.na(pwh) && !is.na(tpwh))
                {
                    plSpecArgs <- list(obj = compounds, formulas = formulas, index = c(pwh, tpwh),
                                       MSPeakLists = MSPeakLists, plotStruct = FALSE)
                }
            }
            
            if (length(plSpecArgs) == 0 && !is.null(formulas) && !is.null(cmpInfoRow[["parent_formula"]]) &&
                !is.null(ctRow[["formula"]]) &&
                all(c(cmpInfoRow$parent_group, ctRow$group) %chin% groupNames(formulas)) &&
                !is.null(MSPeakLists[[cmpInfoRow$parent_group]][["MSMS"]]) &&
                !is.null(MSPeakLists[[ctRow$group]][["MSMS"]]))
            {
                pwh <- match(cmpInfoRow$parent_formula, formulas[[cmpInfoRow$parent_group]]$UID)
                tpwh <- match(ctRow$formula, formulas[[ctRow$group]]$UID)
                if (!is.na(pwh) && !is.na(tpwh))
                    plSpecArgs <- list(obj = formulas, index = c(pwh, tpwh), MSPeakLists = MSPeakLists)
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

getTPCandReactTab <- function(obj, plots)
{
    ret <- rbindlist(lapply(componentTable(obj), function(cmp)
    {
        allc <- rbindlist(cmp$candidates, fill = TRUE, idcol = "groupInd")
        allc[, group := cmp$group[groupInd]][, groupInd := NULL]
        return(allc)
    }), fill = TRUE, idcol = "cmpName")
    setnames(ret, "name", "candidate_name")
    if (!is.null(ret[["InChIKey"]]))
        ret[, TP_structure := plots$structs[InChIKey]]
    return(ret)
}

reportHTMLUtils$methods(
    getTPComponObj = function()
    {
        # remove any data from missing fGroups (i.e. those removed after creating the components)
        gNames <- names(objects$fGroups)
        ret <- objects$components[, gNames]
        ret <- ret[componentInfo(ret)$parent_group %chin% gNames]
        return(ret)
    },
    
    getTPComponIDs = function() names(getTPComponObj()),

    getTPComponNames = function()
    {
        ret <- getTPComponIDs()
        if (parentsFromScreening(objects$components))
            ret <- paste0(ret, " (", componentInfo(getTPComponObj())$parent_name, ")")
        return(ret)
    },

    getTPParentInfoJSON = function()
    {
        cInfo <- componentInfo(getTPComponObj())
        cInfo <- subsetDTColumnsIfPresent(cInfo, c("name", "parent_name", "parent_group", "parent_InChIKey"))
        setnames(cInfo, c("name", "parent_name", "parent_group", "parent_InChIKey"),
                 c("cmpName", "name", "group", "InChIKey"),
                 skip_absent = TRUE)
        splt <- lapply(split(cInfo, by = "cmpName", keep.by = FALSE), as.list)
        return(jsonlite::toJSON(splt, auto_unbox = TRUE))
    },
    
    genMainTableTPsParents = function()
    {
        compObj <- getTPComponObj()
        cInfo <- componentInfo(compObj)

        tab <- if (parentsFromScreening(objects$components))
        {
            ft <- getFGReactTab(objects, settings, collapseSuspects = NULL, onlyHits = TRUE)
            merge(cInfo, ft, by.x = c("parent_group", "parent_name"), by.y = c("group", "susp_name"))
        }
        else
        {
            ft <- getFGReactTab(objects, settings)
            merge(cInfo, tab, by.x = "parent_group", by.y = "group")
        }
        
        setnames(tab, c("name", "parent_group"), c("component", "group"))
        tab <- tab[match(cInfo$name, component, nomatch = 0)] # restore component order

        # used for table row selection by sidebar
        tab[, cmpName := component]
        
        if (!is.null(tab[["parent_InChIKey"]]))
            tab[, parent_structure := plots$structs[parent_InChIKey]]
            
        makeMainResultsReactable(tab, "TPsParents", settings$features$retMin, plots, initView = "TPsParents",
                                 initTabFunc = "initTabTPsParents")
    },
    
    genMainTableTPsByGroup = function()
    {
        ctab <- as.data.table(getTPComponObj())
        setnames(ctab, "name", "cmpName")
        ctab[, cmpIndex := seq_len(.N), by = "cmpName"]
        ftab <- getFGReactTab(objects, settings)
        tab <- merge(ctab, ftab, by = "group", sort = FALSE)
        makeMainResultsReactable(tab, "TPsByGroup", settings$features$retMin, plots, initView = "TPsByGroup",
                                 initTabFunc = "initTabTPs")
    },
    genMainTableTPsCandSuspect = function()
    {
        candTab <- getTPCandReactTab(getTPComponObj(), plots)
        
        if (TPsFromScreening(objects$components))
        {
            scr <- getFGScreeningReactTab(screenInfo(objects$fGroups), plots)
            candTab <- merge(candTab, scr, by.x = c("group", "candidate_name"), by.y = c("group", "susp_name"),
                             sort = FALSE)
        }
        
        makeMainResultsReactable(candTab, "TPsCandSuspect", settings$features$retMin, plots)
    },

    genMainTableTPsBySuspect = function()
    {
        
        candTab <- getTPCandReactTab(getTPComponObj(), plots)
        
        if (TPsFromScreening(objects$components))
        {
            scr <- getFGScreeningReactTab(screenInfo(objects$fGroups), plots)
            candTab <- merge(candTab, scr, by.x = c("group", "candidate_name"), by.y = c("group", "susp_name"),
                             sort = FALSE)
        }
        
        candTab[, candidate_groups := paste0(group, collapse = ", "), by = "candidate_name"][, group := NULL]
        candTab <- unique(candTab, by = "candidate_name")
        makeMainResultsReactable(candTab, "TPsBySuspect", settings$features$retMin, plots, initView = "TPsBySuspect",
                                 initTabFunc = "initTabTPs")
    },
    genMainTableTPsCandGroup = function()
    {
        ctab <- as.data.table(getTPComponObj(), candidates = TRUE)
        setnames(ctab, "name", "cmpName")
        ctab[, cmpIndex := seq_len(.N), by = "cmpName"]

        ftab <- if (TPsFromScreening(objects$components))
            getFGReactTab(objects, settings, collapseSuspects = NULL, onlyHits = TRUE)
        else
            getFGReactTab(objects, settings)
        
        ftab <- getFGReactTab(objects, settings,
                              collapseSuspects = if (TPsFromScreening(objects$components)) NULL else ",",
                              onlyHits = TRUE)
        
        # remove overlapping columns (eg er, mz)
        ctab <- removeDTColumnsIfPresent(ctab, setdiff(names(ftab), "group"))
        
        tab <- if (TPsFromScreening(objects$components))
            merge(ctab, ftab, by.x = c("group", "candidate_name"), by.y = c("group", "susp_name"), sort = FALSE)
        else
            merge(ctab, ftab, by = "group", sort = FALSE)
        
        # HACK: use a different name (and col definition) so that we get a hidden column used for filtering
        setnames(tab, "candidate_name", "candidate_ID")
        makeMainResultsReactable(tab, "TPsCandGroup", settings$features$retMin, plots)
    },
    
    genTPSimTable = function()
    {
        tab <- as.data.table(getTPComponObj())
        tab[, cmpID := paste0(name, "-", group)] # make unique IDs
        tab[, name := NULL]
        cols <- data.table::copy(names(tab))
        
        for (col in cols)
        {
            if (grepl("^specSimilarity", col))
                set(tab, j = col, value = round(tab[[col]], 2))
            else if (col %chin% c("totalFragmentMatches", "totalNeutralLossMatches"))
                set(tab, j = col, value = round(tab[[col]], 0))
            else if (col != "cmpID")
                set(tab, j = col, value = NULL)
        }
        
        ptab <- makePropTab(tab, if (isFGSet(objects$fGroups)) sets(objects$fGroups) else NULL, "cmpID")
        makePropReactable(ptab, "similarityTab", "cmpID", minPropWidth = 150, minValWidth = 85)
    },
    
    genTPGraphs = function()
    {
        compObj <- getTPComponObj()
        cInfo <- componentInfo(compObj)
        pars <- parents(objects$TPs)
        
        # NOTE: bit less height to avoid scrollbar in card
        mainArgs <- list(objects$TPs, components = compObj, width = "100%", height = "85%")
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
    },
    
    genTPsSidebar = function()
    {
        hasParStruct <- "parent_InChIKey" %in% names(componentInfo(objects$components))
        bd <- bslib::card_body(
            padding = 0,
            pruneUI(bslib::accordion,
                    maybeInclUI(settings$features$chromatograms$large, bslib::accordion_panel(
                        "Chromatogram",
                        bslib::card_body_fill(htmltools::img(id = "chrom_view-parent"))
                    )),
                    maybeInclUI(hasParStruct, bslib::accordion_panel(
                        "Structure",
                        bsCardBodyNoFill(
                            htmltools::img(id = "struct_view-parent", style = "min-width: 20%;")
                        )
                    )),
                    maybeInclUI(parentsFromScreening(objects$components), bslib::accordion_panel(
                        "Suspect information",
                        bsCardBodyNoFill(
                            genSuspInfoTable("parentSuspInfoTab"),
                        )
                    )),
                    maybeInclUI(parentsFromScreening(objects$components) && hasSuspAnn(), bslib::accordion_panel(
                        "Suspect annotation",
                        bsCardBodyNoFill(
                            genSuspAnnTable("parentSuspAnnTab")
                        )
                    )),
                    maybeInclUI(settings$features$intensityPlots, bslib::accordion_panel(
                        "Intensities",
                        class = "mt-2",
                        bslib::card_body_fill(htmltools::img(id = "int_plot-parent"))
                    )),
                    maybeInclUI(hasTPGraphs() && hasComponentsFromTPs(), bslib::accordion_panel(
                        "Transformations",
                        bslib::card(
                            full_screen = TRUE,
                            style = "margin-bottom: 0px;",
                            bslib::card_body(height = "200px", genTPGraphs())
                        )
                    ))
            )
        )
        
        parentInfo <- htmltools::HTML(sprintf("const TPComponParentInfo = JSON.parse('%s');", getTPParentInfoJSON()))
        
        htmltools::tagList(
            htmltools::tags$script(parentInfo),
            makeSideBar("TPsParents TPsByGroup TPsBySuspect", "Parent", "TPCompon-select", "updateTPCompon",
                        getTPComponIDs(), getTPComponNames(), "advanceTPCompon", bd)
        )
    },
    
    genDetailsTPsUI = function()
    {
        hasForms <- !is.null(objects$components[[1]]$candidates[[1]][["formula"]])
        groupBy <- if (hasForms)
        {
            list(
                list(value = "", name = "None"),
                list(value = "formula", name = "Formula"),
                list(value = "formulaDiff", name = "\U0394 Formula")
            )
        }
        else
            NULL
        
        ret <- list(
            genTPsSidebar(),
            makeFGTableCard(genMainTableTPsParents(), TRUE, "TPsParents"),
            makeFGTableCard(genMainTableTPsByGroup(), TRUE, "TPsByGroup")
        )
        
        if (objects$components@fromTPs)
        {
            ret <- append(ret, list(
                makeCandTableCard(genMainTableTPsCandSuspect(), FALSE, "TPsByGroup", "TP candidates",
                                  groupBy = groupBy),
                makeCandTableCard(genMainTableTPsBySuspect(), TRUE, "TPsBySuspect", "TP candidates",
                                  groupBy = groupBy),
                makeFGTableCard(genMainTableTPsCandGroup(), FALSE, "TPsBySuspect")
            ))
        }
        
        return(ret)
    }
)

