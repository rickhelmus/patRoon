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
        if (objects$components@fromTPs)
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
    
    genFGTableTPsOld = function()
    {
        # UNDONE: remove function
        
        
        fromTPs <- getTPComponObj()@fromTPs
        
        tabTPsFeat <- getFGTable(objects$fGroups, if (fromTPs) NULL else ",", settings$features$retMin,
                                 settings$features$aggregateConcs, settings$features$aggregateTox)
        
        tabCompon <- as.data.table(getTPComponObj(), candidates = TRUE)
        tabCompon <- tabCompon[parent_group %chin% names(objects$fGroups)]
        tabCompon <- subsetDTColumnsIfPresent(tabCompon, c("name", "parent_name", "parent_group", "group", "TP_retDir",
                                                           "candidate_name", "retDir", "retDiff", "mzDiff", "formulaDiff",
                                                           "specSimilarity", "mergedBy"))
        tabCompon[, cmpIndex := seq_len(.N), by = "name"]
        
        if (fromTPs)
        {
            tabTPs <- merge(tabCompon, tabTPsFeat, by.x = c("group", "candidate_name"), by.y = c("group", "susp_name"),
                            sort = FALSE)
            setnames(tabTPs, c("name", "candidate_name", "parent_name"),
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
                       list(reactable::colGroup("TP", columns = intersect(c("candidate_name", "retDiff", "mzDiff",
                                                                            "formulaDiff", "retDir",
                                                                            "specSimilarity", "mergedBy"),
                                                                          names(tabTPs)),
                                                headerStyle = getMainReactColSepStyle())),
                       groupDefs[seq(3, length(groupDefs))])
        
        colDefs <- getFeatGroupColDefs(tabTPs, rgs)
        
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
        if (!is.null(tabTPs[["candidate_name"]]))
            colDefs$candidate_name <- reactable::colDef("name")
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
        
        makeFGReactable(tabTPs, "detailsTabTPs", plots, settings = settings, objects = objects, groupBy = groupBy,
                        colDefs = colDefs, groupDefs = groupDefs)
    },

    genFGTableTPsParents = function()
    {
        compObj <- getTPComponObj()
        cInfo <- componentInfo(compObj)
        cInfo <- subsetDTColumnsIfPresent(cInfo, c("name", "parent_name", "parent_group", "parent_InChIKey"))
        tab <- getFGTable(objects$fGroups[, cInfo$parent_group], NULL, settings$features$retMin,
                          settings$features$aggregateConcs, settings$features$aggregateTox)
        
        tab <- if (!is.null(cInfo[["parent_name"]]))
            merge(tab, cInfo, by.x = c("group", "susp_name"), by.y = c("parent_group", "parent_name"), sort = FALSE)
        else
            merge(tab, cInfo, by.x = "group", by.y = "parent_group", sort = FALSE)
        
        setnames(tab, "name", "component")

        rgs <- replicateGroups(objects$fGroups)
        groupDefs <- list(reactable::colGroup("component", columns = "component"))
        groupDefs <- append(groupDefs, getFGGroupDefs(tab, NULL, rgs))
        groupDefs <- append(groupDefs, getScrGroupDefs(tab), after = 2)
        colDefs <- getFeatGroupColDefs(tab, rgs)
        colDefs <- modifyList(colDefs, getScrColDefs(tab))
        
        if (!is.null(tab[["parent_InChIKey"]]))
        {
            tab[, structure := plots$structs[parent_InChIKey]]
            tab[, parent_InChIKey := NULL]
            colDefs$structure <- reactable::colDef(cell = getReactImgCell, minWidth = 125)
            groupDefs[[2]]$columns <- append(groupDefs[[2]]$columns, "structure", after = 1) # set after susp name
        }
        
        makeFGReactable(tab, "detailsTabTPsParents", plots, settings = settings, objects = objects, colDefs = colDefs,
                        groupDefs = groupDefs, updateRowFunc = "updateTabRowSelTPsByParents", initView = "TPsParents")
    },
    
    genFGTableTPs = function()
    {
        tabTPsFeat <- getFGTable(objects$fGroups, ",", settings$features$retMin,
                                 settings$features$aggregateConcs, settings$features$aggregateTox)
        
        tabCompon <- as.data.table(getTPComponObj())
        tabCompon <- subsetDTColumnsIfPresent(tabCompon, c("name", "group", "retDir", "retDiff", "mzDiff",
                                                           "specSimilarity"))
        tabCompon[, cmpIndex := seq_len(.N), by = "name"]
        
        tabTPs <- merge(tabCompon, tabTPsFeat, by = "group", sort = FALSE)
        setnames(tabTPs, "name", "component")
        
        if (settings$features$retMin)
            tabTPs[, retDiff := retDiff / 60]
        tabTPs[, retDiff := round(retDiff, 2)]
        tabTPs[, mzDiff := round(mzDiff, 5)]
        if (!is.null(tabTPs[["specSimilarity"]]))
            tabTPs[, specSimilarity := round(specSimilarity, 2)]
        
        rgs <- replicateGroups(objects$fGroups)
        groupDefs <- getFGGroupDefs(tabTPs, NULL, rgs)
        # squeeze in TP column
        groupDefs <- append(groupDefs,
                            list(reactable::colGroup("TP", columns = intersect(c("retDiff", "mzDiff", "retDir",
                                                                                 "specSimilarity"),
                                                                               names(tabTPs)),
                                                     headerStyle = getMainReactColSepStyle())),
                            after = 1)
        
        colDefs <- getFeatGroupColDefs(tabTPs, rgs)
        
        colDefs$retDiff <- reactable::colDef(name = "\U0394 ret")
        colDefs$mzDiff <- reactable::colDef(name = "\U0394 mz")
        
        # internally used variables
        colDefs$component <- reactable::colDef(show = FALSE, filterMethod = reactExactFilter())
        colDefs$cmpIndex <- reactable::colDef(show = FALSE)
        
        makeFGReactable(tabTPs, "detailsTabTPs", plots, settings = settings, objects = objects, colDefs = colDefs,
                        groupDefs = groupDefs, updateRowFunc = "updateTabRowSelTPsByGroup", initView = "TPsByGroup")
    },
    
    genTPCandidatesTable = function()
    {
        # UNDONE: don't call this when !fromTPs
        
        candidatesTab <- rbindlist(lapply(componentTable(getTPComponObj()), function(cmp)
        {
            allc <- rbindlist(cmp$candidates, fill = TRUE, idcol = "groupInd")
            allc[, group := cmp$group[groupInd]][, groupInd := NULL]
            return(allc)
        }), fill = TRUE, idcol = "component")
        
        candidatesTab <- subsetDTColumnsIfPresent(candidatesTab, c("component", "group", "name", "TP_retDir",
                                                                   "formula", "formulaDiff", "fragmentMatches",
                                                                   "neutralLossMatches", "mergedBy", "InChIKey"))
        
        # UNDONE: need this?
        # tabCompon[, cmpIndex := seq_len(.N), by = "name"]
        
        scr <- data.table::copy(screenInfo(objects$fGroups))
        scr <- subsetDTColumnsIfPresent(scr, c("group", "name", "estIDLevel", "d_rt", "d_mz", "sets"))
        
        candidatesTab <- merge(candidatesTab, scr, by = c("group", "name"), sort = FALSE)

        # UNDONE: add susp_ prefix, mainly for getScrGroupDefs()/roundScrTab() --> keep this requirement?
        scols <- setdiff(names(scr), "group")
        setnames(candidatesTab, scols, paste0("susp_", scols))
        
        candidatesTab <- roundScrTab(candidatesTab)

        # add parent intensities & screening info
        # UNDONE: add this to parent side widget
        # rgs <- replicateGroups(objects$fGroups)
        # tabTPsPar <- unique(subsetDTColumnsIfPresent(tabTPsFeat,
        #                                              c("group", rgs,
        #                                                paste0("susp_", c("estIDLevel", "d_rt", "d_mz", "sets",
        #                                                                  "InChIKey")))),
        #                     by = "group")
        # setnames(tabTPsPar, paste0("parent_", names(tabTPsPar)))
        # tabTPs <- merge(tabTPs, tabTPsPar, by = "parent_group", sort = FALSE, all.x = TRUE)

        groupDefs <- getScrGroupDefs(candidatesTab)
        # squeeze in TP column
        groupDefs <- append(groupDefs,
                            list(reactable::colGroup("TP", columns = intersect(c("TP_retDir", "formula", "formulaDiff",
                                                                                 "fragmentMatches",
                                                                                 "neutralLossMatches", "mergedBy"),
                                                                               names(candidatesTab)),
                                                     headerStyle = getMainReactColSepStyle())))
        
        colDefs <- getScrColDefs(candidatesTab)
        
        if (!is.null(candidatesTab[["InChIKey"]]))
        {
            candidatesTab[, structure := plots$structs[InChIKey]]
            colDefs$structure <- reactable::colDef(cell = getReactImgCell, minWidth = 125)
            groupDefs[[1]]$columns <- append(groupDefs[[1]]$columns, "structure", after = 1) # set after TP susp name
        }
        
        formCell <- function(value) htmltools::span(dangerouslySetInnerHTML = list("__html" = subscriptFormulaHTML(value, charges = FALSE)))
        if (!is.null(candidatesTab[["formula"]]))
            colDefs$formula <- reactable::colDef(cell = formCell)
        if (!is.null(candidatesTab[["formulaDiff"]]))
            colDefs$formulaDiff <- reactable::colDef(name = "\U0394 formula", cell = formCell)
        
        # internally used variables
        for (col in c("component", "group", "InChIKey"))
        {
            if (!is.null(candidatesTab[[col]]))
                colDefs[[col]] <- reactable::colDef(show = FALSE)
        }
        colDefs$component$filterMethod <- colDefs$group$filterMethod <- reactExactFilter()
        
        makeMainResultsReactable(candidatesTab, "TPCandidatesTab", colDefs = colDefs, groupDefs = groupDefs,
                                 updateRowFunc = "updateTabRowSelTPsCandidates", meta = list())
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
        bslib::card(
            id = "parentCard",
            class = "detailsSidebar",
            detailsViewOfParent = "TPsParents TPsByGroup TPsBySuspect",
            bslib::card_header("Parent"),
            bsCardBodyNoFill(
                htmltools::withTags({
                    div(
                        style = "display: flex;",
                        select(name = "TPCompon", id = "TPCompon-select", onChange = "updateTPCompon(this.value)",
                               style = "min-width: 0;",
                               Map(getTPComponIDs(), getTPComponNames(), f = function(v, n) option(value = v, n))),
                        
                        div(class = "btn-group btn-group-sm mx-1", role = "group", "aria-label" = "advance TP compon group",
                            style = "align-self: flex-start;",
                            button(type = "button", class = "btn btn-light btn-sm", onClick = "advanceTPCompon(-1)",
                                   bsicons::bs_icon("caret-left-fill")),
                            button(type = "button", class = "btn btn-light btn-sm", onClick = "advanceTPCompon(1)",
                                   bsicons::bs_icon("caret-right-fill"))
                        ),
                    )
                })
            ),
            
            bslib::card_body(
                padding = 0,
                pruneUI(bslib::accordion,
                         maybeInclUI(settings$features$chromatograms$large, bslib::accordion_panel(
                             "Chromatogram",
                             bslib::card_body_fill(htmltools::img(id = "chrom_view-parent"))
                         )),
                         # UNDONE: only include if there actually are structures
                         maybeInclUI(hasComponentsFromTPs(), bslib::accordion_panel(
                             "Structure",
                             bsCardBodyNoFill(
                                 htmltools::img(id = "struct_view-parent", style = "min-width: 20%;")
                             ),
                         )),
                         maybeInclUI(hasComponentsFromTPs(), bslib::accordion_panel(
                             "Screening",
                             bsCardBodyNoFill(
                                 genSuspInfoTable("parentInfoTab"),
                             ),
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
        )
    },
    
    genDetailsTPsUI = function()
    {
        list(
            genTPsSidebar(),
            
            bslib::card(
                class = "detailsMainTableOnly",
                detailsViewOfParent = "TPsParents",
                full_screen = TRUE,
                bslib::card_header("Parents"),
                bsCardBodyNoFill(makeFGToolbar("detailsTabTPsParents")),
                bslib::card_body(genFGTableTPsParents())
            ),
            
            
            bslib::card(
                class = "detailsMainTable",
                detailsViewOfParent = "TPsByGroup",
                full_screen = TRUE,
                bslib::card_header("TPs Feature groups"),
                bsCardBodyNoFill(makeFGToolbar("detailsTabTPs")),
                bslib::card_body(genFGTableTPs())
            ),
            
            bslib::card(
                class = "detailsCandTable",
                detailsViewOfParent = "TPsByGroup TPsBySuspect",
                full_screen = TRUE,
                bslib::card_header("TP candidates"),
                bsCardBodyNoFill(
                    makeToolbar("TPCandidatesTab", groupBy = list(
                        list(value = "", name = "None"),
                        list(value = "formula", name = "Formula"),
                        list(value = "formulaDiff", name = "\U0394 Formula")
                    ), toggleExpand = TRUE)
                ),
                bslib::card_body_fill(genTPCandidatesTable())
            )
        )
    }
)

