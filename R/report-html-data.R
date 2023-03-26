#' @include main.R
#' @include report-html.R
NULL

roundFGTab <- function(ftab, fGroups)
{
    ftab <- copy(ftab)
    for (col in names(ftab)[(sapply(ftab, is.numeric))])
        set(ftab, j = col, value = round(ftab[[col]],
                                         if (col %in% c("mz", "neutralMass", "susp_d_mz")) 5 else 2))
    ftab[, (replicateGroups(fGroups)) := lapply(.SD, round, 0), .SDcols = replicateGroups(fGroups)]
    return(ftab)    
}

getFGTable <- function(fGroups, colSusp)
{
    tab <- if (isScreening(fGroups))
        as.data.table(fGroups, qualities = "score", average = TRUE, collapseSuspects = colSusp)
    else
        as.data.table(fGroups, qualities = "score", average = TRUE)

    tab <- subsetDTColumnsIfPresent(tab, c("group", "ret", "mz", replicateGroups(fGroups), "adduct", "neutralMass",
                                           paste0("susp_", c("name", "estIDLevel", "d_rt", "d_mz", "sets", "InChIKey")),
                                           featureQualityNames(scores = TRUE),
                                           grep("^ISTD_assigned", names(tab), value = TRUE)))
    
    # if (rmdVars$retMin) UNDONE
    tab[, ret := ret / 60]
    if (!is.null(tab[["susp_d_rt"]]))
        tab[, susp_d_rt := susp_d_rt / 60]

    tab <- roundFGTab(tab, fGroups)
    
    if (nrow(internalStandards(fGroups)) > 0)
    {
        wrapISTDs <- function(s) wrapStr(gsub(",", ", ", s, fixed = TRUE), 50)
        if (isFGSet(fGroups))
        {
            for (s in sets(fGroups))
            {
                cn <- paste0("ISTD_assigned-", s)
                tab[, (cn) := sapply(get(cn), wrapISTDs)]
            }
        }
        else if (!is.null(tab[["ISTD_assigned"]]))
            tab[, ISTD_assigned := sapply(ISTD_assigned, wrapISTDs)]
    }
    
    return(tab)
}

getFGColSepStyle <- function() list(borderLeft = "1px solid DarkGrey")

getFGColGrpStartCols <- function(groupDefs) sapply(groupDefs[-1], function(col) col$columns[1])

featGroupTabHasSusps <- function(tab) !is.null(tab[["susp_d_mz"]]) # HACK: this column should always be there if there are (non-collapsed) suspect results

getFeatGroupColDefs <- function(tab)
{
    colDefs <- list()
    
    setCD <- function(col, field, value)
    {
        if (col %chin% names(tab))
        {
            if (is.null(colDefs[[col]]))
                colDefs[[col]] <<- do.call(reactable::colDef, setNames(list(value), field))
            else
                colDefs[[col]][[field]] <<- value
        }
    }

    setCD("mz", "name", "m/z")
    if (featGroupTabHasSusps(tab))
        setCD("susp_name", "name", "suspect")
    else
        setCD("susp_name", "name", "name(s)")
    setCD("susp_d_rt", "name", "\U0394 ret")
    setCD("susp_d_mz", "name", "\U0394 mz")
    setCD("susp_estIDLevel", "name", "estIDLevel") # UNDONE: tool-tip?
    setCD("susp_estIDLevel", "align", "right")
    setCD("susp_sets", "name", "sets")
    # InChIKeys are only there for internal usage
    setCD("susp_InChIKey", "show", FALSE)
    
    featScoreNames <- intersect(featureQualityNames(scores = TRUE), names(tab))
    if (length(featScoreNames) > 0)
    {
        for (col in featScoreNames)
        {
            setCD(col, "name", sub("Score$", "", col))
            setCD(col, "show", FALSE) # hidden by default, controlled by checkbox
        }
    }
    
    return(colDefs)
}

getFGGroupDefs <- function(tab, groupBy, rgs)
{
    colSepStyle <- getFGColSepStyle()
    hasSusp <- featGroupTabHasSusps(tab)
    isGrouped <- !is.null(groupBy)
    featScoreNames <- intersect(featureQualityNames(scores = TRUE), names(tab))
    
    return(pruneList(list(
        # NOTE: the first group doesn't have a headerStyle (left border)
        # workaround for stickies: https://github.com/glin/reactable/issues/236#issuecomment-1107911895
        if (isGrouped) reactable::colGroup("", columns = groupBy, sticky = "left") else NULL,
        reactable::colGroup("feature", columns = intersect(c("group", "ret", "mz", "adduct", "neutralMass",
                                                             grep("^ISTD_assigned", names(tab), value = TRUE)),
                                                           names(tab)),
                            headerStyle = if (isGrouped) colSepStyle else NULL),
        if (hasSusp) reactable::colGroup("screening",
                                         columns = intersect(c("susp_d_rt", "susp_d_mz", "susp_estIDLevel", "susp_sets"),
                                                             names(tab)),
                                         headerStyle = colSepStyle)
        # may still be suspects, but collapsed
        else if (!is.null(tab[["susp_name"]])) reactable::colGroup("suspect", "susp_name", headerStyle = colSepStyle) else NULL,
        reactable::colGroup("intensity", columns = rgs, headerStyle = colSepStyle),
        if (length(featScoreNames) > 0) reactable::colGroup("feature quality scores", columns = featScoreNames,
                                                            headerStyle = colSepStyle) else NULL
    )))
}

reactExactFilter <- function()
{
    htmlwidgets::JS("function(rows, columnId, filterValue)
{
    return rows.filter(row => row.values[columnId] === filterValue)
}")
}

reactSelectFilter <- function(id, values, name)
{
    # from examples
    htmltools::tags$select(
        onchange = sprintf("Reactable.setFilter('%s', '%s', event.target.value || undefined)", id, name),
        tags$option(value = "", "All"),
        lapply(unique(values), tags$option),
        "aria-label" = paste("Filter", name),
        style = "width: 100%; height: 28px;"
    )
}

setReactNumRangeFilters <- function(id, tab, colDefs)
{
    for (col in names(tab))
    {
        if (!is.numeric(tab[[col]]))
            next
        if (is.null(colDefs[[col]]))
            colDefs[[col]] <- reactable::colDef()
        colDefs[[col]]$filterInput <- function(values, name)
        {
            htmltools::tags$button(class = "btn btn-secondary btn-sm", "data-bs-toggle" = "modal",
                                   "data-bs-target" = "#filterRangeModal",
                                   onclick = sprintf("filtRangeModalInit('%s', '%s')", id, name), "Range")
        }
        colDefs[[col]]$filterMethod <- htmlwidgets::JS("function(rows, columnId, filterValue)
        {
            if (filterValue[0] == '')
                filterValue[0] = -Infinity;
            if (filterValue[1] == '')
                filterValue[1] = Infinity;
            return rows.filter(function(row)
            {
                return row.values[columnId] >= filterValue[0] && row.values[columnId] <= filterValue[1];
            })
        }")
    }
    
    return(colDefs)
}

setReactSelRangeFilter <- function(id, colDef)
{
    if (is.null(colDef))
        colDef <- reactable::colDef()
    colDef$filterInput <- function(values, name)
    {
        htmltools::tags$button(class = "btn btn-secondary btn-sm", "data-bs-toggle" = "modal",
                               "data-bs-target" = "#filterSelModal",
                               onclick = sprintf("filtSelModalInit('%s', '%s')", id, name), "Select")
    }
    colDef$filterMethod <- htmlwidgets::JS("function(rows, columnId, filterValue)
    {
        return rows.filter(r => filterValue.has(r.values[columnId]));
    }")
    
    return(colDef)
}

makeReactable <- function(tab, id, bordered = TRUE, pagination = FALSE, ...)
{
    return(reactable::reactable(tab, elementId = id, resizable = TRUE, bordered = bordered, wrap = FALSE,
                                pagination = pagination, showPageSizeOptions = TRUE,
                                pageSizeOptions = c(25, 50, 100, 250, 500), defaultPageSize = 50,...))
}

makeFGReactable <- function(tab, id, colDefs, groupDefs, visible, plots, ...)
{
    # sync column order
    tab <- copy(tab)
    setcolorder(tab, unlist(lapply(groupDefs, "[[", "columns")))
    
    # add EICs
    tab[, c("chrom_small", "chrom_large") := group]
    tabn <- names(tab)
    setcolorder(tab, c(tabn[seq_len(match("group", tabn))], "chrom_small", "chrom_large")) # move after group column
    colDefs$chrom_small <- reactable::colDef("chromatogram", filterable = FALSE, searchable = FALSE,
                                             cell = function(value, index)
    {
        htmltools::img(src = plots$chromsSmall[[value]], style = list("max-height" = "20px"), class = "noZoomImg")
    })
    colDefs$chrom_large <- reactable::colDef("chromatogram", minWidth = 400, show = FALSE, filterable = FALSE,
                                             searchable = FALSE, cell = function(value, index)
    {
        htmltools::img(src = plots$chromsLarge[[value]])
    })
    
    for (i in seq_along(groupDefs))
    {
        if (groupDefs[[i]]$name == "feature")
        {
            groupDefs[[i]]$columns <- c(groupDefs[[i]]$columns, "chrom_small", "chrom_large")
            break
        }
    }
    
    colSepStyle <- getFGColSepStyle()
    grpStartCols <- getFGColGrpStartCols(groupDefs)
    
    bgstyle <- htmlwidgets::JS(sprintf("function(rowInfo, column, state)
{
    let ret = { }
    const gby = state.groupBy;
    if (gby.length !== 0)
    {
        if (rowInfo.level === 0)
        {
            ret.background = 'black';
            ret.color = 'white';
        }
        else if (gby.length > 1 && rowInfo.level === 1 && column.id !== gby[0])
            ret.background = 'LightGrey';
    }
    
    if ([ %s ].includes(column.id))
        ret.borderLeft = '%s';

    return ret;
}", paste0("'", grpStartCols, "'", collapse = ","), colSepStyle))
    
    for (col in grpStartCols)
        colDefs[[col]]$headerStyle <- colSepStyle
    
    colDefs <- lapply(colDefs, function(cd)
    {
        cd$style <- bgstyle
        return(cd)
    })
    
    groupCols <- lapply(groupDefs, "[[", "columns")
    names(groupCols) <- sapply(groupDefs, "[[", "name")
    colToggles <- list(
        group = setdiff(as.vector(groupCols$feature), c("group", "chrom_small", "chrom_large")), # don't toggle group and chroms
        intensities = as.vector(groupCols$intensity),
        qualities = as.vector(groupCols[["feature quality scores"]]),
        chrom_large = "chrom_large"
    )
    
    colDefs <- setReactNumRangeFilters(id, tab, colDefs)

    onClick = htmlwidgets::JS("function(rowInfo, column)
{
    updateFeatTabRowSel(rowInfo.values, rowInfo.index);
}")
    
    headThemeStyle <- list(padding = "2px 4px")
    rt <- makeReactable(tab, id, highlight = TRUE, onClick = onClick, defaultExpanded = TRUE, columns = colDefs,
                        defaultColDef = reactable::colDef(style = bgstyle), columnGroups = groupDefs,
                        filterable = FALSE, pagination = TRUE,
                        theme = reactable::reactableTheme(headerStyle = headThemeStyle,
                                                          groupHeaderStyle = headThemeStyle,
                                                          cellPadding = "2px 4px"),
                        meta = list(selectedRow = 0, plots = plots, colToggles = colToggles),
                        rowStyle = htmlwidgets::JS("function(rowInfo, state)
{
    const sel = state.meta.selectedRow;
    let ret = { cursor: 'pointer' };
    if (sel != null && rowInfo != undefined && rowInfo.index === sel)
    {
        ret.background = '#eee';
        ret.fontWeight = 'bold';
    }
    return ret;
}"), ...)
    
    if (!visible)
        rt <- htmlwidgets::onRender(rt, htmlwidgets::JS("function(el, x) { el.style.display = 'none'; }"))
    
    return(rt)
}

makeReactableCompact <- function(tab, id, ...)
{
    makeReactable(tab, id = id, compact = TRUE, fullWidth = FALSE, striped = TRUE, sortable = FALSE,
                  rowStyle = list(fontSize = "11pt"), ...)
}

makePropTab <- function(tab, sets, idcol = FALSE)
{
    setCols <- if (!is.null(sets))
        grep(paste0("\\-(", paste0(sets, collapse = "|"), ")$"), names(tab), value = TRUE)
    else
        NULL
    haveSets <- !is.null(sets) && length(setCols) > 0
    propTabList <- lapply(split(tab, seq_len(nrow(tab))), function(trow)
    {
        if (!isFALSE(idcol))
            trow <- removeDTColumnsIfPresent(trow, idcol)
        
        if (!haveSets)
            return(setnames(transpose(trow, keep.names = "property"), 2, "value"))
        
        commonRow <- trow[, setdiff(names(trow), setCols), with = FALSE]
        setRows <- sapply(sets, function(s)
        {
            pat <- paste0("\\-", s, "$")
            sr <- trow[, grep(pat, setCols, value = TRUE), with = FALSE]
            if (nrow(sr) == 0) # set lacks data, make non-empty dummy table so it still gets added
                sr <- setnames(data.table(NA), names(trow)[1])
            else
                setnames(sr, sub(pat, "", names(sr)))
            return(sr)
        }, simplify = FALSE)
        ret <- rbindlist(c(setNames(list(commonRow), "common"), setRows), fill = TRUE, idcol = "row")
        return(transpose(ret, keep.names = "property", make.names = "row"))
    })
    
    if (!isFALSE(idcol))
        names(propTabList) <- tab[[idcol]]
    
    return(rbindlist(propTabList, idcol = idcol))
}

makePropReactable <- function(tab, id, idcol = FALSE, minPropWidth = 150, minValWidth = 75, ...)
{
    colDefs <- list()
    if (!isFALSE(idcol))
    {
        # exact match filter
        colDefs[[idcol]] <- reactable::colDef(show = FALSE, filterMethod = reactExactFilter())
    }
    
    for (col in names(tab))
    {
        if (col == idcol)
            next
        colDefs[[col]] <- reactable::colDef(minWidth = if (col == "property") minPropWidth else minValWidth,
                                            align = if (col == "property") "left" else "right", html = TRUE)
        if (col %chin% c("property", "value", "common"))
            colDefs[[col]]$name <- ""
    }
    
    return(makeReactableCompact(tab, id = id, columns = colDefs, ...))
}

titleTab <- function(title, tab)
{
    return(htmltools::div(
        style = list(margin = "10px 20px"),
        htmltools::div(style = list("text-align" = "center", "font-weight" = "bold"), title),
        tab
    ))
}

makeAnnKable <- function(tab, mark = NULL, ...)
{
    withr::local_options(list(knitr.kable.NA = ""))
    kab <- knitr::kable(tab, format = "html", escape = FALSE, ...) %>%
        kableExtra::kable_styling(font_size = 13, bootstrap_options = c("striped", "condensed"))
    # BUG: row_spec needs to be before scroll_box
    if (!is.null(mark))
        kab <- kableExtra::row_spec(kab, mark, italic = TRUE)
    kab <- kableExtra::scroll_box(kab, box_css = "border: 1px solid #ddd;",
                                  extra_css = "overflow-x: hidden; overflow-y: auto; height: 200px;")
    return(htmltools::span(dangerouslySetInnerHTML = list("__html" = kab)))
}

makeAnnDetailsReact <- function(title, tab, sets)
{
    if (nrow(tab) > 0)
    {
        tab <- copy(tab)
        cols <- copy(names(tab))
        for (col in cols)
        {
            if (is.na(tab[[col]]) | (is.character(tab[[col]]) & !nzchar(tab[[col]])))
                set(tab, j = col, value = NULL)
        }
    }
    ptab <- makePropTab(tab, sets, FALSE)
    return(titleTab(title, makeAnnKable(ptab, col.names = sub("common", "", names(ptab), fixed = TRUE))))
}

makeAnnPLReact <- function(apl)
{
    rtab <- if (is.null(apl) || nrow(apl) == 0)
        htmltools::div(align = "center", "No annotation available.")
    else
    {
        isPrec <- apl$precursor
        
        apl[, c("ID", "annotated", "precursor") := NULL]
        apl[, c("mz", "intensity") := .(round(mz, 5), round(intensity))]
        apl[, ion_formula := subscriptFormulaHTML(ion_formula)]
        apl[, neutral_loss := subscriptFormulaHTML(neutral_loss)]
        if (!is.null(apl[["ion_formula_MF"]]))
            apl[, ion_formula_MF := subscriptFormulaHTML(ion_formula_MF)]
        
        colDefs <- pruneList(list(
            ion_formula = reactable::colDef(html = TRUE, minWidth = 125),
            neutral_loss = reactable::colDef(html = TRUE, minWidth = 125),
            ion_formula_MF = if (!is.null(apl[["ion_formula_MF"]])) reactable::colDef(html = TRUE, minWidth = 150) else NULL
        ))
        
        makeAnnKable(apl, mark = if (any(isPrec)) which(isPrec) else NULL)
    }
    
    return(titleTab("Peak list annotations", rtab))
}

makeAnnScoreReact <- function(annRow, sets)
{
    annRow[, (names(annRow)) := lapply(.SD, round, 2), .SDcols = names(annRow)]
    ptab <- makePropTab(annRow, sets, FALSE)
    return(titleTab("Scorings", makeAnnKable(ptab)))
}

getAnnReactImgCell <- function(value) htmltools::img(src = value, style = list("max-height" = "300px"))

makeAnnReactable <- function(tab, id, detailsTabFunc = NULL, annPLTabFunc = NULL, scoreTabFunc = NULL, ...)
{
    details <- if (is.null(detailsTabFunc))
        NULL
    else
    {
        function(index)
        {
            htmltools::div(style = list(margin = "12px 45px", display = "flex", "flex-wrap" = "no-wrap",
                                        background = "#FCFCFC", border = "dashed 1px",
                                        "justify-content" = "space-between", "overflow-x" = "auto"),
                           detailsTabFunc(index), annPLTabFunc(index), scoreTabFunc(index))
        }
    }
    
    return(makeReactable(tab, id, compact = TRUE, details = details, pagination = TRUE,
                         language = reactable::reactableLang(noData = "No annotations available"), ...))
}

reportHTMLUtils$methods(
    genFGTablePlain = function()
    {
        tab <- getFGTable(objects$fGroups, ",")
        groupDefs <- getFGGroupDefs(tab, NULL, replicateGroups(objects$fGroups))
        colDefs <- getFeatGroupColDefs(tab)
        makeFGReactable(tab, "detailsTabPlain", colDefs = colDefs, groupDefs = groupDefs, visible = TRUE, plots = plots)
    },
    genFGTableSuspects = function()
    {
        tab <- getFGTable(objects$fGroups, NULL)[!is.na(susp_name)]
        groupDefs <- getFGGroupDefs(tab, "susp_name", replicateGroups(objects$fGroups))
        colDefs <- getFeatGroupColDefs(tab)
        makeFGReactable(tab, "detailsTabSuspects", colDefs = colDefs, groupDefs = groupDefs, visible = FALSE,
                        plots = plots, groupBy = "susp_name")
    },
    genFGTableISTDs = function()
    {
        istds <- data.table::copy(internalStandards(objects$fGroups))
        istds <- subsetDTColumnsIfPresent(istds, c("name", "group", "InChIKey", "d_rt", "d_mz", "sets"))
        # if (rmdVars$retMin) UNDONE
        istds[, d_rt := d_rt / 60]
        
        # HACK: the ISTD table is essentially the same as what screenInfo() returns for suspects. Rename the columns
        # here, so that groupDefs etc are set like suspects.
        rncols <- setdiff(names(istds), "group")
        setnames(istds, rncols, paste0("susp_", rncols))
        
        ftab <- getFGTable(objects$fGroups, ",")
        ftab <- removeDTColumnsIfPresent(ftab, c("susp_name", "susp_sets", grep("^ISTD_assigned",
                                                                                names(ftab), value = TRUE)))
        tab <- merge(ftab, istds, by = "group")
        tab <- roundFGTab(tab, objects$fGroups)
        
        groupDefs <- getFGGroupDefs(tab, "susp_name", replicateGroups(objects$fGroups))
        colDefs <- getFeatGroupColDefs(tab)
        colDefs$susp_name$name <- "Internal standard" # HACK
        makeFGReactable(tab, "detailsTabISTDs", colDefs = colDefs, groupDefs = groupDefs, visible = FALSE,
                        plots = plots, groupBy = "susp_name")
    },
    genFGTableComponents = function()
    {
        tab <- getFGTable(objects$fGroups, ",")
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
                        plots = plots, groupBy = "component")
    },
    genFGTableTPs = function()
    {
        fromTPs <- objects$components@fromTPs
        
        tabTPsFeat <- getFGTable(objects$fGroups, if (fromTPs) NULL else ",")
        
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
        # UNDONE: convert RT to minutes if needed
        tabTPs[, c("parent_ret", "parent_mz") := groupInfo(objects$fGroups)[parent_group, ]]
        
        for (col in intersect(c("parent_ret", "retDiff", "specSimilarity"), names(tabTPs)))
            tabTPs[, (col) := round(get(col), 2)]
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
        colDefs$formulaDiff <- reactable::colDef(name = "\U0394 formula")
        
        # InChIKeys are only there for internal usage
        if (!is.null(tabTPs[["parent_susp_InChIKey"]]))
            colDefs$parent_susp_InChIKey <- reactable::colDef(show = FALSE)
        # same for cmpIndex
        colDefs$cmpIndex <- reactable::colDef(show = FALSE)

        makeFGReactable(tabTPs, "detailsTabTPs", FALSE, plots, groupBy = groupBy, colDefs = colDefs,
                        groupDefs = groupDefs)
    },

    genSuspInfoTable = function(id)
    {
        tab <- as.data.table(objects$fGroups, collapseSuspects = NULL)
        tab <- subsetDTColumnsIfPresent(tab, paste0("susp_", suspMetaDataCols()))
        
        setnames(tab, sub("^susp_", "", names(tab)))
        
        tab <- unique(tab, by = "name")
        
        for (col in intersect(c("neutralMass", "mz"), names(tab)))
            set(tab, j = col, value = round(tab[[col]], 5))
        if (!is.null(tab[["rt"]]))
            set(tab, j = "rt", value = round(tab$rt, 2))
        if (!is.null(tab[["formula"]]))
            set(tab, j = "formula", value = subscriptFormulaHTML(tab$formula))
        
        ptab <- makePropTab(tab, NULL, "name")
        makePropReactable(ptab, id, "name", minPropWidth = 120, minValWidth = 150)
    },

    genISTDInfoTable = function()
    {
        tab <- data.table::copy(internalStandards(objects$fGroups))
        tab <- subsetDTColumnsIfPresent(tab, suspMetaDataCols())
        tab <- unique(tab, by = "name")
        
        for (col in intersect(c("neutralMass", "mz"), names(tab)))
            set(tab, j = col, value = round(tab[[col]], 5))
        if (!is.null(tab[["rt"]]))
            set(tab, j = "rt", value = round(tab$rt, 2))
        if (!is.null(tab[["formula"]]))
            set(tab, j = "formula", value = subscriptFormulaHTML(tab$formula))
        
        ptab <- makePropTab(tab, NULL, "name")
        makePropReactable(ptab, "ISTDInfoTab", "name", minPropWidth = 120, minValWidth = 150)
    },
    
    genComponentInfoTable = function()
    {
        tab <- componentInfo(objects$components)
        tab <- removeDTColumnsIfPresent(tab, c("links", "size"))
        
        for (col in intersect(c("neutral_mass", "mz_increment"), names(tab)))
            set(tab, j = col, value = round(tab[[col]], 5))
        for (col in intersect(c("cmp_ret", "cmp_retsd", "cmp_ppm", "ret_increment", "ret_min", "ret_max", "ret_range"),
                              names(tab)))
            set(tab, j = col, value = round(tab[[col]], 2))
        
        ptab <- makePropTab(tab, NULL, "name")
        makePropReactable(ptab, "componentInfoTab", "name", minPropWidth = 120, minValWidth = 150)
    },
        
    genFeaturesTable = function()
    {
        tab <- as.data.table(getFeatures(objects$fGroups))
        tab <- removeDTColumnsIfPresent(tab, "adduct") # can already be seen in group table
        tab <- removeDTColumnsIfPresent(tab, featureQualityNames(group = FALSE, scores = FALSE)) # only show scores
      
        for (col in names(tab)[sapply(tab, is.numeric)])
            set(tab, j = col, value = round(tab[[col]], if (col %in% c("mz", "mzmin", "mzmax")) 5 else 2))

        anaInfo <- analysisInfo(objects$fGroups)
        
        # add data for 'missing' analyses
        missingTab <- rbindlist(sapply(unique(tab$group), function(grp)
        {
            mt <- data.table(analysis = setdiff(analyses(objects$fGroups), tab[group == grp]$analysis))
            if (!is.null(tab[["set"]]))
                mt[, set := anaInfo[match(analysis, anaInfo$analysis), "set"]]
            return(mt)
        }, simplify = FALSE), idcol = "group")
        if (nrow(missingTab) > 0)
        {
            tab <- rbind(tab, missingTab, fill = TRUE)
            # make sure analyses retain order of anaInfo
            tab[, anaInd := match(analysis, anaInfo$analysis)]
            setorderv(tab, "anaInd")
            tab[, anaInd := NULL]
        }
        
        tab[, rGroup := anaInfo[match(analysis, anaInfo$analysis), "group"]]
        
        # add EICs
        tab[, chromatogram := ""] # dummy value, not needed
        tabn <- names(tab)
        
        setcolorder(tab, c("analysis", "rGroup", "ID", "chromatogram"))

        colDefs <- list(
            group = reactable::colDef(show = FALSE, filterMethod = reactExactFilter()),
            rGroup = reactable::colDef("replicate group"),
            chromatogram = reactable::colDef(minWidth = 175, cell = function(value, index)
            {
                htmltools::img(src = plots$chromsFeatures[[tab$group[index]]][[tab$analysis[index]]])
            })
        )
        if (!is.null(tab[["set"]]))
        {
            colDefs$set <- reactable::colDef(filterInput = function(values, name) reactSelectFilter("featuresTab",
                                                                                                    values, name))
        }
        
        fqn <- featureQualityNames(group = FALSE, scores = TRUE)
        for (col in fqn)
        {
            if (!is.null(tab[[col]]))
                colDefs[[col]] <- reactable::colDef(show = FALSE)
        }
        
        colDefs <- setReactNumRangeFilters("featuresTab", tab, colDefs)
        
        colDefs$analysis <- setReactSelRangeFilter("featuresTab", reactable::colDef())
        colDefs$rGroup <- setReactSelRangeFilter("featuresTab", colDefs$rGroup)
        
        makeReactable(tab, "featuresTab", compact = TRUE, defaultExpanded = TRUE, columns = colDefs, filterable = FALSE,
                      meta = list(featQualCols = fqn), pagination = TRUE)
    },
    
    genMSPLTable = function(MSLevel)
    {
        MSPeakLists <- objects$MSPeakLists[, names(objects$fGroups)]
        
        id <- if (MSLevel == 2) "MSMSPLTab" else "MSPLTab"
        theight <- 350
        
        if (length(MSPeakLists) == 0)
        {
            # dummy table to show empty results
            return(makeAnnReactable(data.table(ID = integer()), id, height = theight))
        }
        
        tab <- as.data.table(MSPeakLists)
        tab <- tab[type == if (MSLevel == 2) "MSMS" else "MS"][, type := NULL]
        isPrec <- tab$precursor
        tab[, precursor := NULL]
        
        colDefs <- list(
            group = reactable::colDef(show = FALSE, filterMethod = reactExactFilter()),
            mz = reactable::colDef(format = reactable::colFormat(digits = 5)),
            intensity = reactable::colDef(format = reactable::colFormat(digits = 0))
        )
        
        colDefs <- setReactNumRangeFilters(id, tab, colDefs)
        if (!is.null(tab[["set"]]))
            colDefs$set <- reactable::colDef(filterInput = function(values, name) reactSelectFilter(id, values, name))
        
        return(makeAnnReactable(tab, id, columns = colDefs, height = theight, filterable = TRUE,
                                rowClass = function(index) if (isPrec[index]) "fw-bold" else ""))
    },
    
    genFormulasTable = function()
    {
        formulas <- objects$formulas[names(objects$fGroups)]
        
        mcn <- mergedConsensusNames(formulas)
        
        if (length(formulas) == 0)
        {
            # dummy table to show empty results
            return(makeAnnReactable(data.table(formula = character()), "formulasTab"))
        }
        
        tab <- as.data.table(formulas)
        
        # NOTE: for consensus results, duplicate algo columns (eg explainedPeaks) are only shown in details
        tab <- subsetDTColumnsIfPresent(tab, c("group", "neutral_formula", "neutralMass", "explainedPeaks",
                                               "explainedIntensity", "error"))
        
        tab[, candidate := seq_len(.N), by = "group"]
        
        for (col in names(tab))
        {
            if (is.numeric(tab[[col]]))
                tab[, (col) := round(get(col), if (col == "neutralMass") 5 else 2)]
        }
        
        tab[, spectrum := plots$formulas[[group]]$spectra, by = "group"]
        tab[, scorings := plots$formulas[[group]]$scores, by = "group"]

        if (hasSuspects())
        {
            scr <- screenInfo(objects$fGroups)
            if (!is.null(scr[["formula"]]))
            {
                tab[neutral_formula %chin% scr$formula, suspect := {
                    wrapStr(paste0(scr[match(neutral_formula, formula)]$name, collapse = ", "), 50)
                }, by = "neutral_formula"]
            }
        }

        tab[, neutral_formula := subscriptFormulaHTML(neutral_formula)]
        
        getFormDetails <- function(index)
        {
            ft <- formulas[[tab$group[index]]][tab$candidate[index]]
            ft <- ft[, setdiff(names(ft), names(tab)), with = FALSE]
            takeCols <- getAllMergedConsCols(c("neutral_formula", "ion_formula", "neutralMass", "ion_formula_mz",
                                               "error", "error_frag_median", "error_frag_median_abs",
                                               "explainedPeaks", "explainedIntensity"), names(ft), mcn)
            ft <- ft[, takeCols, with = FALSE]
            for (col in getAllMergedConsCols(c("neutral_formula", "ion_formula"), names(ft), mcn))
                set(ft, j = col, value = subscriptFormulaHTML(ft[[col]]))
            massCols <- getAllMergedConsCols(c("neutralMass", "ion_formula_mz"), names(ft), mcn)
            numCols <- names(ft)[sapply(ft, is.numeric)]
            for (col in numCols)
                set(ft, j = col, value = round(ft[[col]], if (col %chin% massCols) 5 else 2))
            return(makeAnnDetailsReact("Formula properties", ft, sets(objects$fGroups)))
        }
        
        getAnnPLDetails <- function(index)
        {
            apl <- annotatedPeakList(formulas, index = tab$candidate[index], groupName = tab$group[index],
                                     MSPeakLists = objects$MSPeakLists, onlyAnnotated = TRUE)
            return(makeAnnPLReact(apl))
        }
        
        getScoreDetails <- function(index)
        {
            fRow <- formulas[[tab$group[index]]][tab$candidate[index]]
            cols <- getAllMergedConsCols(annScoreNames(formulas, FALSE), names(fRow), mcn)
            return(makeAnnScoreReact(fRow[, cols, with = FALSE], sets(objects$fGroups)))
        }
        
        setcolorder(tab, intersect(c("candidate", "suspect"), names(tab)))
        
        colDefs <- pruneList(list(
            group = reactable::colDef(show = FALSE, filterMethod = reactExactFilter()),
            candidate = reactable::colDef("#", minWidth = 15),
            suspect = if (!is.null(tab[["suspect"]])) reactable::colDef("suspect(s)") else NULL,
            neutral_formula = reactable::colDef("formula", html = TRUE),
            neutralMass = reactable::colDef("neutral mass"),
            spectrum = reactable::colDef(cell = getAnnReactImgCell, minWidth = 200),
            scorings = reactable::colDef(cell = getAnnReactImgCell, minWidth = 200)
        ))
        
        colDefs <- setReactNumRangeFilters("formulasTab", tab, colDefs)
        
        return(makeAnnReactable(tab, "formulasTab", columns = colDefs, getFormDetails, getAnnPLDetails,
                                getScoreDetails))
    },
    
    genCompoundsTable = function()
    {
        compounds <- objects$compounds[names(objects$fGroups)]
        
        mcn <- mergedConsensusNames(compounds)
        
        if (length(compounds) == 0)
        {
            # dummy table to show empty results
            return(makeAnnReactable(data.table(compound = character()), "compoundsTab"))
        }

        tab <- as.data.table(compounds)
        
        # NOTE: for consensus results, duplicate algo columns (eg identifier) are only shown in details
        tab <- subsetDTColumnsIfPresent(tab, c("group", "compoundName", "compoundName2", "identifier", "database",
                                               "neutral_formula", "neutralMass", "explainedPeaks", "score", "InChIKey"))
        
        tab[, candidate := seq_len(.N), by = "group"]
        cmpNames2 <- tab[["compoundName2"]]
        if (!is.null(cmpNames2))
            tab[, compoundName2 := NULL]
        
        if (!is.null(tab[["identifier"]]))
            tab[, identifier := sapply(identifier, makeDBIdentLink, db = database[1])][, database := NULL]
        
        tab[, neutral_formula := subscriptFormulaHTML(neutral_formula)]
        tab[, neutralMass := round(neutralMass, 5)]
        if (!is.null(tab[["score"]]))
            tab[, score := round(score, 2)]
        
        tab[, structure := plots$structs[InChIKey]]
        tab[, spectrum := plots$compounds[[group]]$spectra, by = "group"]
        tab[, scorings := plots$compounds[[group]]$scores, by = "group"]
        
        if (hasSuspects())
        {
            scr <- screenInfo(objects$fGroups)
            if (!is.null(scr[["InChIKey"]]))
            {
                tab[InChIKey %chin% scr$InChIKey, suspect := {
                    IK <- InChIKey[1]
                    wrapStr(paste0(scr[match(IK, InChIKey)]$name, collapse = ", "), 50)
                }, by = "InChIKey"]
            }
        }
        
        tab[, InChIKey := NULL]
        
        getCompCell <- function(value, index)
        {
            cn2 <- cmpNames2[index]
            if (!is.null(cn2) && !is.na(cn2))
            {
                if (is.na(value))
                    value <- cn2
                value <- htmltools::tagList(htmltools::strong(value), htmltools::br(), paste0("(", cn2, ")"))
            }
            else
                value <- htmltools::strong(value)
            return(value)
        }
        
        getCompDetails <- function(index)
        {
            ct <- data.table::copy(compounds[[tab$group[index]]][tab$candidate[index]])
            ct <- ct[, setdiff(names(ct), names(tab)), with = FALSE]
            takeCols <- c(
                "compoundName",
                "identifier",
                "relatedCIDs",
                "database",
                "explainedPeaks",
                "neutral_formula",
                "SMILES",
                "InChI",
                "InChIKey",
                "XlogP", "AlogP", "LogP",
                
                # PubChemLite
                "FP", "compoundName2", # UNDONE: update?
                
                # CompTox
                "CASRN", "QCLevel",
                
                # FOR-IDENT
                "tonnage", "categories",
                
                # TP DB
                "parent", "transformation", "enzyme", "evidencedoi"
            )
            ct <- ct[, getAllMergedConsCols(takeCols, names(ct), mcn), with = FALSE]
            
            if (!is.null(tab[["compoundName"]]) && !is.null(tab[["compoundName2"]]))
                set(ct, j = "compoundName2", value = NULL) # already in main table compoundName column
            
            set(ct, j = "neutral_formula", value = subscriptFormulaHTML(ct$neutral_formula))
            
            numCols <- names(ct)[sapply(ct, is.numeric)]
            for (col in numCols)
                set(ct, j = col, value = round(ct[[col]], 2))
            
            for (col in getAllMergedConsCols(c("identifier", "relatedCIDs"), names(ct), mcn))
            {
                # get db: handle consensus results
                cn <- sub("identifier|relatedCIDs", "", col)
                db <- if (nzchar(cn)) ct[[paste0("database", cn)]] else ct$database
                set(ct, j = col, value = makeDBIdentLink(db, ct[[col]]))
            }
            
            return(makeAnnDetailsReact("Compound properties", ct, sets(objects$fGroups)))
        }
        
        getAnnPLDetails <- function(index)
        {
            apl <- annotatedPeakList(compounds, index = tab$candidate[index], groupName = tab$group[index],
                                     MSPeakLists = objects$MSPeakLists, formulas = objects$formulas,
                                     onlyAnnotated = TRUE)
            return(makeAnnPLReact(apl))
        }
        
        getScoreDetails <- function(index)
        {
            cRow <- compounds[[tab$group[index]]][tab$candidate[index]]
            cols <- getAllMergedConsCols(annScoreNames(compounds, FALSE), names(cRow), mcn)
            return(makeAnnScoreReact(cRow[, cols, with = FALSE], sets(objects$fGroups)))
        }
        
        setcolorder(tab, intersect(c("candidate", "compoundName", "suspect", "structure"), names(tab)))
        
        colDefs <- pruneList(list(
            group = reactable::colDef(show = FALSE, filterMethod = reactExactFilter()),
            candidate = reactable::colDef("#", minWidth = 15),
            compoundName = if (!is.null(tab[["compoundName"]])) reactable::colDef("compound", cell = getCompCell) else NULL,
            suspect = if (!is.null(tab[["suspect"]])) reactable::colDef("suspect(s)") else NULL,
            identifier = if (!is.null(tab[["identifier"]])) reactable::colDef(html = TRUE) else NULL,
            neutral_formula = reactable::colDef("formula", html = TRUE),
            neutralMass = reactable::colDef("neutral mass"),
            structure = reactable::colDef(cell = getAnnReactImgCell, minWidth = 125),
            spectrum = reactable::colDef(cell = getAnnReactImgCell, minWidth = 200),
            scorings = reactable::colDef(cell = getAnnReactImgCell, minWidth = 200)
        ))
        
        colDefs <- setReactNumRangeFilters("compoundsTab", tab, colDefs)

        mfWebLinks <- sapply(names(objects$fGroups), function(grp)
        {
            if (is(compounds, "compoundsMF") && grp %chin% names(compounds))
            {
                # make link with used MF settings
                set <- settings(compounds)
            }
            else
                set <- NULL
            
            return(buildMFLandingURL(set, objects$MSPeakLists[[grp]][["MSMS"]],
                                     groupInfo(objects$fGroups)[grp, "mzs"]))
        })
                
        return(makeAnnReactable(tab, "compoundsTab", columns = colDefs, getCompDetails, getAnnPLDetails,
                                getScoreDetails, meta = list(mfWebLinks = mfWebLinks)))
    },
    
    genSuspAnnTable = function()
    {
        tab <- as.data.table(objects$fGroups, collapseSuspects = NULL)
        mcn <- mergedConsensusNames(objects$fGroups)
        tab <- tab[, unique(c("group", "susp_name", getAllSuspCols(paste0("susp_", suspAnnCols()), names(tab), mcn))),
                   with = FALSE]
        tab[, suspID := paste0(susp_name, "-", group)][, c("susp_name", "group") := NULL]
        setnames(tab, sub("^susp_", "", names(tab)))
        
        rcols <- getAllSuspCols(c("annSimForm", "annSimComp", "annSimBoth", "maxFragMatchesRel"), names(tab), mcn)
        if (length(rcols) > 0)
            tab[, (rcols) := lapply(.SD, round, 2), .SDcols = rcols]
        
        ptab <- makePropTab(tab, if (isFGSet(objects$fGroups)) sets(objects$fGroups) else NULL, "suspID")
        makePropReactable(ptab, "suspAnnTab", "suspID", minPropWidth = 150, minValWidth = 100)
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
    }
)

