getMainReactColSepStyle <- function() list(borderLeft = "1px solid DarkGrey")

getReactColGrpStartCols <- function(groupDefs) sapply(groupDefs[-1], function(col) col$columns[1])

getReactImgCell <- function(value) htmltools::img(src = value, style = list("max-height" = "300px"))

reactSelFilterButton <- function(id, name, target, ocFunc, title)
{
    htmltools::tags$button(class = "btn btn-secondary btn-sm", "data-bs-toggle" = "modal",
                           "data-bs-target" = target,
                           onclick = sprintf("%s('%s', '%s')", ocFunc, id, name), title)
}

reactExactFilter <- function()
{
    htmlwidgets::JS("function(rows, columnId, filterValue)
{
    return rows.filter(row => row.values[columnId] === filterValue);
}")
}

reactSuspectFilter <- function()
{
    htmlwidgets::JS("function(rows, columnId, filterValue)
{
    return rows.filter(row => row.values[columnId] && row.values[columnId].split(', ').includes(filterValue));
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
        colDefs[[col]]$filterInput <- function(values, name) reactSelFilterButton(id, name, "#filterRangeModal",
                                                                                  "filtRangeModalInit", "Range")
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
    colDef$filterInput <- function(values, name) reactSelFilterButton(id, name, "#filterSelModal",
                                                                      "filtColSelModalInit", "Select")
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
                sr <- data.table(dummy = NA)
            else
                setnames(sr, sub(pat, "", names(sr)))
            return(sr)
        }, simplify = FALSE)
        ret <- rbindlist(c(setNames(list(commonRow), "common"), setRows), fill = TRUE, idcol = "row")
        if (!is.null(ret[["dummy"]]))
            ret[, dummy := NULL]
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

makeMainResultsReactable <- function(tab, id, colDefs, groupDefs, visible, updateRowFunc, meta, ...)
{
    # sync column order
    tab <- copy(tab)
    setcolorder(tab, unlist(lapply(groupDefs, "[[", "columns")))
    
    colSepStyle <- getMainReactColSepStyle()
    grpStartCols <- getReactColGrpStartCols(groupDefs)
    
    bgstyle <- htmlwidgets::JS(sprintf("function(rowInfo, column, state)
{
    let ret = { }
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

    colDefs <- setReactNumRangeFilters(id, tab, colDefs)
    
    onClick = htmlwidgets::JS(sprintf("function(rowInfo, column)
{
    %s(rowInfo.values, rowInfo.index);
}", updateRowFunc))
    
    headThemeStyle <- list(padding = "2px 4px")
    rt <- makeReactable(tab, id, highlight = TRUE, onClick = onClick, columns = colDefs,
                        defaultColDef = reactable::colDef(style = bgstyle), columnGroups = groupDefs,
                        filterable = FALSE, pagination = TRUE,
                        theme = reactable::reactableTheme(headerStyle = headThemeStyle,
                                                          groupHeaderStyle = headThemeStyle,
                                                          cellPadding = "2px 4px"),
                        meta = modifyList(meta, list(selectedRow = 0)),
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


getHTMLReportPlotPath <- function(outPath)
{
    plotPath <- file.path(outPath, "report_files", "plots")
    mkdirp(plotPath)
    return(plotPath)
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

makeHTMLReportPlot <- function(outPrefix, outPath, func, args, parParams = list(), cache = TRUE, doPrint = FALSE, ...)
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
    
    return(ppath)
}

bsCardBodyNoFill <- function(...)
{
    if (packageVersion("bslib") > "0.4.2")
    {
        # this was the default for bslib 0.4.2
        bslib::card_body(fillable = FALSE, fill = FALSE, ...)
    }
    else
        bslib::card_body(...)
}

reportHTMLUtils$methods(
    makeToolbar = function(tableID, groupBy = NULL, columnToggles = NULL, toggleExpand = FALSE,
                           toggleExpandDisableIfNoGrouping = TRUE, tagsPreButtons = NULL, tagsPostButtons = NULL)
    {
        makeID <- function(id) paste0(tableID, "-", id)
        
        htmltools::withTags({
            ret <- div(class = "pt-2 pb-1 ps-1")
            
            toggleExpandID <- if (toggleExpand) makeID("toggleExpand") else NULL
            
            if (!is.null(groupBy))
            {
                id <- makeID("groupBy")
                onChange <- sprintf("setTabGroupBy('%s', this.value, %s)", tableID,
                                    if (is.null(toggleExpandID)) "undefined" else paste0('"', toggleExpandID, '"'))
                ret <- htmltools::tagAppendChildren(
                    ret,
                    label("for" = makeID("groupBy"), "Group by"),
                    select(id = makeID("groupBy"), onChange = onChange,
                           style = list("margin-right" = "10px", "margin-top" = "3px"),
                           lapply(pruneList(groupBy), function(o) option(value = o$value, o$name))
                    )
                )
            }
            
            if (!is.null(columnToggles))
            {
                ret <- htmltools::tagAppendChildren(
                    ret,
                    lapply(pruneList(columnToggles), function(tc)
                    {
                        id <- makeID(paste0("cols-", tc$value))
                        htmltools::tagList(
                            input(type = "checkbox", id = id,
                                  checked = if (!is.null(tc[["checked"]])) tc$checked else NULL,
                                  onChange = sprintf('showTabCols("%s", "%s", this.checked)', tableID, tc$value)),
                            label("for" = id, tc$name)
                        )
                    })
                )
            }
            
            id <- makeID("filter")
            ret <- htmltools::tagAppendChildren(
                ret,
                
                # UNDONE: make JS function that accepts table ID --> requires metadata that says which columns should be left alone
                input(type = "checkbox", id = id, onChange = 'toggleFGFilters(this.checked)'),
                label("for" = id, style = list("margin-right" = "10px"), "Filters"),
                
                tagsPreButtons,
                
                htmltools::tags$button(type = "button", class = "btn btn-secondary btn-sm",
                                       onClick = sprintf("downloadCSV('%s', '%s')", tableID, paste0(tableID, ".csv")),
                                       bsicons::bs_icon("filetype-csv", size = "1.5em", title = "Download CSV"))
            )
            
            if (toggleExpand)
            {
                ret <- htmltools::tagAppendChildren(
                    ret,
                    # NOTE: button is hidden if toggleExpandDisableIfNoGrouping since the initial selection is assumed
                    # to be no grouping.
                    htmltools::tags$button(type = "button", class = "btn btn-secondary btn-sm", id = toggleExpandID,
                                           style = paste("display:", if (!toggleExpandDisableIfNoGrouping) "" else "none"),
                                           onClick = sprintf("Reactable.toggleAllRowsExpanded('%s')", tableID),
                                           bsicons::bs_icon("caret-right-fill", size = "1.5em", title = "Expand/Collapse"))
                )
            }
            
            ret <- htmltools::tagAppendChild(ret, tagsPostButtons)
            return(ret)
        })
    }
)
